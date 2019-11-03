import luigi
from luigi.contrib.external_program import ExternalProgramTask
from glob import glob
import os
from os.path import join

"""
This pipeline largely based on GATK 3 guidelines for RNA-Seq variant calling
adapted for GATK 4.

Reference: https://software.broadinstitute.org/gatk/documentation/article.php?id=3891
"""

class rnaseq_variant_pipeline(luigi.Config):
    star_bin = luigi.Parameter()
    gatk_bin = luigi.Parameter()

    output_dir = luigi.Parameter()
    tmp_dir = luigi.Parameter()

cfg = rnaseq_variant_pipeline()

class ProduceGenome(luigi.Task):
    """
    Produce a reference genome sequence.
    """
    def output(self):
        return luigi.LocalTarget('genomes/hg19_ensembl98/Homo_sapiens.GRCh37.dna.primary_assembly.fa')

class ProduceGenomeDict(ExternalProgramTask):
    def requires(self):
        return ProduceGenome()

    def program_args(self):
        return [cfg.gatk_bin, 'CreateSequenceDictionary', '-R', self.input().path, '-O', self.output().path]

    def output(self):
        return luigi.LocalTarget('genomes/hg19_ensembl98/Homo_sapiens.GRCh37.dna.primary_assembly.dict')

class ProduceAnnotations(luigi.Task):
    """
    Produce a reference annotations.
    """
    def output(self):
        return luigi.LocalTarget('genomes/hg19_ensembl98/Homo_sapiens.GRCh37.87.gtf')

class ProduceReference(ExternalProgramTask):
    """
    Produce a reference.
    TODO: generate the reference
    """

    resources = {'cpu': 8, 'mem': 32}

    def requires(self):
        return [ProduceGenome(), ProduceAnnotations()]

    def program_args(self):
        return [cfg.star_bin,
                '--runThreadN', self.resources['cpu'],
                '--runMode', 'genomeGenerate',
                '--genomeDir', os.path.dirname(self.output().path),
                '--genomeFastaFiles', self.input()[0].path,
                '--sjdbGTFfile', self.input()[1].path]

    def run(self):
        self.output().makedirs()
        return super(ProduceReference, self).run()

    def output(self):
        return luigi.LocalTarget('references/hg19_ensembl98/Genome')

class ProduceSampleFastqs(luigi.Task):
    """
    Produce the FASTQs that relate to a sample.
    """
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()
    def output(self):
        return [luigi.LocalTarget(f) for f in sorted(glob(join(cfg.output_dir, 'data', self.experiment_id, self.sample_id, '*.fastq.gz')))]

class AlignSample(ExternalProgramTask):
    """
    Align a sample on a reference.
    """
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()

    resources = {'cpu': 8, 'mem': 32}

    def requires(self):
        return [ProduceReference(), ProduceSampleFastqs(self.experiment_id, self.sample_id)]

    def program_args(self):
        args = [cfg.star_bin,
                '--genomeDir', os.path.dirname(self.input()[0].path),
                '--outFileNamePrefix', os.path.dirname(self.output().path) + '/',
                '--outSAMtype', 'BAM', 'SortedByCoordinate',
                '--runThreadN', self.resources['cpu'],
                # FIXME: '--readStrand', 'Forward',
                '--readFilesCommand', 'zcat']

        args.append('--readFilesIn')
        args.extend(fastq.path for fastq in self.input()[1])

        return args

    def run(self):
        self.output().makedirs()
        return super(AlignSample, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'aligned', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.bam'))

class PrepareSampleReference(ExternalProgramTask):
    """
    Prepare a personalized reference for the sample.
    """
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()

    resources = {'cpu': 8, 'mem': 32}

    def requires(self):
        return [ProduceGenome(), AlignSample(self.experiment_id, self.sample_id)]

    def program_args(self):
        return [cfg.star_bin,
                '--runMode', 'genomeGenerate',
                '--genomeDir', os.path.dirname(self.output().path),
                '--genomeFastaFiles', self.input()[0].path,
                '--sjdbFileChrStartEnd', join(os.path.dirname(self.input()[1].path), 'SJ.out.tab'),
                '--sjdbOverhang', 75,
                '--runThreadN', self.resources['cpu']]

    def run(self):
        self.output().makedirs()
        return super(PrepareSampleReference, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'aligned-reference', self.experiment_id, self.sample_id, 'Genome'))

class AlignStep2Sample(ExternalProgramTask):
    """
    Perform an alignment (2-step)
    """
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()

    resources = {'cpu': 8, 'mem': 32}

    def requires(self):
        return [PrepareSampleReference(self.experiment_id, self.sample_id), ProduceSampleFastqs(self.experiment_id, self.sample_id)]

    def program_args(self):
        args = [cfg.star_bin,
                '--genomeDir', os.path.dirname(self.input()[0].path),
                '--outFileNamePrefix', os.path.dirname(self.output().path) + '/',
                '--outSAMtype', 'BAM', 'SortedByCoordinate',
                '--runThreadN', self.resources['cpu'],
                # FIXME: '--readStrand', 'Forward',
                '--readFilesCommand', 'zcat']

        args.append('--readFilesIn')
        args.extend(fastq.path for fastq in self.input()[1])

        return args

    def run(self):
        self.output().makedirs()
        return super(AlignStep2Sample, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'aligned-step2', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.bam'))

class AddOrReplaceReadGroups(ExternalProgramTask):
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()

    resources = {'cpu': 1, 'mem': 24}

    def requires(self):
        return AlignStep2Sample(self.experiment_id, self.sample_id)

    def program_args(self):
        return [cfg.gatk_bin,
                'AddOrReplaceReadGroups',
                '-I', self.input().path,
                '-O', self.output().path,
                '-SO', 'coordinate',
                '--RGID', 'id',
                '--RGLB', 'library',
                '--RGPL', 'platform',
                '--RGPU', 'machine',
                '--RGSM', 'sample',
                '--TMP_DIR', cfg.tmp_dir]

    def run(self):
        self.output().makedirs()
        return super(AddOrReplaceReadGroups, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'grouped', self.experiment_id, self.sample_id + '.bam'))

class MarkDuplicates(ExternalProgramTask):
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()

    resources = {'cpu': 1, 'mem': 24}

    def requires(self):
        return AddOrReplaceReadGroups(self.experiment_id, self.sample_id)

    def program_args(self):
        return [cfg.gatk_bin,
                'MarkDuplicates',
                '--INPUT', self.input().path,
                '--OUTPUT', self.output().path,
                '--CREATE_INDEX',
                '--TMP_DIR', cfg.tmp_dir,
                '--VALIDATION_STRINGENCY', 'SILENT',
                '--METRICS_FILE', join(os.path.dirname(self.output().path), '{}.metrics'.format(self.sample_id))]

    def run(self):
        self.output().makedirs()
        return super(MarkDuplicates, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'duplicates-marked', self.experiment_id, '{}.bam'.format(self.sample_id)))

class SplitNCigarReads(ExternalProgramTask):
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()

    resources = {'cpu': 1, 'mem': 24}

    def requires(self):
        return [ProduceGenome(), ProduceGenomeDict(), MarkDuplicates(self.experiment_id, self.sample_id)]

    def program_args(self):
        genome, genome_dict, duplicates = self.input()
        return [cfg.gatk_bin,
                'SplitNCigarReads',
                '-R', genome.path,
                '-I', duplicates.path,
                '-O', self.output().path,
                '--tmp-dir', cfg.tmp_dir]

    def run(self):
        self.output().makedirs()
        return super(SplitNCigarReads, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'ncigar-splitted', self.experiment_id, self.sample_id + '.bam'))

class CallVariants(ExternalProgramTask):
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()

    resources = {'cpu': 8, 'mem': 24}

    def requires(self):
        return [ProduceGenome(), SplitNCigarReads(self.experiment_id, self.sample_id)]

    def program_args(self):
        genome, splitted = self.input()
        return [cfg.gatk_bin,
                'HaplotypeCaller',
                '--native-pair-hmm-threads', self.resources['cpu'],
                '-R', genome.path,
                '-I', splitted.path,
                '-O', self.output().path,
                '--dont-use-soft-clipped-bases',
                '--standard-min-confidence-threshold-for-calling', 30.0,
                '--tmp-dir', cfg.tmp_dir]

    def run(self):
        self.output().makedirs()
        return super(CallVariants, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'called', self.experiment_id, '{}.vcf.gz'.format(self.sample_id)))

class FilterVariants(ExternalProgramTask):
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()

    resources = {'cpu': 1, 'mem': 24}

    def requires(self):
        return [ProduceGenome(), CallVariants(self.experiment_id, self.sample_id)]

    def program_args(self):
        genome, variants = self.input()
        return [cfg.gatk_bin,
                'VariantFiltration',
                '-R', genome.path,
                '-V', variants.path,
                '-O', self.output().path,
                '--cluster-window-size', 35,
                '--cluster-size', 3,
                '--filter-name', 'FS',
                '--filter-expression', 'FS > 30.0',
                '--filter-name', 'QD',
                '--filter-expression', 'QD < 2.0',
                '--tmp-dir', cfg.tmp_dir]

    def run(self):
        self.output().makedirs()
        return super(FilterVariants, self).run()

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'filtered', self.experiment_id, '{}.vcf.gz'.format(self.sample_id)))

class FilterVariantsFromExperiment(luigi.WrapperTask):
    experiment_id = luigi.Parameter()
    def requires(self):
        return [FilterVariants(self.experiment_id, os.path.basename(sample_id))
                for sample_id in glob(join(cfg.output_dir, 'data', self.experiment_id, '*'))]
