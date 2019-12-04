from glob import glob
import os
from os.path import join

import luigi
from luigi.contrib.external_program import ExternalProgramTask
from luigi.util import requires

"""
This pipeline largely based on GATK 3 guidelines for RNA-Seq variant calling
adapted for GATK 4.

Reference: https://software.broadinstitute.org/gatk/documentation/article.php?id=3891
"""

class rnaseq_variant_pipeline(luigi.Config):
    star_bin = luigi.Parameter()
    gatk_bin = luigi.Parameter()

    genome = luigi.Parameter()
    annotations = luigi.Parameter()
    star_index_dir = luigi.Parameter()

    output_dir = luigi.Parameter()
    tmp_dir = luigi.Parameter()

cfg = rnaseq_variant_pipeline()

class ProduceGenome(luigi.Task):
    """
    Produce a reference genome sequence.
    """
    def output(self):
        return luigi.LocalTarget(cfg.genome)

@requires(ProduceGenome)
class ProduceGenomeDict(ExternalProgramTask):
    """
    Produce a sequence dictionnary to query efficiently a genome.
    """

    def program_args(self):
        return [cfg.gatk_bin, 'CreateSequenceDictionary', '-R', self.input().path, '-O', self.output().path]

    def output(self):
        return luigi.LocalTarget(os.path.splitext(self.input().path)[0] + '.dict')

class ProduceAnnotations(luigi.Task):
    """
    Produce a reference annotations.
    """
    def output(self):
        return luigi.LocalTarget(cfg.annotations)

@requires(ProduceGenome, ProduceAnnotations)
class ProduceReference(ExternalProgramTask):
    """
    Produce a reference.
    TODO: generate the reference
    """

    resources = {'cpu': 8, 'mem': 32}

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
        return luigi.LocalTarget(join(cfg.star_index_dir, 'Genome'))

class ProduceSampleFastqs(luigi.Task):
    """
    Produce the FASTQs that relate to a sample.
    """
    experiment_id = luigi.Parameter()
    sample_id = luigi.Parameter()
    def output(self):
        return [luigi.LocalTarget(f) for f in sorted(glob(join(cfg.output_dir, 'data', self.experiment_id, self.sample_id, '*.fastq.gz')))]

@requires(ProduceReference, ProduceSampleFastqs)
class AlignSample(ExternalProgramTask):
    """
    Align a sample on a reference.
    """

    resources = {'cpu': 8, 'mem': 32}

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

@requires(ProduceGenome, AlignSample)
class PrepareSampleReference(ExternalProgramTask):
    """
    Prepare a personalized reference for the sample.
    """

    resources = {'cpu': 8, 'mem': 32}

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

@requires(PrepareSampleReference, ProduceSampleFastqs)
class AlignStep2Sample(ExternalProgramTask):
    """
    Perform an alignment (2-step)
    """

    resources = {'cpu': 8, 'mem': 32}

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

@requires(AlignStep2Sample)
class AddOrReplaceReadGroups(ExternalProgramTask):

    resources = {'cpu': 1, 'mem': 24}

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

@requires(AddOrReplaceReadGroups)
class MarkDuplicates(ExternalProgramTask):

    resources = {'cpu': 1, 'mem': 24}

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

@requires(ProduceGenome, ProduceGenomeDict, MarkDuplicates)
class SplitNCigarReads(ExternalProgramTask):

    resources = {'cpu': 1, 'mem': 24}

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

@requires(ProduceGenome, SplitNCigarReads)
class CallVariants(ExternalProgramTask):

    resources = {'cpu': 8, 'mem': 24}

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

@requires(ProduceGenome, CallVariants)
class FilterVariants(ExternalProgramTask):
    """Filter variants with GATK VariantFiltration tool."""

    resources = {'cpu': 1, 'mem': 24}

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