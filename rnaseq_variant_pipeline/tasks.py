"""
This pipeline largely based on GATK 3 guidelines for RNA-Seq variant calling
adapted for GATK 4.

Reference: https://software.broadinstitute.org/gatk/documentation/article.php?id=3891
"""

import datetime
from glob import glob
import logging
import os
from os.path import join

import luigi
from luigi.task import flatten_output, flatten, getpaths
from luigi.util import requires
from bioluigi.tasks import cutadapt, fastqc
from bioluigi.tasks.utils import RemoveTaskOutputOnFailureMixin, CreateTaskOutputDirectoriesBeforeRunMixin
from bioluigi.scheduled_external_program import ScheduledExternalProgramTask

logger = logging.getLogger('luigi-interface')

luigi.namespace('rnaseq_variant_pipeline')

class core(luigi.Config):
    genome = luigi.Parameter()
    annotations = luigi.Parameter()
    star_index_dir = luigi.Parameter()

    output_dir = luigi.Parameter(description='Output directory')
    filtered_dir = luigi.Parameter(description='Output directory for filtered variants')
    tmp_dir = luigi.Parameter()

cfg = core()

class ProduceGenome(luigi.Task):
    """
    Produce a reference genome sequence.
    """
    def output(self):
        return luigi.LocalTarget(cfg.genome)

@requires(ProduceGenome)
class ProduceGenomeDict(ScheduledExternalProgramTask):
    """
    Produce a sequence dictionnary to query efficiently a genome.
    """

    def program_args(self):
        return ['gatk', 'CreateSequenceDictionary', '-R', self.input().path, '-O', self.output().path]

    def output(self):
        return luigi.LocalTarget(os.path.splitext(self.input().path)[0] + '.dict')

class ProduceAnnotations(luigi.Task):
    """
    Produce a reference annotations.
    """
    def output(self):
        return luigi.LocalTarget(cfg.annotations)

@requires(ProduceGenome, ProduceAnnotations)
class ProduceReference(CreateTaskOutputDirectoriesBeforeRunMixin, ScheduledExternalProgramTask):
    """
    Produce a reference.
    TODO: generate the reference
    """

    cpus = 8
    memory = 32

    def program_args(self):
        return ['STAR',
                '--runThreadN', self.cpus,
                '--runMode', 'genomeGenerate',
                '--genomeDir', os.path.dirname(self.output().path),
                '--genomeFastaFiles', self.input()[0].path,
                '--sjdbGTFfile', self.input()[1].path]

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

@requires(ProduceSampleFastqs)
class TrimSample(luigi.Task):
    """
    Trim Illumina universal adapters from reads.
    """

    walltime = datetime.timedelta(days=2)

    def run(self):
        r1, r2 = self.input()
        r1_out, r2_out = self.output()
        r1_out.makedirs()
        r2_out.makedirs()
        yield cutadapt.TrimPairedReads(
                r1.path, r2.path,
                r1_out.path, r2_out.path,
                adapter_3prime='AGATCGGAAGAGC',
                minimum_length=25,
                cpus=8)

    def output(self):
        return [luigi.LocalTarget(join(cfg.output_dir, 'data-trimmed', self.experiment_id, self.sample_id, os.path.basename(i.path)))
                for i in self.input()]

class DynamicWrapperTask(luigi.Task):
    """
    Similar to WrapperTask but for dynamic dependencies yielded in the body of
    the run() method.
    """
    def output(self):
        tasks = []
        if all(req.complete() for req in flatten(self.requires())):
            try:
                tasks = list(self.run())
            except:
                logger.exception('%s failed at run() step; the exception will not be raised because Luigi is still building the graph.', repr(self))

        # FIXME: conserve task structure: the generator actually create an
        # implicit array level even if a single task is yielded.
        # For now, we just handle the special singleton case.
        if len(tasks) == 1:
            tasks = tasks[0]

        return getpaths(tasks)

    def complete(self):
        # ensure that all static dependencies are satisfied
        if not all(req.complete() for req in flatten(self.requires())):
            return False

        # check that all yielded tasks are completed
        return all(req.complete() for chunk in self.run()
                   for req in flatten(chunk))

@requires(TrimSample)
class QualityControlSample(DynamicWrapperTask):
    """Generate FastQC reports for each mates of a given sample."""
    def run(self):
        destdir = join(cfg.output_dir, 'data-qc', self.experiment_id, self.sample_id)
        os.makedirs(destdir, exist_ok=True)
        yield [fastqc.GenerateReport(fastq.path, destdir)
               for fastq in self.input()]

class QualityControlExperiment(luigi.WrapperTask):
    """Generate FastQC reports for each sample in a given experiment."""
    experiment_id = luigi.Parameter()
    def requires(self):
        return [QualityControlSample(self.experiment_id, os.path.basename(sample_id))
                for sample_id in glob(join(cfg.output_dir, 'data', self.experiment_id, '*'))]

@requires(ProduceReference, TrimSample, QualityControlSample)
class AlignSample(CreateTaskOutputDirectoriesBeforeRunMixin, ScheduledExternalProgramTask):
    """
    Align a sample on a reference.
    """

    cpus = 8
    memory = 64

    def program_args(self):
        genome, fastqs, _ = self.input()
        args = ['STAR',
                '--genomeDir', os.path.dirname(genome.path),
                '--outFileNamePrefix', os.path.dirname(self.output().path) + '/',
                '--outSAMtype', 'BAM', 'SortedByCoordinate',
                '--runThreadN', self.cpus,
                # FIXME: '--readStrand', 'Forward',
                '--readFilesCommand', 'zcat']

        args.append('--readFilesIn')
        args.extend(fastq.path for fastq in fastqs)

        return args

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'aligned', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.bam'))

@requires(ProduceGenome, AlignSample)
class PrepareSampleReference(CreateTaskOutputDirectoriesBeforeRunMixin, ScheduledExternalProgramTask):
    """
    Prepare a personalized reference for the sample.
    """

    cpus = 8
    memory = 64

    def program_args(self):
        return ['STAR',
                '--runMode', 'genomeGenerate',
                '--genomeDir', os.path.dirname(self.output().path),
                '--genomeFastaFiles', self.input()[0].path,
                '--sjdbFileChrStartEnd', join(os.path.dirname(self.input()[1].path), 'SJ.out.tab'),
                '--sjdbOverhang', 75,
                '--runThreadN', self.cpus]

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'aligned-reference', self.experiment_id, self.sample_id, 'Genome'))

@requires(PrepareSampleReference, TrimSample)
class AlignStep2Sample(CreateTaskOutputDirectoriesBeforeRunMixin, ScheduledExternalProgramTask):
    """
    Perform an alignment (2-step)
    """

    cpus = 8
    memory = 64

    def program_args(self):
        args = ['STAR',
                '--genomeDir', os.path.dirname(self.input()[0].path),
                '--outFileNamePrefix', os.path.dirname(self.output().path) + '/',
                '--outSAMtype', 'BAM', 'SortedByCoordinate',
                '--runThreadN', self.cpus,
                # FIXME: '--readStrand', 'Forward',
                '--readFilesCommand', 'zcat']

        args.append('--readFilesIn')
        args.extend(fastq.path for fastq in self.input()[1])

        return args

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'aligned-step2', self.experiment_id, self.sample_id, 'Aligned.sortedByCoord.out.bam'))

@requires(AlignStep2Sample)
class AddOrReplaceReadGroups(CreateTaskOutputDirectoriesBeforeRunMixin, RemoveTaskOutputOnFailureMixin, ScheduledExternalProgramTask):

    cpus = 1
    memory = 32

    def program_args(self):
        return ['gatk',
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

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'grouped', self.experiment_id, self.sample_id + '.bam'))

@requires(AddOrReplaceReadGroups)
class MarkDuplicates(CreateTaskOutputDirectoriesBeforeRunMixin, RemoveTaskOutputOnFailureMixin, ScheduledExternalProgramTask):

    cpus = 1
    memory = 32

    def program_args(self):
        return ['gatk',
                'MarkDuplicates',
                '--INPUT', self.input().path,
                '--OUTPUT', self.output().path,
                '--CREATE_INDEX',
                '--TMP_DIR', cfg.tmp_dir,
                '--VALIDATION_STRINGENCY', 'SILENT',
                '--METRICS_FILE', join(os.path.dirname(self.output().path), '{}.metrics'.format(self.sample_id))]

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'duplicates-marked', self.experiment_id, '{}.bam'.format(self.sample_id)))

@requires(ProduceGenome, ProduceGenomeDict, MarkDuplicates)
class SplitNCigarReads(CreateTaskOutputDirectoriesBeforeRunMixin, RemoveTaskOutputOnFailureMixin, ScheduledExternalProgramTask):

    cpus = 1
    memory = 32
    walltime = datetime.timedelta(days=10)

    def program_args(self):
        genome, genome_dict, duplicates = self.input()
        return ['gatk',
                'SplitNCigarReads',
                '-R', genome.path,
                '-I', duplicates.path,
                '-O', self.output().path,
                '--tmp-dir', cfg.tmp_dir]

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'ncigar-splitted', self.experiment_id, self.sample_id + '.bam'))

@requires(ProduceGenome, SplitNCigarReads)
class CallVariants(CreateTaskOutputDirectoriesBeforeRunMixin, RemoveTaskOutputOnFailureMixin, ScheduledExternalProgramTask):

    cpus = 8
    memory = 32
    walltime = datetime.timedelta(days=30)

    def program_args(self):
        genome, splitted = self.input()
        return ['gatk',
                'HaplotypeCaller',
                '--native-pair-hmm-threads', self.cpus,
                '-R', genome.path,
                '-I', splitted.path,
                '-O', self.output().path,
                '--dont-use-soft-clipped-bases',
                '--standard-min-confidence-threshold-for-calling', 30.0,
                '--tmp-dir', cfg.tmp_dir]

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, 'called', self.experiment_id, '{}.vcf.gz'.format(self.sample_id)))

@requires(ProduceGenome, CallVariants)
class FilterVariants(CreateTaskOutputDirectoriesBeforeRunMixin, RemoveTaskOutputOnFailureMixin, ScheduledExternalProgramTask):
    """Filter variants with GATK VariantFiltration tool."""

    cpus = 1
    memory = 32

    def program_args(self):
        genome, variants = self.input()
        return ['gatk',
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

    def output(self):
        return luigi.LocalTarget(join(cfg.output_dir, cfg.filtered_dir, self.experiment_id, '{}.vcf.gz'.format(self.sample_id)))

class FilterVariantsFromExperiment(luigi.WrapperTask):
    experiment_id = luigi.Parameter()
    def requires(self):
        return [FilterVariants(self.experiment_id, os.path.basename(sample_id))
                for sample_id in glob(join(cfg.output_dir, 'data', self.experiment_id, '*'))]
