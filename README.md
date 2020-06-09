# RNA-Seq Variant Pipeline

Pipeline for calling and analyzing variants from RNA-Seq data

## Usage

We provide a simpe `luigi-wrapper` script to set the `$PYTHONPATH` and `--module` flag for you.

```bash
./luigi-wrapper <task> <task_args>
```

## Prerequisites

Minimally configure `luigi.cfg` as described in the [example configuration](example.luigi.cfg) we provide.

You only need to provide a path to a FASTA file containing a primary assembly and a GTF for gene annotations
which is in turn used by STAR for read alignment.

Setup the Conda environment:

```
conda env create -f environment.yml
conda activate rnaseq-variant-pipeline
python setup.py install
```

Add your FASTQs under `pipeline-output/data/<experiment_id>/<sample_id>/` and run a task described in the
following section.

## Tasks

| Task                       | Description |
|--------------------------- |-------------|
| AlignSample                | First aligment with genomic reference |
| PrepareSampleReference     | Prepare a personalized genomic reference for the sample |
| AlignStep2Sample           | Second alignment with sample reference |
| AddOrReplaceReadGroups     | Replace read groups |
| MarkDuplicates             | Mark duplicated alignments |
| SplitNCigarReads           | Split N-CIGAR reads |
| CallVariants               | Call variants |
| FilterVariants             | Filter variants |
| FilterVariantsInExperiment | Meta-task to run all the sample of an experiment |
