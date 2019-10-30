# RNA-Seq Variant Pipeline

Pipeline for calling and analyzing variants from RNA-Seq data

## Usage

We provide a simpe `luigi-wrapper` script to set the `$PYTHONPATH` and `--module` flag for you.

```bash
./luigi-wrapper <task> <task_args>
```

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
