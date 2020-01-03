This README describes the scripts used for analyzing the result of a deep mutational scanning experiment that aims to identify mutations that allow A/Panama/2007/1999 (H3N2) to escape from antibody 2G12.

### REQUIREMENTS
* [Python](https://www.python.org/) version 2.7
* [R](https://www.r-project.org) version 3.6.1

### FILES
* All raw sequencing reads, which can be downloaded from NIH SRA database [PRJNA598865](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA598865), should be placed in fastq/ folder. The filename for read 1 should match those described in [./data/SampleID.tsv](./data/SampleID.tsv). The filename for read 2 should be the same as read 1 except "R1" is replaced by "R2".
* [./Fasta/Pan99_refseq.fa](./Fasta/Pan99_refseq.fa): Reference sequence (wild-type sequence) for Pan99 amplicon.
* [./doc/SampleID.tsv](./doc/SampleID.tsv): Describes the sample identity for each fastq file.
* [./doc/mut_list.tsv](./doc/mut_list.tsv): A list of the mutations of interests.

### ANALYSIS PIPELINE
1. [./Pan99_read_to_count.py](./Pan99_read_to_count.py): Converts raw reads to variant counts.
    - Input file:
      - ./fastq/\*_R1_\*.fastq
    - Output file:
      - [./result/Pan99_mut_count.tsv](./result/Pan99_mut_count.tsv)
2. [./script/Pan99_count_to_fit.py](./script/Pan99_count_to_fit.py): Computes relative fitness based on read counts.
    - Input file:
      - [./result/Pan99_mut_count.tsv](./result/Pan99_mut_count.tsv)
      - [./doc/mut_list.tsv](./doc/mut_list.tsv)
    - Output file:
      - [./result/Pan99_mut_fit.tsv](./result/Pan99_mut_fit.tsv)

### PLOTTING
1. [./script/Pan99_plot_fit_compare.R](./script/Pan99_plot_fit_compare.R): Plots correlation between replicates.
    - Input file:
      - [./result/Pan99_mut_fit.tsv](./result/Pan99_mut_fit.tsv)
    - Output file:
      - [./graph/compare_rep_Pan99_noAb.png](./graph/compare_rep_Pan99_noAb.png)
      - [./graph/compare_rep_Pan99_2G12.png](./graph/compare_rep_Pan99_2G12.png)
2. [./script/Pan99_plot_escape.R](./script/Pan99_plot_escape.R): Compare the fitness profiles with and without 2G12 selection
    - Input file:
      - [./result/Pan99_mut_fit.tsv](./result/Pan99_mut_fit.tsv)
    - Output file:
      - [./graph/compare_condition_Pan99.png](./graph/compare_condition_Pan99.png)
      - [./result/Pan99_esp.tsv](./result/Pan99_esp.tsv)
      - [./result/Pan99_fit.tsv](./result/Pan99_fit.tsv)
