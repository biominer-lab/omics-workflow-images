## RNA-Seq Workflow

### Softwares

| Name              | Version   |
| ----------------- | --------- |
| sentieon-genomics | 202112.04 |
| strelka           | 2.9.10    |
| pindel            | 0.2.5b9   |
| qualimap          | 2.2.2d    |
| fastq-screen      | 0.15.2    |
| fastqc            | 0.11.9    |
| fastp             | 0.23.2    |

### Installation and Testing

```
mamba create -n rnaseq-workflow hisat2=2.2.1 samtools=1.14 bioconductor-ballgown=2.26.0 bioconductor-genefilter=1.76.0 qualimap=2.2.2d fastq-screen=0.15.2 fastqc=0.11.9 fastp=0.23.2 stringtie=2.2.1

conda activate rnaseq-workflow
```