# read-mapping

[![GitHub Actions CI Status](https://github.com/timrozday/read-mapping/actions/workflows/nf-test.yml/badge.svg)](https://github.com/timrozday/read-mapping/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/timrozday/read-mapping/actions/workflows/linting.yml/badge.svg)](https://github.com/timrozday/read-mapping/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**read-mapping** is a Nextflow DSL2 pipeline for mapping sequencing reads to reference sequences (genomes or databases) using multiple alignment tools. Each sample can be mapped against external reference databases and against its own assembly in the same run. Reads are standardised and optionally quality-filtered before alignment, and any combination of aligners can be enabled per run.

**Pipeline steps:**

1. Read standardisation — de-interleave, verify pairing, trim descriptions ([BBMap reformat.sh](https://sourceforge.net/projects/bbmap/))
2. Quality control (optional) — adapter trimming and quality filtering ([fastp](https://github.com/OpenGEN/fastp))
3. Alignment (one or more):
   - [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2) _(default)_ — builds index from assembly at runtime; uses pre-built index for databases
   - [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) — builds index from assembly at runtime; uses pre-built index for databases
   - [Minimap2](https://github.com/lh3/minimap2) — aligns directly to FASTA
   - [BBMap](https://sourceforge.net/projects/bbmap/) — aligns directly to FASTA (minid=0.95)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet CSV:

```csv
sample,fastq_1,fastq_2,fasta,single_end
sample1,reads_R1.fastq.gz,reads_R2.fastq.gz,assembly.fasta,false
sample2,reads_interleaved.fastq.gz,,assembly.fasta,false
sample3,reads.fastq.gz,,assembly.fasta,true
```

| Column | Description |
|---|---|
| `sample` | Sample name (no spaces) |
| `fastq_1` | Path to FASTQ file (gzipped). For interleaved paired-end reads, provide only this column with `single_end` set to `false`. |
| `fastq_2` | Path to second FASTQ file for paired-end reads (optional). |
| `fasta` | Path to the sample's own assembly FASTA for self-mapping. |
| `single_end` | `true` for single-end reads, `false` for paired-end or interleaved. |

Run the pipeline:

```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheet.csv \
  --outdir results/
```

> [!WARNING]
> Provide pipeline parameters via the CLI or a `-params-file` YAML/JSON file. Do not use `-c` to pass parameters — custom config files passed with `-c` are for resource tuning only.

Configure reference databases in `nextflow.config` (or a params file). Each aligner has its own database block:

```nextflow
params {
    // Used by BBMap and Minimap2 (FASTA only)
    databases {
        mydb {
            files { fasta = '/path/to/ref.fasta' }
        }
    }

    // Used by BWA-MEM2 (pre-built index required)
    bwa_databases {
        mydb {
            files {
                index = '/path/to/bwa-mem2/index'
                fasta = '/path/to/ref.fasta'
            }
        }
    }

    // Used by Bowtie2 (pre-built index required)
    bowtie_databases {
        mydb {
            files {
                index = '/path/to/bowtie2/index'
                fasta = '/path/to/ref.fasta'
            }
        }
    }
}
```

Enable aligners:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results/ \
  --use_bwamem2 true \
  --use_minimap2 true \
  --use_bbmap false \
  --use_bowtie2 false
```

For more details, see [docs/usage.md](docs/usage.md).

## Pipeline output

For each sample, BAM files are written to `{outdir}/{sample}/{aligner}/`, named `{sample}_{database_id}.bam`. Assembly self-mapping results use the database ID `{sample}_assembly`.

See [docs/output.md](docs/output.md) for full details.

## Credits

read-mapping was written by Tim Rozday.

## Citations

An extensive list of references for the tools used by the pipeline can be found in [`CITATIONS.md`](CITATIONS.md).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
