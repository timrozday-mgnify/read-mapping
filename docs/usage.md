# read-mapping: Usage

## Samplesheet input

Prepare a comma-separated samplesheet and pass it with `--input`:

```bash
--input '[path to samplesheet file]'
```

### Samplesheet format

```csv
sample,fastq_1,fastq_2,fasta,single_end
sample1,reads_R1.fastq.gz,reads_R2.fastq.gz,assembly.fasta,false
sample2,reads_interleaved.fastq.gz,,assembly.fasta,false
sample3,reads.fastq.gz,,assembly.fasta,true
```

| Column | Required | Description |
|---|---|---|
| `sample` | Yes | Sample identifier — no spaces. |
| `fastq_1` | Yes | Gzipped FASTQ (`.fastq.gz` or `.fq.gz`). For interleaved paired-end reads, supply only this column and set `single_end` to `false`. |
| `fastq_2` | No | Second FASTQ for paired-end reads. Leave blank for single-end or interleaved. |
| `fasta` | Yes | Sample assembly FASTA (`.fa`, `.fna`, `.fasta`, optionally `.gz`). Used for self-mapping alongside any configured databases. |
| `single_end` | Yes | `true` = single-end reads; `false` = paired-end or interleaved. |

An [example samplesheet](../assets/samplesheet.csv) is included with the pipeline.

## Configuring reference databases

Databases are configured in `nextflow.config` (or a `-params-file`). Each aligner reads from its own parameter block because BBMap and Minimap2 work directly from FASTA while BWA-MEM2 and Bowtie2 require pre-built indexes.

```nextflow
params {
    // Used by BBMap and Minimap2
    databases {
        human {
            files { fasta = '/path/to/human.fasta' }
        }
    }

    // Used by BWA-MEM2
    bwa_databases {
        human {
            files {
                index = '/path/to/bwa-mem2/human'   // index prefix
                fasta = '/path/to/human.fasta'
            }
        }
    }

    // Used by Bowtie2
    bowtie_databases {
        human {
            files {
                index = '/path/to/bowtie2/human'    // index prefix
                fasta = '/path/to/human.fasta'
            }
        }
    }
}
```

Each sample is also automatically mapped against its own assembly (the `fasta` column). Assembly results appear in the output with the database ID `{sample}_assembly`.

## Selecting aligners

| Parameter | Default | Aligner |
|---|---|---|
| `--use_bwamem2` | `true` | BWA-MEM2 |
| `--use_minimap2` | `false` | Minimap2 |
| `--use_bbmap` | `false` | BBMap |
| `--use_bowtie2` | `false` | Bowtie2 |

Multiple aligners can be enabled simultaneously:

```bash
nextflow run main.nf -profile docker \
  --input samplesheet.csv \
  --outdir results/ \
  --use_bwamem2 true \
  --use_minimap2 true
```

## Running the pipeline

```bash
nextflow run main.nf \
  -profile docker \
  --input samplesheet.csv \
  --outdir results/
```

### Params file

Repeated parameters can be collected in a YAML params file:

```bash
nextflow run main.nf -profile docker -params-file params.yaml
```

```yaml
input: './samplesheet.csv'
outdir: './results/'
use_bwamem2: true
use_minimap2: false
```

> [!WARNING]
> Do not use `-c <file>` to specify parameters. Custom config files passed with `-c` are for resource tuning only (process labels, executor settings, etc.).

### Skipping QC

Pass `--skip_qc true` to skip fastp trimming and run alignment on raw reads. Samples with zero reads after filtering are dropped automatically when QC is enabled.

## Core Nextflow arguments

### `-profile`

Available container/environment profiles: `docker`, `singularity`, `apptainer`, `podman`, `shifter`, `charliecloud`, `conda`, `mamba`.

Additional utility profiles: `test` (minimal test data), `debug` (verbose logging), `wave` (Seqera Wave containers), `arm64` (Apple Silicon / ARM HPC nodes), `gpu`.

Profiles can be combined: `-profile test,docker`. Order matters — later profiles override earlier ones.

### `-resume`

Resume a previous run from cached results:

```bash
nextflow run main.nf -resume [run-name]
```

Use `nextflow log` to list previous run names.

## Resource requirements

Default resource allocations per process:

| Process | CPUs | Memory |
|---|---|---|
| BBMAP_REFORMAT_STANDARDISE | 1 | 16 GB |
| FASTP | 4 | 4–8 GB |
| BWAMEM2_INDEX / BOWTIE2_BUILD | 4 | 32 GB |
| BWAMEM2_MEM | 4 | 32 GB |
| BBMAP_ALIGN / BOWTIE2_ALIGN / MINIMAP2_ALIGN | 4 | 64 GB |

To override resources, add a custom config:

```nextflow
// custom.config
process {
    withName: 'BWAMEM2_MEM' {
        memory = 64.GB
        cpus   = 8
    }
}
```

and pass it with `-c custom.config`.

See the [nf-core documentation](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) for full details on tuning resources.

## Running in the background

Use `-bg` to detach Nextflow from your terminal:

```bash
nextflow run main.nf -profile docker --input samplesheet.csv --outdir results/ -bg
```

Alternatively, run inside `screen` or `tmux`, or submit the Nextflow process itself to a cluster scheduler.

To limit Nextflow JVM memory:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```
