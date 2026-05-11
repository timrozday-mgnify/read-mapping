# read-mapping: Output

## Output directory structure

All paths are relative to `--outdir`.

```
{outdir}/
├── {sample}/
│   ├── qc/
│   │   └── {sample}_qc.fastp.json
│   ├── bwamem2/           # if --use_bwamem2
│   │   ├── {sample}_{database_id}.bam
│   │   └── {sample}_{sample}_assembly.bam
│   ├── bowtie2/           # if --use_bowtie2
│   │   ├── {sample}_{database_id}.bam
│   │   └── {sample}_{sample}_assembly.bam
│   ├── minimap2/          # if --use_minimap2
│   │   ├── {sample}_{database_id}.bam
│   │   └── {sample}_{sample}_assembly.bam
│   └── bbmap/             # if --use_bbmap
│       ├── {sample}_{database_id}.bam
│       └── {sample}_{sample}_assembly.bam
└── pipeline_info/
    ├── execution_timeline_{timestamp}.html
    ├── execution_report_{timestamp}.html
    ├── execution_trace_{timestamp}.txt
    └── pipeline_dag_{timestamp}.html
```

The `{database_id}` in BAM file names corresponds to the key used in the database configuration block (e.g. `human`, `tomato`). Assembly self-mapping results use `{sample}_assembly` as the database ID.

## Per-step outputs

### Read standardisation (always runs)

BBMap `reformat.sh` de-interleaves interleaved paired-end reads, verifies pairing of explicitly paired reads, and trims read name descriptions. Output is not published — it feeds directly into the QC or alignment steps.

### Quality control

<details markdown="1">
<summary>Output files</summary>

- `{sample}/qc/{sample}_qc.fastp.json` — fastp summary statistics (read counts before and after filtering, quality metrics, adapter content).

</details>

[fastp](https://github.com/OpenGEN/fastp) trims adapters and low-quality bases. Samples with zero reads remaining after filtering are dropped and not aligned. Skip this step with `--skip_qc true`.

### Alignment

<details markdown="1">
<summary>Output files</summary>

- `{sample}/{aligner}/{sample}_{database_id}.bam` — sorted BAM for each sample × database combination.

</details>

Each enabled aligner produces one BAM per (sample, reference) pair. References include all configured databases plus the sample's own assembly. BAM files are not indexed by the pipeline — index with `samtools index` if needed downstream.

**BBMap** (`--use_bbmap`): aligns with `minid=0.95` and `maxindel=80` to reduce cross-mapping. Reference FASTA is provided at run time (no index pre-building).

**BWA-MEM2** (`--use_bwamem2`, default): assembly indexes are built at runtime; database indexes must be pre-built and configured under `params.bwa_databases`.

**Bowtie2** (`--use_bowtie2`): assembly indexes are built at runtime; database indexes must be pre-built and configured under `params.bowtie_databases`.

**Minimap2** (`--use_minimap2`): aligns directly to FASTA; configured under `params.databases`.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/execution_timeline_{timestamp}.html` — per-process timing.
- `pipeline_info/execution_report_{timestamp}.html` — resource usage and run summary.
- `pipeline_info/execution_trace_{timestamp}.txt` — tab-separated process-level trace.
- `pipeline_info/pipeline_dag_{timestamp}.html` — workflow DAG.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) generates these reports automatically for every run.
