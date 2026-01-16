include { FASTP                      } from '../modules/nf-core/fastp/main'
include { BWAMEM2_MEM                } from '../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX              } from '../modules/nf-core/bwamem2/index/main'
include { MINIMAP2_ALIGN             } from '../modules/nf-core/minimap2/align/main'
include { BOWTIE2_BUILD              } from '../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN              } from '../modules/nf-core/bowtie2/align/main'
include { BBMAP_ALIGN                } from '../modules/nf-core/bbmap/align/main'
include { BBMAP_REFORMAT_STANDARDISE } from '../modules/local/bbmap/reformat_standardise/main'
include { samplesheetToList          } from 'plugin/nf-schema'

workflow READMAPPING {
    main:
    // Parse samplesheet and fetch reads
    samplesheet = channel.fromList(samplesheetToList(params.input, "${workflow.projectDir}/assets/schema_input.json"))
        .map { sample, fq1, fq2, fasta, single_end ->
            def single_file = (fq2 == [])
            [
                [
                    id: sample,
                    single_end: single_end,
                    interleaved: (!single_end) && single_file,
                ],
                single_file ? [file(fq1)] : [file(fq1), file(fq2)],
                file(fasta)
            ]
        }
    // samplesheet.view{ it -> "samplesheet — ${it}" }

    ch_versions = channel.empty()


    // De-interleave interleaved paired-end reads
    BBMAP_REFORMAT_STANDARDISE(
        samplesheet.map{ meta, reads, _fasta -> [meta, reads] }, 
        'fastq.gz'
    )
    reads = BBMAP_REFORMAT_STANDARDISE.out.reformated
        .map { it -> [it.head(), it.tail()] }

    // QC
    if (params.skip_qc) {
        reads.set { qc_reads }
        qc_stats = channel.empty()
    }
    else {
        input_reads = reads.map{ meta, reads_ ->
            [meta, (reads_.size()==1) ? reads_[0] : reads_, []]
        }

        FASTP(
            input_reads,
            false,
            false,
            false,
        )
        ch_versions = ch_versions.mix(FASTP.out.versions)

        qc_reads = FASTP.out.reads.map{ meta, reads_ ->
            [meta, (reads_ instanceof Collection) ? reads_ : [reads_]]
        }
        qc_stats = FASTP.out.json

        qc_read_counts = qc_stats.map {
            meta, json_file ->
            def json = new groovy.json.JsonSlurper().parseText(json_file.text)
            return tuple(
                meta,
                tuple(
                    json["summary"]["before_filtering"]["total_reads"],
                    json["summary"]["after_filtering"]["total_reads"],
                )
            )
        }

        qc_reads = qc_reads.join(qc_read_counts)
            .map{ meta, reads_, counts ->
                tuple(
                    meta + ['read_count': counts[0], 'qc_read_count': counts[1]],
                    reads_
                )
            }
            .filter{ meta, _reads -> meta.qc_read_count > 0 }
    }
    

    if (params.use_bbmap) {
        bbmap_assembly_mapping_ch = samplesheet.map{ meta, _reads, fasta -> [[id: meta.id], fasta] }
            .combine(qc_reads.map{ meta, reads_ -> [[id: meta.id], reads_] })
            .map { db_meta, fasta, reads_meta, reads_ -> [reads_meta + [db_id: "${db_meta.id}_assembly"], reads_, fasta] }
        bbmap_db_mapping_ch = qc_reads
            .combine(bwa_db_ch.map{ meta -> [meta, file(meta.files.fasta)] })
            .map { reads_meta, reads_, db_meta, db -> [reads_meta + [db_id: db_meta.id], reads_, db] }
        bbmap_mapping_ch = bbmap_assembly_mapping_ch.mix(bbmap_db_mapping_ch)
            .multiMap{ meta, reads_, fasta ->
                reads: [meta, reads_]
                fasta: fasta
            }
        BBMAP_ALIGN(bbmap_mapping_ch.reads, bbmap_mapping_ch.fasta)
    } 

    if (params.use_minimap2) {
        minimap_assembly_mapping_ch = samplesheet.map{ meta, _reads, fasta -> [[id: meta.id], fasta] }
            .combine(qc_reads.map{ meta, reads_ -> [[id: meta.id], reads_] })
            .map { db_meta, fasta, reads_meta, reads_ -> [reads_meta + [db_id: "${db_meta.id}_assembly"], reads_, fasta] }
        minimap_db_mapping_ch = qc_reads
            .combine(bwa_db_ch.map{ meta -> [meta, file(meta.files.fasta)] })
            .map { reads_meta, reads_, db_meta, db -> [reads_meta + [db_id: db_meta.id], reads_, db] }
        minimap_mapping_ch = minimap_assembly_mapping_ch.mix(minimap_db_mapping_ch)
            .multiMap{ meta, reads_, fasta ->
                reads: [meta, reads_]
                fasta: fasta
            }
        MINIMAP2_ALIGN(minimap_mapping_ch.reads, minimap_mapping_ch.fasta, true, false, false, false)
    } 

    if (params.use_bowtie2) {
        // Parse databases from parameters
        bowtie_db_ch = channel
            .from(
                params.bowtie_databases.collect { k, v ->
                    if (v instanceof Map) {
                        if (v.containsKey('files')) {
                            return [id: k] + v
                        }
                    }
                }.flatten()
            )
            .filter { it -> it }
        // bowtie_db_ch.view{ it -> "bowtie_db_ch — ${it}" }

        // Generate BWA-MEM2 indexes from FASTA files
        fasta_ch = samplesheet.map{ meta, _reads, fasta -> [meta, fasta] }
        // fasta_ch.view{ it -> "fasta_ch — ${it}" }
        BOWTIE2_BUILD(fasta_ch)

        // Run mapping
        bowtie_assembly_mapping_ch = BOWTIE2_BUILD.out.index.map{ meta, index -> [[id: meta.id], index] }
            .join(fasta_ch.map{ meta, fasta -> [[id: meta.id], fasta] })
            .combine(qc_reads.map{ meta, reads_ -> [[id: meta.id], reads_] })
            .map { db_meta, index, fasta, reads_meta, reads_ -> [reads_meta + [db_id: "${db_meta.id}_assembly"], reads_, index, fasta] }
        // bowtie_assembly_mapping_ch.view{ it -> "assembly_mapping_ch — ${it}" }
            
        bowtie_index_ch = bowtie_db_ch.map{ meta -> [meta, file(meta.files.index), file(meta.files.fasta)] }
        bowtie_db_mapping_ch = qc_reads
            .combine(bowtie_index_ch)
            .map { reads_meta, reads_, db_meta, index, fasta -> [reads_meta + [db_id: db_meta.id], reads_, index, fasta] }
        // bowtie_db_mapping_ch.view{ it -> "db_mapping_ch — ${it}" }

        bowtie_mapping_ch = bowtie_assembly_mapping_ch.mix(bowtie_db_mapping_ch)
            .multiMap{ meta, reads_, index, fasta ->
                reads: [meta, reads_]
                index: [[id: meta.db_id], index]
                fasta: [[id: meta.db_id], fasta]
            }
        BOWTIE2_ALIGN(bowtie_mapping_ch.reads, bowtie_mapping_ch.index, bowtie_mapping_ch.fasta, false, false)
    }

    if (params.use_bwamem2) {
        // Parse databases from parameters
        bwa_db_ch = channel
            .from(
                params.bwa_databases.collect { k, v ->
                    if (v instanceof Map) {
                        if (v.containsKey('files')) {
                            return [id: k] + v
                        }
                    }
                }.flatten()
            )
            .filter { it -> it }
        // bwa_db_ch.view{ it -> "bwa_db_ch — ${it}" }

        // Generate BWA-MEM2 indexes from FASTA files
        fasta_ch = samplesheet.map{ meta, _reads, fasta -> [meta, fasta] }
        // fasta_ch.view{ it -> "fasta_ch — ${it}" }
        BWAMEM2_INDEX(fasta_ch)

        // Run mapping
        assembly_mapping_ch = BWAMEM2_INDEX.out.index.map{ meta, index -> [[id: meta.id], index] }
            .join(fasta_ch.map{ meta, fasta -> [[id: meta.id], fasta] })
            .combine(qc_reads.map{ meta, reads_ -> [[id: meta.id], reads_] })
            .map { db_meta, index, fasta, reads_meta, reads_ -> [reads_meta + [db_id: "${db_meta.id}_assembly"], reads_, index, fasta] }
        // assembly_mapping_ch.view{ it -> "assembly_mapping_ch — ${it}" }
            
        bwa_index_ch = bwa_db_ch.map{ meta -> [meta, file(meta.files.index), file(meta.files.fasta)] }
        db_mapping_ch = qc_reads
            .combine(bwa_index_ch)
            .map { reads_meta, reads_, db_meta, index, fasta -> [reads_meta + [db_id: db_meta.id], reads_, index, fasta] }
        // db_mapping_ch.view{ it -> "db_mapping_ch — ${it}" }

        mapping_ch = assembly_mapping_ch.mix(db_mapping_ch)
            .multiMap{ meta, reads_, index, fasta ->
                reads: [meta, reads_]
                index: [[id: meta.db_id], index]
                fasta: [[id: meta.db_id], fasta]
            }
        BWAMEM2_MEM(mapping_ch.reads, mapping_ch.index, mapping_ch.fasta, false)
    }


    emit:
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}
