include { FASTP                      } from '../modules/nf-core/fastp/main'
include { BWAMEM2_MEM                } from '../modules/nf-core/bwamem2/mem/main'
include { BWAMEM2_INDEX              } from '../modules/nf-core/bwamem2/index/main'
include { BBMAP_REFORMAT_STANDARDISE } from '../modules/local/bbmap/reformat_standardise/main'
include { samplesheetToList          } from 'plugin/nf-schema'

workflow READMAPPING {
    main:
    // Parse databases from parameters
    bwa_db_ch = channel
        .from(
            params.databases.collect { k, v ->
                if (v instanceof Map) {
                    if (v.containsKey('files')) {
                        return [id: k] + v
                    }
                }
            }.flatten()
        )
        .filter { it -> it }
    bwa_db_ch.view{ it -> "bwa_db_ch — ${it}" }

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
    samplesheet.view{ it -> "samplesheet — ${it}" }

    ch_versions = channel.empty()


    // De-interleave interleaved paired-end reads
    BBMAP_REFORMAT_STANDARDISE(
        samplesheet.map{ meta, reads, _fasta -> [meta, reads] }, 
        'fastq.gz'
    )
    reads = BBMAP_REFORMAT_STANDARDISE.out.reformated

    // QC
    if (params.skip_qc) {
        reads.set { qc_reads }
        qc_stats = channel.empty()
    }
    else {
        input_reads = reads.map{ meta, reads_ ->
            [meta, (reads_.size()==1) ? reads_[0] : reads_]
        }

        FASTP(
            input_reads,
            [],
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
    

    // Generate BWA-MEM2 indexes from FASTA files
    fasta_ch = samplesheet.map{ meta, _reads, fasta -> [meta, fasta] }
    fasta_ch.view{ it -> "fasta_ch — ${it}" }
    BWAMEM2_INDEX( fasta_ch )


    // Run mapping
    assembly_mapping_ch = BWAMEM2_INDEX.out.index
        .join(fasta_ch)
        .join(qc_reads)
        .map { meta, index, fasta, reads_ -> [meta + [db_id: 'assembly'], reads_, index, fasta] }
    assembly_mapping_ch.view{ it -> "assembly_mapping_ch — ${it}" }
        
    bwa_index_ch = bwa_db_ch.map{ meta -> [meta, [file(meta.files.index), file(meta.files.fasta)]] }
    db_mapping_ch = qc_reads
        .combine(bwa_index_ch)
        .map { reads_, db -> [reads_[0] + [db_id: db[0].id], reads_[1], db[1][0], db[1][1]] }
    db_mapping_ch.view{ it -> "db_mapping_ch — ${it}" }

    mapping_ch = assembly_mapping_ch.mix(db_mapping_ch)
        .multiMap{ meta, reads_, index, fasta ->
            reads: [meta, reads_]
            index: [[id: meta.db_id], index]
            fasta: [[id: meta.db_id], fasta]
        }
    BWAMEM2_MEM(mapping_ch.reads, mapping_ch.index, mapping_ch.fasta, false)


    emit:
    versions = ch_versions                 // channel: [ path(versions.yml) ]
}
