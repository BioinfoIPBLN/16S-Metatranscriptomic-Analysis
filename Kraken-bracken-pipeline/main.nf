def helpMessage() {
    log.info"""
     Name: nf-differential-abundance-kraken2
     Author: Bioinformatica IPPBLN <bioinformatica@ipb.csic.es>
    =========================================
    Mandatory arguments:
      --reads                       Path to input data (if paired end sequences must be a regular expression such as *{1,2}.fastq.gz)
      --krakendb                    Path to kraken database
    Settings:
      --kraken_mem                  Necesary memory to load kraken database. Default = ${params.kraken_mem}
      --confidence                  Confidence score threshold (0-1). Default = ${params.confidence}
      --pairedEnd                   Specifies if reads are paired-end (true | false). Default = ${params.pairedEnd}
      --skip_bracken_build          Skip building bracken database (true | false). Default = ${params.skip_bracken_build}
      --skip_krona                  Skip generating krona reports (true | false). Default = ${params.skip_krona}
        --krona_dir                   Path to krona directory. Default = ${params.krona_dir}
      --taxonomy_filter             Specifies the taxonomic level to filter by bracken. Defaults to ${params.taxonomy_filter}
      --kmer_len                    Specifies the kmer length. Default = ${params.kmer_len}
      --read_len                    Specifies the read length of the input data (needed for Bracken). Default = ${params.read_len}
      --b_threshold                 Specifies threshold for bracken filter. Default = ${params.b_threshold}
      --skip_diff                   Skip metagenome-Seq differential abundance analysis (true | false). Default = ${params.skip_diff}
        --targets                     Metadata file that contains a tab-delimited table with filenames and contrast condition in the 3rd column.
        --contrast                    File with contrasts to perform in the abundance analysis.
    Options:
      --outdir                      The output directory where the results will be saved. Defaults to ${params.outdir}
      --help  --h                   Shows this help page
    
    Usage example:
        nextflow run bioinfo/krakenbracken-nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz' \
        --krakendb '/path/to/krakendb/' --krona_dir '/path/to/ktImportTaxonomy' \
        --targets '/path/to/targets.txt' --contrast '/path/to/contrast.txt'
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

read_ch = Channel.fromFilePairs(params.reads, checkIfExists: true, size: params.pairedEnd ? 2 : 1 )
    .ifEmpty { exit 1, "No reads matching: ${params.reads}\n" }
 
krakendb_ch=Channel.fromPath(params.krakendb, checkIfExists: true)
    .ifEmpty { exit 1, "No KrakenDB matching: ${params.krakendb}\n" }.first()


process kraken2 {
    label 'parallel'
    tag "kraken2 on $sample_id"

    input:
        path db from krakendb_ch
        tuple val(sample_id), path(reads) from read_ch

    output:
        tuple val(sample_id), file('*_kr.out') into kraken_out
        tuple val(sample_id), path('*.report.txt') into kr_report_ch

    script:
    out= sample_id+"_kr.out"
    repo_name= sample_id+".report.txt"
    if (params.pairedEnd) {
        """
        kraken2 \
--use-names --report-zero-counts --db $db --threads ${task.cpus} \
--output $out --gzip-compressed --paired ${reads[0]} ${reads[1]} \
--report $repo_name --confidence ${params.confidence}
	 """
    } else{
        """
        kraken2 \
--use-names --report-zero-counts --db $db --threads ${task.cpus} \
--output $out --report $repo_name --gzip-compressed ${reads[0]} \
--confidence ${params.confidence}
        """
    }
}


kr_report_ch
    .into {kraken2krona_ch; kraken2bracken_ch} 

process krona {
    publishDir "${params.outdir}/krona_reports", mode: 'copy', pattern: '*.html'

    input:
        tuple val(sample_id), path(report) from kraken2krona_ch

    output:
        path "*.html"

    when:
        !params.skip_krona
    script:
    """
    "${params.krona_dir}ktImportTaxonomy" -m 3 -t 5 $report -o "${sample_id}.html"
    """
}
if (!params.skip_bracken_build){
    process bracken_build {
        label 'building'

        input:
            path db from krakendb_ch
        
        output:
            env db_name into brackendb_ch

        script:
        """
        db_name=\$(readlink $db)
        bracken-build -d $db -t ${task.cpus} -k ${params.kmer_len} -l ${params.read_len}
        """
    }
} else {
    brackendb_ch = krakendb_ch
}

process bracken {
    tag "bracken on $sample_id"
    input:
        tuple val(sample_id), path(report) from kraken2bracken_ch
        path db from brackendb_ch

    output:
        path "*_bracken.tsv" into bracken_report

    script:
    """
    bracken -d $db -i $report -o "${sample_id}_bracken.tsv" \
-r ${params.read_len} -l ${params.taxonomy_filter} -t ${params.b_threshold}
    """
}

process counts_parse {
    publishDir params.outdir, mode: 'copy'
    input:
        path abundances from bracken_report.collect()

    output:
        path "*.tab" into counts_ch
        path "*_OTU.tsv" into otus_ch

    script:
    """
    counts_parser_nxt.pl -dir "${abundances}" -label "Experiment"
    """
}

if (!params.skip_diff){ 
    targets_ch = Channel.fromPath(params.targets, checkIfExists: true)
        .ifEmpty { exit 1, "No targets file matching: ${params.targets}\n" }

    contrast_ch = Channel.fromPath(params.contrast, checkIfExists: true)
        .ifEmpty { exit 1, "No contrast file matching: ${params.contrast}\n" }

    process metagenome_seq {

        publishDir "${params.outdir}/Metagenome_seq", mode: 'copy'
        input:
            path counts from counts_ch
            path otus from otus_ch
            path targets from targets_ch
            path contrast from contrast_ch
        
        output:
            path "*.{pdf,xls}"
        
        script:
            """
            Rscript '$baseDir/bin/metagenome_seq.R' $counts $otus $targets $contrast
            """
    }
}
