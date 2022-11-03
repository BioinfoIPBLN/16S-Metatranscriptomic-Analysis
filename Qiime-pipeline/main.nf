def helpMessage() {
    log.info"""
     Name: nf-differential-abundance-qiime
     Author: Bioinformatica IPBLN <bioinformatica@ipb.csic.es>
    =========================================
    Mandatory arguments:
      --data_dir                    Path to input data (must be surrounded with quotes)
      --metadata                    Path to Kraken2 Database
      --classifier                  Path to classifier file (.qza) used for feature classification (if skip_build_classifier = false)
    Settings:
      --casava                      Specifies if the input data is in casava format (true | false). Default to ${params.casava}
      --pairedEnd                   Specifies if reads are paired-end (true | false). Default = ${params.pairedEnd}
      --demultiplexed               Specified if the input data is demultiplexed (true | false). Default = ${params.demultiplexed}
        --brcds_col                   If not demultiplexed, name of the column with barcodes.
      --chimera_method              Method used by dada2 to remove chimeras ('none' | 'consensus' | 'pooled'). Default = ${params.chimera_method}
      --trunc_len                   Position to truncate reads at (if pairedEnd, the forward read). Default = ${params.trunc_len}
      --trunc_len_r                 Position to truncate the reverse reads at (when pairedEnd). Default = ${params.trunc_len_r}
      --trim                        Position to trim reads at (if pairedEnd, the forward read). Default = ${params.trim}
      --trim_r                      Position to trim the reverse reads at (when pairedEnd). Default = ${params.trim_r}
      --depth                       Frequency that each sample should be rarefied (1:Inf). Default = ${params.depth}
      --level                       Specifies the taxonomic level to filter by (1:15). Default = ${params.level}
      --beta_div                    Make a beta diversity analysis (true | false). Default = ${params.beta_div}
        --beta_col                    Name of the variable for the beta diversity analysis.
      --skip_build_classifier       Skip the classifier building. Default = ${params.skip_build_classifier}
        --taxa                        Path to input data.
        --levels                      Path to input data.
      --skip_diff                   Skip metagenome-Seq differential abundance analysis (true | false). Default = ${params.skip_diff}
        --targets                     Metadata file that contains a tab-delimited table with filenames and contrast condition in the 3rd column.
        --contrast                    File with contrasts to perform in the abundance analysis.
    Options:
      --outdir                      The output directory where the results will be saved. Defaults to ${params.outdir}
      --help                        Shows this help page

    Usage example:
      nextflow run bioinfo/qiime-nf --data_dir '/path/to/qiime_data/' --metadata '/path/to/metadata.tsv' --classifier '/path/to/classifier.qza'
    """.stripIndent()
}

if (params.help){
    helpMessage()
    exit 0
}

datadir_ch = Channel.fromPath(params.data_dir, checkIfExists: true )
    .ifEmpty { exit 1, "No directory matching: ${params.data_dir}\n" }
 
metadata_ch=Channel.fromPath(params.metadata, checkIfExists: true)
    .ifEmpty { exit 1, "No metadata matching: ${params.metadata}\n" }.first()


process qiime_imp {
  input:
    path data from datadir_ch

  output:
    path "*.qza" into qdata_ch

  script:
  if (params.casava) {
    if (params.pairedEnd) {
      type = "SampleData[PairedEndSequencesWithQuality]"
    } else {
      type = "SampleData[SequencesWithQuality]"
    }
  } else {
    if (params.pairedEnd) {
      if (params.demultiplexed){
        type = "MultiplexedPairedEndBarcodeInSequence"

      } else {
        type = "EMPPairedEndSequences"

      }
    } else {
      if (params.demultiplexed){
        type = "MultiplexedSingleEndBarcodeInSequence"

      } else {
        type = "EMPSingleEndSequences"

      }
    }
  }
    if  (params.casava){
      """
      qiime tools import \
      --type '${type}' \
      --input-path $data \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt \
      --output-path data.qza
      """

    } else{
      """
      qiime tools import \
      --type $type \
      --input-path $data \
      --output-path data.qza
      """
    }
}

if (!params.demultiplexed){
  process demultiplex {
    publishDir "${params.outdir}/Visualizations", mode: 'copy', pattern: '*.qzv'

    input:
      path seqs from qdata_ch
      path meta from metadata_ch

    output:
      path "demux.qza" into demux_ch
      path "*.qzv"

    script:
    if (params.pairedEnd) {
      """
      qiime demux emp-paired \
        --i-seqs $seqs \
        --m-barcodes-file $meta \
        --m-barcodes-column ${params.brcds_col} \
        --o-per-sample-sequences demux.qza \
        --o-error-correction-details demux-details.qza
      
      qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
      """
    } else {
      """
      qiime demux emp-single \
        --i-seqs $seqs \
        --m-barcodes-file $meta \
        --m-barcodes-column ${params.brcds_col} \
        --o-per-sample-sequences demux.qza \
        --o-error-correction-details demux-details.qza
      
      qiime demux summarize --i-data demux.qza --o-visualization demux.qzv
      """
    }
      
  }
} else{
  demux_ch = qdata_ch
}

process dada2 {
  publishDir "${params.outdir}/Visualizations", mode: 'copy', pattern: '*.qzv'

  input:
    path dseqs from demux_ch
    path meta from metadata_ch

  output:
    path "rep-seqs.qza" into repseqs_ch
    path "table.qza" into tab_ch
    path "*.qzv"
    env median_int into median_ch

  script:
  if (params.pairedEnd) {
    """
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs $dseqs \
      --p-trim-left-f ${params.trim} \
      --p-trim-left-r ${params.trim_r} \
      --p-trunc-len-f ${params.trunc_len} \
      --p-trunc-len-r ${params.trunc_len_r} \
      --p-chimera-method ${params.chimera_method} \
      --o-representative-sequences rep-seqs.qza \
      --o-table table.qza \
      --p-n-threads ${task.cpus} \
      --o-denoising-stats stats.qza
    
    qiime metadata tabulate --m-input-file stats.qza --o-visualization stats.qzv
    qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file $meta
    qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
    
    qiime tools export --input-path table.qzv --output-path ./
    
    median=\$(cat index.html | grep -A1 'Median' | grep td | head -1 | sed -r 's/[t><d/,]//g')
    median_int="\${median%.*}"
    """

  } else {
    """
    qiime dada2 denoise-single \
      --i-demultiplexed-seqs $dseqs \
      --p-trim-left ${params.trim} \
      --p-trunc-len ${params.trunc_len} \
      --p-chimera-method ${params.chimera_method} \
      --o-representative-sequences rep-seqs.qza \
      --o-table table.qza \
      --p-n-threads ${task.cpus} \
      --o-denoising-stats stats.qza

    qiime metadata tabulate --m-input-file stats.qza --o-visualization stats.qzv
    qiime feature-table summarize --i-table table.qza --o-visualization table.qzv --m-sample-metadata-file $meta
    qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
    
    qiime tools export --input-path table.qzv --output-path ./
    
    median=\$(cat index.html | grep -A1 'Median' | grep td | head -1 | sed -r 's/[t><d/,]//g')
    median_int="\${median%.*}"
    """
  }
}


tab_ch
    .into {tabdiv_ch; tabalpha_ch; tab_counts_ch; tab_barplot_ch} 
repseqs_ch
    .into {repseqs_phylo_ch; repseqs_class_ch} 

process phylogenetic_div {

  input:
    path rseqs from repseqs_phylo_ch
    path meta from metadata_ch

  output:
    path "aligned-rep-seqs.qza" into aliseqs_ch
    path "rooted-tree.qza" into tree_ch

  script:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences $rseqs \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza
    """
}

tree_ch
    .into {treediv_ch; treealpha_ch} 

process diversity {
  publishDir "${params.outdir}/Visualizations", mode: 'copy', pattern: '*.qzv'

  input:
    path tree from treediv_ch
    path table from tabdiv_ch
    path meta from metadata_ch

  output:
    path "core-metrics-results/faith*" into faith_ch
    path "core-metrics-results/evenness*" into evenness_ch
    path "core-metrics-results/unweighted_unifrac_distance_matrix*" into d_matrix_ch
    path "core-metrics-results/*.qzv"


  script:
    """
    qiime diversity core-metrics-phylogenetic \
    --i-phylogeny $tree \
    --i-table $table \
    --p-sampling-depth ${params.depth} \
    --m-metadata-file $meta \
    --output-dir core-metrics-results

    """
}

process alpha_diversity {
  publishDir "${params.outdir}/Visualizations", mode: 'copy', pattern: '*.qzv'

  input:
    path faith from faith_ch
    path meta from metadata_ch


  output:
    path "*.qzv" 

  script:
    """
    qiime diversity alpha-group-significance \
      --i-alpha-diversity $faith \
      --m-metadata-file $meta \
      --o-visualization faith-pd-group-significance.qzv
    """
}

process beta_diversity {
  publishDir "${params.outdir}/Visualizations", mode: 'copy', pattern: '*.qzv'

  input:
    path d_matrix from d_matrix_ch
    path meta from metadata_ch


  output:
    path "*.qzv" 
  
  when:
    params.beta_div

  script:
    """
    qiime diversity beta-group-significance \
      --i-distance-matrix $d_matrix \
      --m-metadata-file $meta \
      --m-metadata-column ${params.beta_col} \
      --o-visualization unweighted-unifrac-body-site-significance.qzv \
      --p-pairwise

    """
}

process alpha_rarefaction {
  publishDir "${params.outdir}/Visualizations", mode: 'copy', pattern: '*.qzv'
  
  input:
    path tree from treealpha_ch
    path table from tabalpha_ch
    path meta from metadata_ch
    val median from median_ch

  output:
    path "*.qzv" 

  script:
    """
    qiime diversity alpha-rarefaction \
    --i-table $table \
    --i-phylogeny $tree \
    --p-max-depth $median \
    --m-metadata-file $meta \
    --o-visualization alpha-rarefaction.qzv

    """
}
 
if (!params.skip_build_classifier){
  taxfasta_ch = Channel.fromPath(params.taxa, checkIfExists: true)
    .ifEmpty { exit 1, "No taxonomy.fna matching: ${params.taxa}\n" }
  tax_levels_ch = Channel.fromPath(params.levels, checkIfExists: true)
    .ifEmpty { exit 1, "No levels matching: ${params.levels}\n" }

  process build_classifier {
  
    input:
      path taxa from taxfasta_ch
      path levels from tax_levels_ch

    output:
      path "classifier.qza" into classifier_ch
    
    script:
      """
      qiime tools import --type 'FeatureData[Sequence]' \
        --input-path $taxa --output-path taxa.qza

      qiime tools import --type 'FeatureData[Taxonomy]' \
        --input-format HeaderlessTSVTaxonomyFormat \
        --input-path $levels \
        --output-path ref-taxonomy_all_levels.qza 
      
      qiime feature-classifier fit-classifier-naive-bayes \
        --i-reference-reads taxa.qza \
        --i-reference-taxonomy ref-taxonomy_all_levels.qza \
        --o-classifier classifier.qza --output-dir classifier-naive-bayes

      """
  }

} else {
  classifier_ch = Channel.fromPath(params.classifier, checkIfExists: true)
    .ifEmpty { exit 1, "No classifier matching: ${params.classifier}\n"}
}

process feat_classify {
  publishDir "${params.outdir}/Visualizations", mode: 'copy', pattern: '*.qzv'

  input:
    path classifier from classifier_ch
    path rseqs from repseqs_class_ch

  output:
    path "taxonomy.qza" into taxonomy_ch
    path "*.qzv"

  script:
    """
    mkdir -p $baseDir/temp
    export TMPDIR="$baseDir/temp/"
    qiime feature-classifier classify-sklearn \
      --i-classifier $classifier \
      --i-reads $rseqs \
      --o-classification taxonomy.qza

    qiime metadata tabulate \
      --m-input-file taxonomy.qza \
      --o-visualization taxonomy.qzv

    """
}

taxonomy_ch
    .into { taxonomy_barplot_ch; taxonomy_counts_ch}

process barplot {
  publishDir "${params.outdir}/Visualizations", mode: 'copy', pattern: '*.qzv'

  input:
    path taxa from taxonomy_barplot_ch
    path table from tab_barplot_ch
    path meta from metadata_ch

  output:
    path "*.qzv"

  script:
    """
    qiime taxa barplot \
      --i-table $table \
      --i-taxonomy $taxa \
      --m-metadata-file $meta \
      --o-visualization taxa-bar-plots.qzv
    """
}

process export_counts {
  publishDir "${params.outdir}", mode: 'copy', pattern: '*.tsv'

  input:
    path table from tab_counts_ch
    path mytaxa from taxonomy_counts_ch

  output:
    path "*.tsv" into counts_ch

  script:
    """
    qiime taxa collapse \
      --i-table $table \
      --i-taxonomy $mytaxa \
      --p-level ${params.level} \
      --o-collapsed-table table_export.qza

    qiime composition add-pseudocount \
      --i-table table_export.qza \
      --o-composition-table comp-table_export.qza

    qiime tools export --input-path comp-table_export.qza --output-path ./ 

    biom convert -i feature-table.biom -o feature_table.tsv --to-tsv
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
          path targets from targets_ch
          path contrast from contrast_ch
      
      output:
          path "*.{pdf,xls}"
      
      script:
          """
          metagenome_seq.R $counts $targets $contrast
          """
  }
}
