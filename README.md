# 16S-Metatranscriptomic-Analysis
NextFlow Pipelines for Qiime2 &amp; Kraken2 metagenomic/metatranscriptomic analyses

## Workflow
This repository contains 2 alternative pipelines that use different softwares to perform a metagenomic or metatranscriptomic analyses. 
The input for both pipelines are the reads coming from the sample of interest (not aligned with the host genome), 
already processed (quality control, adapters removal, etc.). The details of each of the pipelines can be found in the README.md files 
in each of the corresponding folders.   

Both pipelines generate a matrix with the abundance of each microorganism detected in each of the analysis samples. 
This matrix can be used, alternatively, for differential abundance analysis using the R package metagenomeSeq.
For this last step, it will be required two files containg the condition of each sample and the contrast to perform.  

Next diagram represents the workflow of both pipelines in this repository, with main inputs and outputs. 
The box marked with '*' is not included in the Nextflow code.

![Pipeline workflow](https://user-images.githubusercontent.com/116556564/199960696-5276126f-213d-4656-b426-a6e4c4e09fb0.png)


## Dependencies  

* [Nextflow](https://www.nextflow.io/)

According to the pipelines and steps (defined in the parameters of each pipeline) to be performed, dependencies of this workflow will be different.  

For the Qiime pipeline:  
* [Conda](https://docs.conda.io/en/latest/)
* [Qiime2](https://qiime2.org/) (v >= 2021.8)  
  
  
For the Kraken/bracken pipeline:
* [Kraken](https://github.com/DerrickWood/kraken2) (v >= 2.1.2)
* [Bracken](https://github.com/jenniferlu717/Bracken) (v >= 2.7)  
  If generating Krona reports:
  * Krona (v >= 2.7)
  
  
If differential abundance analysis have to be performed:
* R
* Rstudio
* R packages:
    * [metagenomeSeq](https://www.bioconductor.org/packages/release/bioc/html/metagenomeSeq.html)
    * [ggfortify](https://cran.r-project.org/web/packages/ggfortify/index.html)
    * [sva](https://bioconductor.org/packages/release/bioc/html/sva.html)

## Choosing between kraken and Qiime pipelines

## References
1. Di Tommaso, P.; Chatzou, M.; Floden, E.W.; Barja, P.P.; Palumbo, E.; Notredame, C. Nextflow 
Enables Reproducible Computational Workflows. Nat. Biotechnol. 2017, 35, 316–319, https://doi.org/10.1038/nbt.3820.
2. Bolyen, E.; Rideout, J.R.; Dillon, M.R.; Bokulich, N.A.; Abnet, C.C.; Al-Ghalith, G.A.; Alexander, H.; 
Alm, E.J.; Arumugam, M.; Asnicar, F.; et al. Reproducible, Interactive, Scalable and Extensible Microbiome 
Data Science Using QIIME 2. Nat. Biotechnol. 2019, 37, 852–857, https://doi.org/10.1038/s41587-019-0209-9.
3. Wood, D.E.; Salzberg, S.L. Kraken: Ultrafast Metagenomic Sequence Classification Using Exact Alignments. 
Genome Biol. 2014, 15, R46, https://doi.org/10.1186/gb-2014-15-3-r46.
4. Lu, J.; Breitwieser, F.P.; Thielen, P.; Salzberg, S.L. Bracken: Estimating Species Abundance in Metagenomics 
Data. PeerJ Comput. Sci. 2017, 3, e104, https://doi.org/10.7717/peerj-cs.104.
5. Lee, C.; Lee, S.; Park, T. A Comparison Study of Statistical Methods for the Analysis Metagenome Data. 
In Proceedings of the 2017 IEEE International Conference on Bioinformatics and Biomedicine (BIBM); November 2017; pp. 1777–1781.
