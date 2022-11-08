# Qiime Pipeline

Metagenomics (DNA) data analysis pipeline using qiime2 for OTU detection and abundance quantification. This pipeline is mainly recommended for 16S genomic data.

## Input

There are some important inputs that are requested and other inputs and options that have a default value, which can also be modified. 
You can check all of these parameters by using the --help command of the pipeline.  

All inputs and options can be modified either by command line or directly by changing their default value in the nextflow.config file. 
To modify them by command line you have to type "--name of the variable desired value" right after the execution command.  

### data_dir  
The path to the directory containing your input files. This name must be **enclosed in quotes**, 
and the directory must only contain the input files that have to be analysed.

It is possible to specify the format of the input files by modifying other options. Depending on the format specified, more input data may be required (consult --help command). 
By default, this pipeline will expect paired reads with Cassava 1.8 format.   

### metadata  
Location of a file that contains all the metadata of the samples analysed. This file must meet qiime 2 requirements for this input 
file (consult [qiime2 documentation](https://docs.qiime2.org/2022.8/))   

### classifier  
Location of the classifier in a .qza format that will be used to make the taxonomic classification. If this file is not available, 
it can be generated through this pipeline, providing 2 additional files (taxa and levels). However, it is worth mentioning 
that this process may require considerable computational time. 

### targets

A tab-separated values file that contains the names and condition of every sample in our study. An example of the format of this file is the following:

| Filename |     Name     | Type | Covariable |
|:--------:|:------------:|:----:|:----------:|
| GSMXXXX1 | ConditionA 1 |   A  |    Male    |
| GSMXXXX2 | ConditionA 2 |   A  |   Female   |
| GSMXXXX3 | ConditionB 1 |   B  |   Female   |
  
It is required that the third column contains the condition that will be compared in the differential abundance analysis. 
It is also possible (but optional), to add a 4th column that correspond to a determined covariable to be fitted in the model (in the example, sex variable).

### contrast

A text file where every line correspond to a contrast to make in the differential abundance analysis. 
The format of writting down the contrast is the following one:  

**Name_of_contrast = (CaseCondition - ControlCondition)**  

It is worth to mention that the first line must contain no contrast and some word such as "Name" or "Contrasts".


## Output

The main output obtained from this pipeline is a tab-separated values file that contains the abundance of each OTU detected in each sample. It also 
generates many visualization plots that can be displayed using the online [qiime2 viewer](https://view.qiime2.org/)

Additionally, if differential abundance analisys is performed, a file with differentially abundant OTUs with some interest statistics (logFC, 
p-value, etc.) will be obtained, as different clustering plots of the different samples in the experiment(in a unique pdf file).


## Help
```
-bash-4.2$ nextflow run qiime.nf --help
NOTE: Nextflow is not tested with Java 1.8.0_242 -- It's recommended the use of version 11 up to 18

N E X T F L O W  ~  version 22.04.5
Launching `qiime.nf` [cheesy_aryabhata] DSL1 - revision: 7378317877

 Name: nf-differential-abundance-qiime
 Author: Bioinformatica IPBLN <bioinformatica@ipb.csic.es>
=========================================
Mandatory arguments:
  --data_dir                    Path to input data (must be surrounded with quotes)
  --metadata                    Path to Kraken2 Database
  --classifier                  Path to classifier file (.qza) used for feature classification (if skip_build_classifier = false)
Settings:
  --casava                      Specifies if the input data is in casava format (true | false). Default to true
  --pairedEnd                   Specifies if reads are paired-end (true | false). Default = true
  --demultiplexed               Specified if the input data is demultiplexed (true | false). Default = true
    --brcds_col                   If not demultiplexed, name of the column with barcodes.
  --chimera_method              Method used by dada2 to remove chimeras ('none' | 'consensus' | 'pooled'). Default = consensus
  --trunc_len                   Position to truncate reads at (if pairedEnd, the forward read). Default = 0
  --trunc_len_r                 Position to truncate the reverse reads at (when pairedEnd). Default = 0
  --trim                        Position to trim reads at (if pairedEnd, the forward read). Default = 0
  --trim_r                      Position to trim the reverse reads at (when pairedEnd). Default = 0
  --depth                       Frequency that each sample should be rarefied (1:Inf). Default = 1
  --level                       Specifies the taxonomic level to filter by (1:15). Default = 7
  --beta_div                    Make a beta diversity analysis (true | false). Default = false
    --beta_col                    Name of the variable for the beta diversity analysis.
  --skip_build_classifier       Skip the classifier building. Default = true
    --taxa                        Path to input data.
    --levels                      Path to input data.
  --skip_diff                   Skip metagenome-Seq differential abundance analysis (true | false). Default = true
    --targets                     Metadata file that contains a tab-delimited table with filenames and contrast condition in the 3rd column.
    --contrast                    File with contrasts to perform in the abundance analysis.
Options:
  --outdir                      The output directory where the results will be saved. Defaults to ./
  --help                        Shows this help page

Usage example:
  nextflow run bioinfo/qiime-nf --data_dir '/path/to/qiime_data/' --metadata '/path/to/metadata.tsv' --classifier '/path/to/classifier.qza'

```
