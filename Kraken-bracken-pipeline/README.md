# Kraken-Bracken Pipeline

Metatranscriptomics (RNA) data analysis pipeline using kraken for OTU detection and bracken for correction of detected abundance values.

## Input

There are some important inputs that are requested and other inputs and options that have a default value, which can also be modified. 
You can check all of these parameters by using the --help command of the pipeline.  

All inputs and options can be modified either by command line or directly by changing their default value in the nextflow.config file. 
To modify them by command line you have to type "--name of the variable desired value" right after the execution command.  

### reads  
Location where your input fastq files are. Example:  

```--reads '/path/to/reads/*.{1,2}.fastq.gz' ```

The name of the path must be **enclosed in quotes** and it must contain at least **one '*'** wildcard character, in order to catch all samples required.  

It is possible to specify whether the reads are paired or not (parameter pairedEnd). 
If they are paired, it is necessary to use **{1,2}** notation to specify read pairs.  

### krakendb  
Path to Kraken database. This database should be downloaded and builded before executing this pipeline, however, 
there is no need for an specific database to perform.

Nevertheless, it is important to reserve enough memory to load the chosen database. The default memory to load the database is 5 GB, but it 
can be tunned by changing the parameter "kraken_mem".  

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

The main output obtained from this pipeline is a file with the extension ".tab" that contains the abundance of each OTU detected in each sample, 
and other file with the extension ".tsv" that contains the names of the different OTUs detected. Krona plots are also obtained (if skip_krona=false)
in a unique directory called "krona_reports"

Additionally, if differential abundance analisys is performed, a file with differentially abundant OTUs with some interest statistics (logFC, 
p-value, etc.) will be obtained, as different clustering plots of the different samples in the experiment(in a unique pdf file).


## Help
```
$ nextflow run kraken_braken.nf --help
N E X T F L O W  ~  version 22.04.5

 Name: nf-differential-abundance-kraken2
 Author: Bioinformatica IPPBLN <bioinformatica@ipb.csic.es>
=========================================
Mandatory arguments:
  --reads                       Path to input data (if paired end sequences must be a regular expression such as *{1,2}.fastq.gz)
  --krakendb                    Path to kraken database
Settings:
  --kraken_mem                  Necesary memory to load kraken database. Default = 5 GB
  --confidence                  Confidence score threshold (0-1). Default = 0
  --pairedEnd                   Specifies if reads are paired-end (true | false). Default = true
  --skip_bracken_build          Skip building bracken database (true | false). Default = true
  --skip_krona                  Skip generating krona reports (true | false). Default = true
    --krona_dir                   Path to krona directory. Default = /usr/bin/Krona-master/KronaTools/bin/
  --taxonomy_filter             Specifies the taxonomic level to filter by bracken. Defaults to S
  --kmer_len                    Specifies the kmer length. Default = 35
  --read_len                    Specifies the read length of the input data (needed for Bracken). Default = 75
  --b_threshold                 Specifies threshold for bracken filter. Default = 10
  --skip_diff                   Skip metagenome-Seq differential abundance analysis (true | false). Default = true
    --targets                     Metadata file that contains a tab-delimited table with filenames and contrast condition in the 3rd column.
    --contrast                    File with contrasts to perform in the abundance analysis.
Options:
  --outdir                      The output directory where the results will be saved. Defaults to ./
  --help  --h                   Shows this help page

Usage example:
    nextflow run main.nf --reads '/path/to/paired_end_reads_*.{1,2}.fastq.gz' --krakendb '/path/to/krakendb/' --krona_dir '/path/to/ktImportTaxonomy'         --targets '/path/to/targets.txt' --contrast '/path/to/contrast.txt'

```
