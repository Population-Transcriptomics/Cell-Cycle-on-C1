

DDBJ data management
====================

The following description is only suitable for the download and renaming of fastq file stored at [DDBJ](http://trace.ddbj.nig.ac.jp/DRASearch/submission?acc=DRA002399). However, the same data can be obtained from [NCBI](http://www.ncbi.nlm.nih.gov/Traces/sra/?study=DRP002435) and [EBI](https://www.ebi.ac.uk/ena/data/view/DRP002435) 

First create a directory in which to download all the 945 compressed fastq files.


```sh
# create a directory for data download, go to that directory and start the download of all fastq files from the 
mkdir 'yourdirectoryname'; cd 'yourdirectoryname'
lftp ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA002/DRA002399
mget -c */*.bz2
```

The following script creates a table that has 6 columns: "Library",  "Well", "DDBJ", "Index", "ExperimentAccession", "cell_id", and 472 rows. Each row represents a single cell. 

- **"Library"**: the sequencing flowcell run ID (e.g. "RNhi10371")
- **"Well"**: the coordinate of a cell on the 96-well plate after transfer of cDNA from the C1 capture array 
- **"DDBJ"**: the unique cell ID given to each cell by DDBJ
- **"Index"**: the left and right read (2 x 8 bases) barcode used to pool 96 cells per run
- **"ExperimentAccession"**: the name given by DDBJ that is the equivalent of the "Library" run ID 
- **"Run"**: the C1 capture array ID 
- **"cell_id"**: a unique identifier for each single cell, that is created by combining the C1 IFC ID and the "Well" with an underscore delimeter


```r
# load required packages
library(XML)
```

```
## Warning: package 'XML' was built under R version 3.1.2
```

```r
library(reshape)
# parse the xml tree structure into R
DDBJ <- xmlTreeParse("ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA002/DRA002399/DRA002399.experiment.xml", useInternal=T)
# make the nodes of the tree readable with common R commands
top <- xmlRoot(DDBJ)
# finds and returns values of a specific node as a list
RunID <- xpathApply(top, "//LIBRARY_NAME", xmlValue)
# get DDBJ Run accession number
Experiment  <- xpathApply(top, "//EXPERIMENT", xmlAttrs)
Accession <- lapply(1:5, function(x) Experiment[[x]][1])
# get the cell barcodes
Barcodes <- xpathApply(top, "//LIBRARY_CONSTRUCTION_PROTOCOL", xmlValue)
# 
filter <- gsub("\n", "", as.character(Barcodes))
filter <- strsplit(filter, split=',')
names(filter) <- paste(RunID, Accession, sep="_")
# load table that links well number with barcodes
# the barcode sequences were taken from the Fluidigm protocol 'PN 100-5950 A1' 
WellID <- read.csv("WellToBarcodes.csv", header=T, stringsAsFactors=F)
# make list of string vectors
filter <- lapply(1:5, function(x) as.vector(filter[x]))
# create combined table in long format with run ID and DDBJ ID and barcode
link <- melt(filter)[,1:2]; names(link) <- c("DDBJIndex", "Library")
# remove white spaces
link <- as.data.frame(apply(link, 2, function(x) gsub('\\s+', '',x)))
# separate DDBJ name and Index
intermediate <- read.table(text=as.character(link$DDBJIndex), sep=':')
final <- data.frame(link[,-1], intermediate); names(final) <- c("Library", "DDBJ", "Index")
# associate each Index from the DDBJ XML with the Well coordinates
final$Well <- WellID$Well[match(as.character(final$Index), as.character(WellID$Index))]
# split RunID and ExperimentAccession in column 1 of final
intermediate <- read.table(text=as.character(final$Library), sep='_'); names(intermediate) <- c("Library", "ExperimentAccession")
final <- cbind(final[,-1], intermediate)
# link Library with Run
LibraryToC1ID <- c(RNhi10371="1772-062-248", RNhi10372="1772-062-249", RNhi10395="1772-064-103", RNhi10396="1772-067-038", RNhi10397="1772-067-039")
# get cell IDs for DDBJ names
final$Run <- as.character(sapply(final$Library, function(x) LibraryToC1ID[x]))
final$cell_id <- paste(final$Run, final$Well, sep="_")
head(final)
```

```
##        DDBJ             Index Well   Library ExperimentAccession          Run          cell_id
## 1 DRR028133 AAGAGGCA-AAGGAGTA  G11 RNhi10371           DRX019711 1772-062-248 1772-062-248_G11
## 2 DRR028134 AAGAGGCA-ACTGCATA  F11 RNhi10371           DRX019711 1772-062-248 1772-062-248_F11
## 3 DRR028135 AAGAGGCA-AGAGTAGA  D11 RNhi10371           DRX019711 1772-062-248 1772-062-248_D11
## 4 DRR028136 AAGAGGCA-CTAAGCCT  H11 RNhi10371           DRX019711 1772-062-248 1772-062-248_H11
## 5 DRR028137 AAGAGGCA-CTCTCTAT  B11 RNhi10371           DRX019711 1772-062-248 1772-062-248_B11
## 6 DRR028138 AAGAGGCA-GTAAGGAG  E11 RNhi10371           DRX019711 1772-062-248 1772-062-248_E11
```


```r
# write final output table to current working directory
write.csv(final, "DDBJLink.csv", row.names=F)
```

The next part utilises the file "DDBJLink.csv" to rename downloaded fastq.bz2 files with the unique "cell_id".
Adjust the below script to the directory in which the "DDBJLink.csv" is located. Furthermore, your working directory should be 'yourdirectoryname' where all the dowloaded fastq.bz2 files are saved. 


```sh
# the code below replaces the "DDBJ" name of each fastq.bz2 file pair with the corresponding "cell_id"
cut -f4,8 DDBJLink.csv | sed 1d | while read from to ; do mv ${from}_1.fastq.bz2 ${to}.1.fastq.bz2 ; mv ${from}_2.fastq.bz2 ${to}.2.fastq.bz2; done
```

The file fastq_md5sum.csv contains the md5 checksums for all fastq.bz2 files after renaming and can be used for validity checks by comparing the md5sums of your local renamed files.


```sh
# execute in the directory with the renamed files
ls -lh | awk '{print $9}' | xargs md5sum >> ../fastq_md5sums.txt
```
