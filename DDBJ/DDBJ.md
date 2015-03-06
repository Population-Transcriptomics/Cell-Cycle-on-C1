DDBJ data management
====================

Please install lftp if you are not already using it.
First create a directory in which to download all the 944 compressed fastq files.


```sh
mkdir 'yourdirectoryname'; cd 'yourdirectoryname'
lftp ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA002/DRA002399
mget */*.bz2
```

The following Script creates a table that has 6 columns: "Library",  "Well", "DDBJ", "Barcode", "ExperimentAccession", "cell_id" and 472 rows. Each row represents a single cell. 

- **"Library"**: the Flowcell run ID (e.g. "RNhi10371")
- **"Well"**: the coordinate of a cell on the 96-well plate after transfer of cDNA from the C1 IFC 
- **"DDBJ"**: the unique cell ID given to each cell by DDBJ
- **"Barcode"**: the left and right read (2 x 8 bases) barcode used to pool 96 cells per run
- **"ExperimentAccession"**: the name given by DDBJ that is the equivalent of the "Library" run ID 
- **"Run"**: the C1 IFC ID 
- **"cell_id"**: a unique identifier for each single cell, that is created by combining the C1 IFC ID and the "Well" with an underscore delimeter


```r
library(XML)
```

```
## Warning: package 'XML' was built under R version 3.1.2
```

```r
library(reshape2)
```

```
## Warning: package 'reshape2' was built under R version 3.1.2
```

```r
DDBJ <- xmlTreeParse("ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRA002/DRA002399/DRA002399.experiment.xml", useInternal=T)
top <- xmlRoot(DDBJ)
RunID <- xpathApply(top, "//LIBRARY_NAME", xmlValue)
# get DDBJ Run accession number
Experiment  <- xpathApply(top, "//EXPERIMENT", xmlAttrs)
Accession <- lapply(1:5, function(x) Experiment[[x]][1])
Barcodes <- xpathApply(top, "//LIBRARY_CONSTRUCTION_PROTOCOL", xmlValue)
filter1 <- gsub("\n", "", as.character(Barcodes))
filter2 <- strsplit(filter1, split=',')
names(filter2) <- paste(RunID, Accession, sep="_")
# load table that links well number with barcodes
WellID <- read.table("WellToBarcodes.txt", header=T, stringsAsFactors=F)

# make list of string vectors
filter3 <- lapply(1:5, function(x) as.vector(filter2[x]))
# create combined table in long format with run ID and DDBJ ID and barcode
long1 <- c()
y <- data.frame()
x <- c()
for (i in filter3) {
long1 <- melt(i)
y <- rbind(y, long1)
x <- x+1
}
link <- y; names(link) <- c("DDBJBarcode", "Library")
# remove white spaces
link <- as.data.frame(apply(link, 2, function(x) gsub('\\s+', '',x)))
# separate DDBJ name and barcode
intermediate <- as.data.frame(do.call(rbind, sapply(1:length(link$DDBJBarcode), function(x) strsplit(as.character(link$DDBJBarcode)[x], split=":"))))
final <- cbind(link, intermediate)
final <- final[,-1]; names(final) <- c("Library", "DDBJ", "Barcode")
# associate each barcode from the DDBJ XML with the well coordinate
final$Well <- WellID$Well[match(as.character(final$Barcode), as.character(WellID$Barcode))]
# split RunID and ExperiemntAccession
intermediate2 <- as.data.frame(do.call(rbind, sapply(1:length(final$Library), function(x) strsplit(as.character(final$Library)[x], split="_")))); names(intermediate2) <- c("Library", "ExperimentAccession")
final <- cbind(final[,-1], intermediate2)
# load combined metadata table to link to cell IDs
LibraryToC1ID <- read.table("LibraryToC1.txt", header=T, stringsAsFactors=F)
# get cell IDs for DDBJ names
final <- merge(final, LibraryToC1ID, by=c("Library", "Well"))
final$cell_id <- paste(final$Run, final$Well, sep="_")
head(final)
```

```
##     Library Well      DDBJ           Barcode ExperimentAccession
## 1 RNhi10371  A01 DRR028211 TAAGGCGA-TAGATCGC           DRX019711
## 2 RNhi10371  A02 DRR028171 CGTACTAG-TAGATCGC           DRX019711
## 3 RNhi10371  A03 DRR028147 AGGCAGAA-TAGATCGC           DRX019711
## 4 RNhi10371  A04 DRR028227 TCCTGAGC-TAGATCGC           DRX019711
## 5 RNhi10371  A05 DRR028195 GGACTCCT-TAGATCGC           DRX019711
## 6 RNhi10371  A06 DRR028219 TAGGCATG-TAGATCGC           DRX019711
##            Run          cell_id
## 1 1772-062-248 1772-062-248_A01
## 2 1772-062-248 1772-062-248_A02
## 3 1772-062-248 1772-062-248_A03
## 4 1772-062-248 1772-062-248_A04
## 5 1772-062-248 1772-062-248_A05
## 6 1772-062-248 1772-062-248_A06
```

```r
# write final output table to current working directory
write.table(final, "DDBJLink.txt", sep="\t", row.names=T, col.names=T, quote=F)
```

The next part utilises the file "DDBJLink.txt" to rename downloaded fastq.bz2 files with the unique "cell_id".
Adjust the below script to the directory in which the "DDBJLink.txt" is located. Furthermore your working directory should be 'yourdirectoryname' where all the dowloaded fastq.bz2 files are saved. 


```sh
# the code below replaces the "DDBJ" name of each fastq.bz2 fiel pair with the corresponding "cell_id"
cut -f4,8 DDBJLink.txt | sed 1d | while read from to ; do mv ${from}_1.fastq.bz2 ${to}.1.fastq.bz2 ; mv ${from}_2.fastq.bz2 ${to}.2.fastq.bz2; done
```
