Cell-Cycle-on-C1
================

_Single-cell transcriptome analysis by fluorescence labeling_.

This repository contains source code and explanations on how to control quality
in transcriptome analysis of fluorescent cells in the [Fluidigm C1
platform](https://www.fluidigm.com/products/c1-system).  These explanations
cover in particular:

 - Fluorescence [measurement](fluorescence/Fluorescence-measured-in-ImageJ.md)
   and [normalisation](Intensity_correction/BackgroundCorrection.md).

 - Quantification and control of the [yield](cDNA_concentration/cDNA_concentration.md)
   of the single-cell cDNA amplification taking place in the machine.

 - Controls on [sequencing yield](HiSeq/HiSeq.md).

 - [Aggregation and consistency checks](combine_all/combined.md) of these data.

As supporting data for these scripts, we provide a RNA-seq analysis of _Fucci_
cells, available in [public repositories](DDBJ/DDBJ.md).  Deeper analysis of
these data is ongoing and will be the topic of a separate publication.

The [authors](AUTHORS.md) of the source code in this repository dedicate it to
the public domain under the Creative Commons Zero [License](LICENSE).


Required software
-----------------

As a system to produce and reproduce formatted reports, we use the
[`knitr`](http://yihui.name/knitr/) engine.

### `R` libraries:

`gdata`, `ggplot2`, `knitr`, `reshape`, `XML`, `moments`, `lattice`, `flexmix`, `limma`, `MASS`, `rms`, `contrast`

In Debian, some of these libraries are provided by the packages `r-cran-gdata`,
`r-cran-ggplot2`, `r-cran-reshape`

### Other software:

`lftp`, `Fiji ImageJ`

