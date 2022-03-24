
# ArchRtoSignac : an Object Conversion Package from ArchR to Signac


ArchRtoSignac is an R package for performing object conversion from ArchRProject [(ArchR)](https://www.archrproject.com/index.html) to Signac SeuratObject [(Signac)](https://satijalab.org/signac/index.html).

ArchR and Signac are both commonly used scATAC-seq analysis packages with a comparable set of features and under development, which means they are likely to change over time. Researchers are welcome to stay only with Signac or ArchR pipeline during the whole analyzing process. But for the time being, Swarup lab would like to provide an option to help with data formatting from ArchRProject to Signac SeuratObject. ArchRtoSignac is a wrapper function utilizing existing ArchR and Seurat functions for an easier implementation.
The functions mentioned in this protocol, by all means, are only here for additional support if researchers want to use a fixed-width peak matrix for its advantage in computation, but work with Signac or want to access additional functions, packages, or resources than the ones provided in ArchR.


## Installation

We recommend creating an R [conda environment](https://docs.conda.io/en/latest)
environment for the scATAC-seq analysis with ArchR (v1), Signac (v1.5), Seurat (v4), since they provide several dependencies for ArchRtoSignac.

```bash
# create new conda environment for R
conda create -n scATAC -c conda-forge r-base r-essentials

# activate conda environment
conda activate scATAC

```

Next, open up R and install the required dependencies:

* [ArchR](https://www.archrproject.com/index.html), a general-purpose toolkit for single-cell ATAC sequencing analysis.
* [Seurat](https://satijalab.org/seurat/index.html), a general-purpose toolkit for single-cell RNA sequencing analysis.
* [Signac](https://satijalab.org/signac/index.html), a general-purpose toolkit for single-cell ATAC sequencing analysis.
* [devtools](https://devtools.r-lib.org/), a package for package development in R.
* [stringr](https://cran.r-project.org/web/packages/stringr/readme/README.html), a package for data cleaning and preparation in R.

```r
# install BiocManager
install.packages("BiocManager")

# install Bioconductor core packages
BiocManager::install()

# install additional packages:
install.packages(c("ArchR","Signac","Seurat","devtools","stringr"))

```

Now you can install the ArchRtoSignac package using `devtools`.

```r
devtools::install_github('swaruplabUCI/ArchRtoSignac')

```
