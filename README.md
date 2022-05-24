
# ArchRtoSignac : an Object Conversion Package for ArchR to Signac


**ArchRtoSignac** is an R package to convert an ArchRProject [(ArchR)](https://www.archrproject.com/index.html) to a Signac SeuratObject [(Signac)](https://satijalab.org/signac/index.html).

ArchR and Signac are both commonly used scATAC-seq analysis packages with comparable sets of features and are currently under development, which means they are likely to change over time. You can choose to use only one of these packages; however, you may want to use both packages for your analysis. For example, we use ArchR to generate a fixed-width peak matrix due to its computational advantage, and we use Signac for reference mapping to assist in cell-type annotation. Here we provide an option to help with the data formatting from an ArchRProject to a Signac SeuratObject: **ArchRtoSignac**, a wrapper function that allows easier implementation of both pipelines. In addition, conversion to a SeuratObject allows the use of other packages available through SeuratWrappers.

---

## Installation

We recommend creating an R [conda environment](https://docs.conda.io/en/latest) specifically for scATAC-seq analysis to install the required packages. This ensures that software versions required here do not conflict with software required for other projects, and several dependencies for ArchRtoSignac will be automatically installed.

```bash
# create new conda environment for R
conda create -n scATAC -c conda-forge r-base r-essentials

# activate conda environment
conda activate scATAC

```

Next, open up R and install ArchRtoSignac using `devtools`.

```r
# install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# install ArchRtoSignac
devtools::install_github("swaruplabUCI/ArchRtoSignac")

```

When installing ArchRtoSignac, the following required dependencies should be automatically installed.

* [ArchR](https://www.archrproject.com/index.html), a general-purpose toolkit for single-cell ATAC sequencing analysis.
* [Seurat](https://satijalab.org/seurat/index.html), a general-purpose toolkit for single-cell RNA sequencing analysis.
* [Signac](https://satijalab.org/signac/index.html), a general-purpose toolkit for single-cell ATAC sequencing analysis.
* [devtools](https://devtools.r-lib.org/), a package for package development in R.
* [biovizBase](https://www.bioconductor.org/packages/release/bioc/html/biovizBase.html), basic graphic utilities for visualization of genomic data in R.
* [stringr](https://cran.r-project.org/web/packages/stringr/readme/README.html), a package for data cleaning and preparation in R.

However, if there are issues with installation, please try the following:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# install Bioconductor core packages
BiocManager::install()

# install additional packages including ArchR, Signac Seurat and etc:
if (!requireNamespace("biovizBase", quietly = TRUE)) BiocManager::install("biovizBase")
if (!requireNamespace("ArchR", quietly = TRUE)) devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
if (!requireNamespace("Signac", quietly = TRUE)) install.packages("Signac")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")

```

---

## Usage

Check all required dependencies have been installed and load them automatically.
```r
packages <- c("ArchR","Seurat", "Signac","stringr") # required packages
loadinglibrary(packages)

```

Obtain ArchRProject peak matrix for object conversion.

```r
pkm <- getPeakMatrix(proj) # proj is an ArchRProject

```

Extract appropriate Ensembl gene annotation and convert to UCSC style.

```r
library(EnsDb.Hsapiens.v86) # Ensembl database to convert to human hg38. Install what is appropriate for your analysis

annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38") # "UCSC" is the default style to change to but can be changed with argument seqStyle

```

Convert ArchRProject to Signac SeuratObject.

```r
fragments_dir <- "path_to_cellranger_atac_output" # the directory before "/outs/" for all samples

seurat_atac <- ArchR2Signac(
  ArchRProject = proj,
  refversion = "hg38",
  #samples = samplelist, # list of samples in the ArchRProject (default will use ArchRProject@cellColData$Sample but another list can be provided)
  fragments_dir = fragments_dir,
  pm = pkm, # peak matrix from getPeakMatrix()
  #output_dir = "/outs/", # folder with "fragments.tsv.gz" ("/outs/" is the default)
  annotation = annotations # annotation from getAnnotation()
)

```

Transfer ArchRProject gene score matrix to Signac SeuratObject.

```r
gsm <- getGeneScoreMatrix(ArchRProject = proj, SeuratObject = seurat_atac)

seurat_atac[['RNA']] <- CreateAssayObject(counts = gsm)

```

Transfer ArchRProject dimension reduction ("IterativeLSI" or "Harmony") and UMAP to Signac SeuratObject.

```r
seurat_atac <- addDimRed(ArchRProject = proj, 
			 SeuratObject = seurat_atac, 
			 reducedDims = "IterativeLSI") # default is "IterativeLSI"
			 #add both 'Harmony' and ‘IterativeLSI’:
			 #reducedDims = c('IterativeLSI', 'Harmony')


```
