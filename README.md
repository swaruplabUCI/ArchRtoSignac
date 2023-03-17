
# ArchRtoSignac : an Object Conversion Package for ArchR to Signac

[![DOI](https://zenodo.org/badge/473458154.svg)](https://zenodo.org/badge/latestdoi/473458154)

**ArchRtoSignac** is an R package to convert an ArchRProject [(ArchR)](https://www.archrproject.com/index.html) to a Signac SeuratObject [(Signac)](https://satijalab.org/signac/index.html).

ArchR and Signac are both commonly used scATAC-seq analysis packages with comparable sets of features and are currently under development, which means they are likely to change over time. You can choose to use only one of these packages; however, you may want to use both packages for your analysis. For example, we use ArchR to generate a fixed-width peak matrix due to its computational advantage, and we use Signac for reference mapping to assist in cell-type annotation. Here we provide an option to help with the data formatting from an ArchRProject to a Signac SeuratObject: **ArchRtoSignac**, a wrapper function that allows easier implementation of both pipelines. In addition, conversion to a SeuratObject allows the use of other packages available through SeuratWrappers.

---
## How to cite

Shi, Zechuan; Das, Sudeshna; Morabito, Samuel; Miyoshi, Emily; Swarup, Vivek. (2022). Protocol for single-nucleus ATAC sequencing and bioinformatic analysis in frozen human brain tissue, STAR Protocols, Volume 3, Issue 3, DOI: https://doi.org/10.1016/j.xpro.2022.101491.

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
# load ArchRtoSignac
library(ArchRtoSignac)

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

* **STEP 0 - Check all required dependencies have been installed and load them automatically.**

```r
packages <- c("ArchR","Seurat", "Signac","stringr") # required packages
loadinglibrary(packages)

```

* **STEP 1 - Obtain ArchRProject peak matrix for object conversion.**

```r
pkm <- getPeakMatrix(proj) # proj is an ArchRProject

```

* **STEP 2 - Extract appropriate Ensembl gene annotation and convert to UCSC style.**

```r
library(EnsDb.Hsapiens.v86) # Ensembl database to convert to human hg38. Install what is appropriate for your analysis

annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38") # "UCSC" is the default style to change to but can be changed with argument seqStyle

```

* **STEP 3 - Convert ArchRProject to Signac SeuratObject.**

Option1: Fragments Files using for `fragments_fromcellranger` from 10X Genomics Cellranger ATAC output

Please select Yes for `fragments_fromcellranger`. Example `fragments_fromcellranger = "Yes"`

```r
# Option 1a: Set one directory containing the cellranger output for each sample
fragments_dir <- "path_to_cellranger_atac_output" # the directory before "/outs/" for all samples

seurat_atac <- ArchR2Signac(
  ArchRProject = proj,
  refversion = "hg38",
  #samples = samplelist, # list of samples in the ArchRProject (default will use ArchRProject@cellColData$Sample but another list can be provided)
  fragments_dir = fragments_dir,
  pm = pkm, # peak matrix from getPeakMatrix()
  fragments_fromcellranger = "Yes", # fragments_fromcellranger This is an Yes or No selection ("NO" | "N" | "No" or "YES" | "Y" | "Yes")
  fragments_file_extension = NULL, # Default - NULL: File_Extension for fragments files (typically they should be '.tsv.gz' or '.fragments.tsv.gz')
  annotation = annotations # annotation from getAnnotation()
)

# Option 1b: Set a list of directories containing the cellranger output for each sample 
# (this newly added code to take in a list of fragments' path work both for fragments from cellranger and fragments not from cellranger, and when fragments are not from cellranger, please provide fragments_file_extension)
#
# Also PLEASE MAKE SURE the order of the fragment_dirs for samples have the same order as samplelist 
# or the order of list from ArchRProject@cellColData$Sample 
fragments_dirs <- list(
  "/path/to/sample1/cellranger/output",
  "/path/to/sample2/cellranger/output",
  "/path/to/sample3/cellranger/output"
)

# # Optional: use when fragments_fromcellranger = "NO", please set the file extension for the fragments file
# fragments_file_extension <- ".fragments.tsv.gz"

# Call the ArchR2Signac function with the provided arguments
SeuratObject <- ArchR2Signac(
  ArchRProject = ArchRProject,
  refversion = refversion,
  samples = samples,
  fragments_dir = fragments_dirs,
  pm = pm,
  fragments_fromcellranger = "YES",
  annotation = annotation
)

```
Option2: Fragments Files using for `fragments_fromcellranger` from **NON** Cellranger ATAC output, ie: SnapATAC tools

Please select No for `fragments_fromcellranger`. Example `fragments_fromcellranger = "NO"`, Also remember to provide the `fragments_file_extension`, for example `fragments_fromcellranger = '.tsv.gz'` or `fragments_fromcellranger = '.fragments.tsv.gz'`.

```r
fragments_dir <- "/ArchR/HemeFragments/" # please see the fragments format provided by ArchR examples

```
Above is the directory accessing the fragments files.

For eample, Fragments files in the folder HemeFragments, which we can check them in terminal

```r
tree /ArchR/HemeFragments/

/ArchR/HemeFragments/
├── scATAC_BMMC_R1.fragments.tsv.gz
├── scATAC_BMMC_R1.fragments.tsv.gz.tbi
├── scATAC_CD34_BMMC_R1.fragments.tsv.gz
├── scATAC_CD34_BMMC_R1.fragments.tsv.gz.tbi
├── scATAC_PBMC_R1.fragments.tsv.gz
└── scATAC_PBMC_R1.fragments.tsv.gz.tbi

```
** **Possible issue** due to the fragments format if fragments files are not from cellranger actac out:
Reported in [#Issue3](https://github.com/swaruplabUCI/ArchRtoSignac/issues/3)
Please check out: Signac snATAC-seq fragment file [Format](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/output/fragments)

** **Solution**
https://github.com/stuart-lab/signac/issues/748

Now back in R

```r
## NOTE: steps before the the conversion from ArchRProject to Signac SeuratObject.

#BiocManager::install("EnsDb.Hsapiens.v75")
#library(EnsDb.Hsapiens.v75)
#annotations <- getAnnotation(seqStyle = 'UCSC', refversion = 'hg19', reference = EnsDb.Hsapiens.v75)
#pm <- getPeakMatrix(ArchRProject= proj)

# Conversion function
seurat_atac <- ArchR2Signac(
  ArchRProject = proj,
  # samples = samples, # Provide a list of unique sample
  fragments_dir = fragments_dir, # the folder that contains all fragments samples in '.fragments.tsv.gz' or '.tsv.gz'
  pm = pm, # geting peak martix
  fragments_fromcellranger = "NO",
  fragments_file_extension = '.fragments.tsv.gz',
  refversion = 'hg19', # write the EnsDb version
  annotation = annotations
)

```

* **STEP 4 - Transfer ArchRProject gene score matrix to Signac SeuratObject.**

```r
gsm <- getGeneScoreMatrix(ArchRProject = proj, SeuratObject = seurat_atac)

seurat_atac[['RNA']] <- CreateAssayObject(counts = gsm)

```

* **STEP 5 - Transfer ArchRProject dimension reduction ("IterativeLSI", "IterativeLSI2" or "Harmony") and UMAP to Signac SeuratObject.**

```r
seurat_atac <- addDimRed(
  ArchRProject = proj,
  SeuratObject = seurat_atac,
  addUMAPs = "UMAP",
  reducedDims = "IterativeLSI"
) # default is "IterativeLSI"

#add both 'Harmony' and ‘IterativeLSI’:
seurat_atac <- addTwoDimRed(
  ArchRProject = proj,
  SeuratObject = seurat_atac,
  addUMAPs = "UMAP",
  reducedDims1 = "IterativeLSI",
       # Please limit your reducedDims to one of the following: IterativeLSI, IterativeLSI2 or Harmony
  reducedDims2 = "Harmony" # IterativeLSI2 or Harmony
)

```
