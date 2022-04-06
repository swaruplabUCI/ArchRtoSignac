#' loadinglibrary
#'
#' This function a quick way to check if required packages for ArchRtoSignac package are installed and load those available packages automatically.
#'
#' @param packages_x Libraries required for the use of ArchRtoSignac
#' @export
#' @examples
#' loadinglibrary(dependencies)
loadinglibrary <- function(
  packages_x = NULL
) {
  for (i in 1:length(packages_x)) {
    loading_package <- packages_x[i]
    if (!require(loading_package, character.only = TRUE)) {
      print(paste("Package", loading_package, 'not found. Please Installing Package!'))
    }
    else {
      print(paste0("Loading Package: ", loading_package))
      suppressPackageStartupMessages(library(loading_package, character.only = TRUE))
      print(paste0("Package: ", loading_package, " -- Loaded Successfully"))
    }
  }
}

#' getPeakMatrix
#'
#' This function gets fixed-width peak matrix from ArchR project and change the row names of peak matrix to their matched chromosome range
#'
#' @param ArchRProject A ArchRProject
#' @export
#' @examples
#' pm <- getPeakMatrix(proj)
getPeakMatrix <- function(
  ArchRProject = NULL
){
  print("In Progress:")
  print("get Matrix From ArchRProject")
  peak_matrix <- ArchR::getMatrixFromProject(ArchRProject, useMatrix='PeakMatrix')
  pm <- assays(peak_matrix)$PeakMatrix # peaks sparse martix
  rownames(pm) <- paste0(as.character(seqnames(ArchRProject@peakSet)), '-', as.character(start(ArchRProject@peakSet)), '-', as.character(end(ArchRProject@peakSet)))

  print("Return peak matrix")
  pm
}

#' getAnnotation
#'
#' This function gets the gene annotation from Ensembl Database in GRanges Object for the Seurat object, which includes information related to genomic locations and their associated annotations
#' Then it changes the annotation to the UCSC style and hg38, since Cellranger ATAC count input FastQ files were mapped to GRCh38.
#'
#' @param reference An Ensembl genome reference used for function GetGRangesFromEnsDb to extract gene annotations from EnsDb, for example: EnsDb.Hsapiens.v86
#' @param seqStyle A default sequence style changes the annotation extracted from EnsDb to ‘UCSC’ since Signac maps to hg38
#' @param refversion The assembly release and versions of UCSC genome reference
#' @export
#' @examples
#' annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, seqStyle = 'UCSC', refversion = 'hg38')

getAnnotation<- function(
  reference = NULL, # GetGRangesFromEnsDb requires an EnsDb, for example: EnsDb.Hsapiens.v86
  seqStyle = 'UCSC', # change to UCSC style
  refversion = NULL # write the EnsDb version for example: 'hg38'
){
  print("In Progress:")
  print("Extract genomic ranges from EnsDb object and prepare annotation")
  #annotations <- Signac::GetGRangesFromEnsDb(ensdb = reference)
  annotations <- suppressWarnings(Signac::GetGRangesFromEnsDb(ensdb = reference))
  # Warnings:
  # 'In .Seqinfo.mergexy(x, y) :
  #     The 2 combined objects have no sequence levels in common. (Use
  #     suppressWarnings() to suppress this warning.)')
  # # https://github.com/timoast/signac/issues/429

  ## Change the seqlevels to UCSC style.
  seqlevelsStyle(annotations) <- seqStyle
  genome(annotations) <- refversion

  print("Return Annotation")
  annotations
}

#' ArchR2Signac
#'
#' This function converts the ArchRProject to SeuratObject by creating a list of seurat objects for each sample with their corresponding peak matrix and then merge objects from each sample in the list created using Swarup lab GitHub wrapper function
#'
#' @param ArchRProject A ArchRProject
#' @param reference A reference used for function GetGRangesFromEnsDb to extract gene annotations from EnsDb
#' @param seqStyle A default seq style changes the annotation extracted to UCSC style since the ArchRProject was mapped to hg38
#' @param refversion The assembly release and versions of UCSC genome reference
#' @param samples A unique sample List for all processed samples from ArchRProject
#' @param fragments_dir A PATH to the cellranger output, the folder that contains all samples folders, not the one with '/outs/fragments.tsv.gz'.
#' @param pm A peak matrix. getMatrixFromProject function extracts peak matrix from ArchRProject and it needs to be reformat into Signac style
#' @param output_dir '/outs/' the PATH before 'fragments.tsv.gz'
#' @param annotation annotation got from getAnnotation() for SeuratObject
#' @export
#' @examples
#' seurat_atac <- ArchR2Signac(ArchRProject = proj1, samples = samples, fragments_dir = fragments_dir, pm = pm, output_dir = '/outs/', seqStyle = 'UCSC', refversion = 'hg38', reference = EnsDb.Hsapiens.v86, annotation = annotations)

ArchR2Signac <- function(
  ArchRProject = NULL,
  samples = "UniqueSampleList", # Provide a list of unique sample
  fragments_dir = "PATHtoCellRangerOutput", # directory of the cellranger output, the folder that contains all samples
  pm = pm, # geting peak martix
  output_dir = '/outs/',
  seqStyle = 'UCSC',
  refversion = 'hg38', # write the EnsDb version
  reference = EnsDb.Hsapiens.v86, # choose the EnsDb as the reference
  annotation = annotations # annotation from getAnnotation()
 ){
   print("In Progress:")
   print("Prepare Seurat list for each sample")
   seurat_list <- lapply(samples, function(cur_sample){
     print(cur_sample)
     #print out the sample name in progress
     cur_fragments <- paste0(fragments_dir, cur_sample, output_dir, 'fragments.tsv.gz')
     # seeking the pattern matched in colnames(pm); metadata of the corresponding sample
     cur_pm <- pm[,grepl(paste0(cur_sample, '#'), colnames(pm))]
     cur_meta <- ArchRProject@cellColData %>% as.data.frame %>% subset(Sample == cur_sample)

     # change colnames and rowname format:
     colnames(cur_pm) <- do.call(rbind, str_split(colnames(cur_pm), '#'))[,2]
     rownames(cur_meta) <- do.call(rbind, str_split(rownames(cur_meta), '#'))[,2]
     print(dim(cur_pm))
     # create chromatin assay
     cur_chromatin <- Signac::CreateChromatinAssay(
       counts=cur_pm, # should we add data instead counts
       sep = c('-', '-'),
       fragments=cur_fragments, # do we need this?
       ranges=ArchRProject@peakSet,
       genome=refversion,
       annotation = annotation
     )

     # create a new Seurat obj with only the archR peaks:
     cur_atac <- Seurat::CreateSeuratObject(
       cur_chromatin,
       assay='peaks',
       meta.data = cur_meta,
     )

   })

   print("In Progress:")
   print("Merge Seurat list")
   # merge objects
   SeuratObject <- merge(
     x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)],
     add.cell.ids = samples
   )

   print("Return SeuratObject")
   SeuratObject

}

#' getGeneScoreMatrix
#'
#' This function gets gene score matrix from ArchR project and change the row names of gene score matrix to their matched gene features
#'
#' @param ArchRProject A ArchRProject
#' @param SeuratObject A Seurat object
#' @export
#' @examples
#' gsm <- getGeneScoreMatrix(ArchRProject = proj, SeuratObject = seurat_atac)
getGeneScoreMatrix <- function(
  ArchRProject = NULL,
  SeuratObject = NULL
){
  print("In Progress:")
  print("get Gene Score Matrix From ArchRProject")
  GeneScore_matrix <- ArchR::getMatrixFromProject(ArchRProject, useMatrix='GeneScoreMatrix')
  gsm <- assays(GeneScore_matrix)$GeneScoreMatrix # peaks sparse martix
  print("get Gene Features From ArchRProject")
  GeneFeatures <-getFeatures(
    ArchRProj = ArchRProject,
    useMatrix = "GeneScoreMatrix",
    select = NULL,
    ignoreCase = TRUE
  )
  # set the column names for gsm
  colnames(gsm) <- gsub("#", "_", colnames(gsm))
  ix <- match(colnames(SeuratObject), colnames(gsm))
  gsm <- gsm[,ix]

  print("Saving Gene Features From ArchRProject into Gene Score Matrix")
  rownames(gsm) <- GeneFeatures

  print("Return Gene Score Matrix")
  gsm
}

#' addDimRed
#'
#' This function adds dimension reduction, 'Harmony' or 'IterativeLSI' and UMAP to SeuratObject
#'
#' @param ArchRProject A ArchRProject
#' @param SeuratObject A Seurat object
#' @param reducedDims A reduction dimension can be transfered from from ArchRProject to Signac SeuratObject
#' @export
#' @examples
#' seurat_atac <- addDimRed(ArchRProject = proj, SeuratObject = seurat_atac,reducedDims = 'IterativeLSI')
addDimRed <- function(
  ArchRProject = NULL,
  SeuratObject = NULL,
  reducedDims = 'IterativeLSI'
){
  print("In Progress:")
  print("add UMAP From ArchRProject to SeuratObject")
  umap_df <- ArchRProject@embeddings$UMAP$df %>% as.matrix
  rownames(umap_df) <- colnames(SeuratObject) # make the rowname the same format as seurat
  colnames(umap_df) <- c('UMAP_1', 'UMAP_2')

  SeuratObject@reductions$umap <- Seurat::CreateDimReducObject(
      embeddings=umap_df,
      assay="peaks"
  )

  print("In Progress:")
  print("add reduction From ArchRProject to SeuratObject")
  if(reducedDims == 'Harmony'){
    harmony_matrix <- ArchRProject@reducedDims$Harmony$matDR
    rownames(harmony_matrix) <- colnames(SeuratObject)
    colnames(harmony_matrix) <- paste0('LSI_', 1:ncol(harmony_matrix))

    SeuratObject@reductions$harmony <- Seurat::CreateDimReducObject(
        embeddings=harmony_matrix,
        assay="peaks"
    )
  } else if(reducedDims == 'IterativeLSI'){
    LSI_matrix <- ArchRProject@reducedDims$IterativeLSI$matSVD
    rownames(LSI_matrix) <- colnames(SeuratObject)
    colnames(LSI_matrix) <- paste0('LSI_', 1:ncol(LSI_matrix))

    SeuratObject@reductions$IterativeLSI <- Seurat::CreateDimReducObject(
        embeddings=LSI_matrix,
        assay="peaks"
    )
  }

  print("Return SeuratObject")
  SeuratObject

}
