#' loadinglibrary
#'
#' This function is a quick way to check if the required packages for ArchRtoSignac package are installed. If the packages are available, they are loaded automatically.
#'
#' @param packages Vector of libraries required for ArchRtoSignac (most dependencies should be installed when installing ArchRtoSignac)
#' @export
#' @examples
#' loadinglibrary(dependencies)
loadinglibrary <- function(
  packages
) {
  for (i in 1:length(packages)) {
    loading_package <- packages[i]
    if (!require(loading_package, character.only = TRUE)) {
      print(paste("Package", loading_package, 'not found. Please Install the Package!!'))
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
#' This function gets the fixed-width peak matrix from the ArchRProject and changes the row names of the peak matrix to their matched chromosome range.
#'
#' @param ArchRProject An ArchRProject
#' @export
#' @examples
#' pm <- getPeakMatrix(proj)
getPeakMatrix <- function(
  ArchRProject
){
  print("In Progress:")
  print("Get Matrix From ArchRProject")
  peak_matrix <- ArchR::getMatrixFromProject(ArchRProject, useMatrix='PeakMatrix')
  pm <- assays(peak_matrix)$PeakMatrix # peaks sparse martix
  rownames(pm) <- paste0(as.character(seqnames(ArchRProject@peakSet)), '-', as.character(start(ArchRProject@peakSet)), '-', as.character(end(ArchRProject@peakSet)))

  print("Return peak matrix")
  pm
}

#' getAnnotation
#'
#' This function gets the gene annotation, which includes information related to genomic locations and their associated annotations, from Ensembl Database in GRanges Object for the Seurat object.
#' Then it changes the annotation to the UCSC style (default).
#'
#' @param reference An Ensembl genome reference used for the Signac function GetGRangesFromEnsDb to extract gene annotations from EnsDb (for example: EnsDb.Hsapiens.v86)
#' @param seqStyle Sequence style to change the annotation extracted from EnsDb to (default is ‘UCSC’)
#' @param refversion The assembly release and versions of UCSC genome reference (for example: 'hg38')
#' @export
#' @examples
#' annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, seqStyle = 'UCSC', refversion = 'hg38')

getAnnotation<- function(
  reference, # GetGRangesFromEnsDb requires an EnsDb, for example: EnsDb.Hsapiens.v86
  seqStyle = 'UCSC', # change to UCSC style
  refversion # write the EnsDb version for example: 'hg38'
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
#' This function converts the ArchRProject to a SeuratObject by creating a list of Seurat objects for each sample with their corresponding peak matrix and then merging all of the objects in the list.
#'
#' @param ArchRProject An ArchRProject
#' @param refversion The assembly release and versions of UCSC genome reference
#' @param samples List of all the samples from the ArchRProject
#' @param fragments_dir PATH to the cellranger-atac output--the folder that contains all samples folders, not the one with '/outs/fragments.tsv.gz'.
#' @param pm Peak matrix (output) from the function getPeakMatrix
#' @param fragments_fromcellranger This is an Yes or No selection ("NO" | "N" | "No" or "YES" | "Y" | "Yes")
#' @param fragments_file_extension File_Extension for fragments files (typically they should be '.tsv.gz' or '.fragments.tsv.gz')
#' @param annotation annotation from the function getAnnotation()
#' @export
#' @examples
#' seurat_atac <- ArchR2Signac(ArchRProject = proj1, refversion = 'hg38', samples = samples, fragments_dir = fragments_dir, pm = pm, fragments_fromcellranger = "Yes", fragments_file_extension = NULL, annotation = annotations)

ArchR2Signac <- function(
  ArchRProject,
  refversion, # write the EnsDb version
  samples = NULL, # Provide a list of unique sample
  fragments_dir = NULL, # directory of the cellranger output, the folder that contains all samples
  pm, # geting peak martix
  fragments_fromcellranger = NULL, # "NO" | "N" | "No" or "YES" | "Y" | "Yes"
  fragments_file_extension = NULL, #  '.tsv.gz' or '.fragments.tsv.gz'
  # output_dir = '/outs/', # removal due to the input format for snapATAC, added when fragments_fromcellranger == "NO" | "N" | "No"
  annotation # annotation from getAnnotation()
 ){
   if (is.null(samples)){
     samples <- unique(ArchRProject@cellColData$Sample)
   }

   if(fragments_fromcellranger == "YES" | fragments_fromcellranger == "Y" | fragments_fromcellranger == "Yes") {
     print("In Progress:")
     print("Prepare Seurat list for each sample")

     output_dir = '/outs/'

     seurat_list <- lapply(samples, function(cur_sample){
       print(cur_sample)
       #print out the sample name in progress
       cur_fragments <- ifelse(is.list(fragments_dir),
                               paste0(fragments_dir[[which(samples == cur_sample)]], output_dir, 'fragments.tsv.gz'),
                               paste0(fragments_dir, cur_sample, output_dir, 'fragments.tsv.gz'))

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
   }

   if(fragments_fromcellranger == "NO" | fragments_fromcellranger == "N" | fragments_fromcellranger == "No") {

     print("IF selecting NO, please make sure to provide fragments_file_extension")
     print("In Progress:")
     print("Prepare Seurat list for each sample")

     seurat_list <- lapply(samples, function(cur_sample){
       print(cur_sample)
       #print out the sample name in progress
       cur_fragments <- ifelse(is.list(fragments_dir),
                        paste0(fragments_dir[[which(samples == cur_sample)]], fragments_file_extension),
                        paste0(fragments_dir, cur_sample, fragments_file_extension))
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

   }


   print("In Progress:")
   print("Merge Seurat list")
   # merge objects
   # SeuratObject <- merge(
   #   x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)],
   #   add.cell.ids = samples
   # )
   # merge objects
   SeuratObject <- if (length(seurat_list) > 1) {
     merge(
       x = seurat_list[[1]],
       y = seurat_list[2:length(seurat_list)],
       add.cell.ids = samples
     )
   } else {
  seurat_list[[1]]
  }

   print("Return SeuratObject")
   SeuratObject

}



#' getGeneScoreMatrix
#'
#' This function gets the gene score matrix from an ArchRProject and changes the row names of gene score matrix to their matched gene features
#'
#' @param ArchRProject An ArchRProject
#' @param SeuratObject A Seurat object
#' @export
#' @examples
#' gsm <- getGeneScoreMatrix(ArchRProject = proj, SeuratObject = seurat_atac)
getGeneScoreMatrix <- function(
  ArchRProject,
  SeuratObject
){
  print("In Progress:")
  print("Get Gene Score Matrix From ArchRProject")
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
#' This function adds dimension reduction ('Harmony' and/or 'IterativeLSI') and UMAP to a SeuratObject
#'
#' @param ArchRProject An ArchRProject
#' @param SeuratObject A Seurat object
#' @param addUMAPs  add UMAPs by default
#' @param reducedDims 'Harmony' or 'IterativeLSI' or 'IterativeLSI2', the dimension reduction to be transfered from ArchRProject to Signac SeuratObject (default is 'IterativeLSI')
#' @export
#' @examples
#' seurat_atac <- addDimRed(ArchRProject = proj, SeuratObject = seurat_atac, reducedDims = 'IterativeLSI') #reducedDims == c('IterativeLSI', 'Harmony')
addDimRed <- function(
  ArchRProject,
  SeuratObject,
  addUMAPs = "UMAP",
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
  } else if(reducedDims == 'IterativeLSI2'){
    LSI_matrix2 <- ArchRProject@reducedDims$IterativeLSI2$matSVD
    rownames(LSI_matrix2) <- colnames(SeuratObject)
    colnames(LSI_matrix2) <- paste0('LSI_', 1:ncol(LSI_matrix2))

    SeuratObject@reductions$IterativeLSI2 <- Seurat::CreateDimReducObject(
        embeddings=LSI_matrix2,
        assay="peaks"
    )
  }
  # else if(reducedDims == c('IterativeLSI', 'Harmony') | reducedDims == c('Harmony', 'IterativeLSI')){
  #
  #   print("In Progress: add IterativeLSI From ArchRProject to SeuratObject")
  #   LSI_matrix <- ArchRProject@reducedDims$IterativeLSI$matSVD
  #   rownames(LSI_matrix) <- colnames(SeuratObject)
  #   colnames(LSI_matrix) <- paste0('LSI_', 1:ncol(LSI_matrix))
  #
  #   SeuratObject@reductions$IterativeLSI <- Seurat::CreateDimReducObject(
  #       embeddings=LSI_matrix,
  #       assay="peaks"
  #   )
  #
  #   print("In Progress: add Harmony From ArchRProject to SeuratObject")
  #   harmony_matrix <- ArchRProject@reducedDims$Harmony$matDR
  #   rownames(harmony_matrix) <- colnames(SeuratObject)
  #   colnames(harmony_matrix) <- paste0('LSI_', 1:ncol(harmony_matrix))
  #
  #   SeuratObject@reductions$harmony <- Seurat::CreateDimReducObject(
  #       embeddings=harmony_matrix,
  #       assay="peaks"
  #   )
  #
  # }

  print("Return SeuratObject")
  SeuratObject

}

#' addCustomizeDimRed
#'
#' This function adds dimension reduction ('Harmony' and/or 'IterativeLSI') and UMAP to a SeuratObject with customized name of the dimension
#'
#' @param ArchRProject An ArchRProject
#' @param SeuratObject A Seurat object
#' @param addUMAPs  add UMAPs of the selected UMAPs
#' @param reducedDims the selected dimension reduction to be transfered from ArchRProject to Signac SeuratObject 
#' @param reducedDimsType 'Harmony' or 'IterativeLSI', this parameter tells the type of the reduced dimension (default is 'IterativeLSI')
#' @export
#' @examples
#' seurat_atac <- addCustomizeDimRed(ArchRProject = proj, SeuratObject = seurat_atac, addUMAPs = "UMAP", reducedDims = 'IterativeLSI', reducedDimsType = 'IterativeLSI') #reducedDims == c('IterativeLSI', 'Harmony')
addCustomizeDimRed <- function(
    ArchRProject,
    SeuratObject,
    addUMAPs = "UMAP",
    reducedDims = 'IterativeLSI',
    reducedDimsType = 'IterativeLSI'
){
  print("In Progress:")
  print("add UMAP From ArchRProject to SeuratObject")
  umap_df <- ArchRProject@embeddings[[addUMAPs]]$df %>% as.matrix
  rownames(umap_df) <- colnames(SeuratObject) # make the rowname the same format as seurat
  colnames(umap_df) <- c('UMAP_1', 'UMAP_2')
  
  SeuratObject@reductions$umap <- Seurat::CreateDimReducObject(
    embeddings=umap_df,
    assay="peaks"
  )
  
  print("In Progress:")
  print("add reduction From ArchRProject to SeuratObject")
  if(reducedDimsType == 'Harmony'){
    harmony_matrix <- ArchRProject@reducedDims[[reducedDims]][[1]]
    rownames(harmony_matrix) <- colnames(SeuratObject)
    colnames(harmony_matrix) <- paste0('LSI_', 1:ncol(harmony_matrix))
    
    SeuratObject@reductions$harmony <- Seurat::CreateDimReducObject(
      embeddings=harmony_matrix,
      assay="peaks"
    )
  } else if(reducedDimsType == 'IterativeLSI'){
    LSI_matrix <- ArchRProject@reducedDims[[reducedDims]][[1]]
    rownames(LSI_matrix) <- colnames(SeuratObject)
    colnames(LSI_matrix) <- paste0('LSI_', 1:ncol(LSI_matrix))
    
    SeuratObject@reductions$IterativeLSI <- Seurat::CreateDimReducObject(
      embeddings=LSI_matrix,
      assay="peaks"
    )
  } else {
    stop("Unsupported reduced dimension type, only 'Harmony' and 'IterativeLSI' are supported.")
  }
  
  print("Return SeuratObject")
  SeuratObject
  
}

#' addTwoDimRed
#'
#' This function adds dimension reduction ('Harmony' and/or 'IterativeLSI') and UMAP to a SeuratObject
#'
#' @param ArchRProject An ArchRProject
#' @param SeuratObject A Seurat object
#' @param addUMAPs  add UMAPs by default
#' @param reducedDims1 'Harmony'or 'IterativeLSI' or 'IterativeLSI2', the dimension reduction to be transfered from ArchRProject to Signac SeuratObject (default is 'IterativeLSI')
#' @param reducedDims2 'Harmony'or 'IterativeLSI' or 'IterativeLSI2', the dimension reduction to be transfered from ArchRProject to Signac SeuratObject (default is 'IterativeLSI')
#' @export
#' @examples
#' seurat_atac <- addTwoDimRed(ArchRProject = proj, SeuratObject = seurat_atac, reducedDims1 = 'IterativeLSI', reducedDims2= 'Harmony')
addTwoDimRed <- function(
  ArchRProject,
  SeuratObject,
  addUMAPs = "UMAP",
  reducedDims1 = NULL,
  reducedDims2 = NULL
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
  print("add reducedDims1 From ArchRProject to SeuratObject")
  if(reducedDims1 == 'Harmony'){
    harmony_matrix <- ArchRProject@reducedDims$Harmony$matDR
    rownames(harmony_matrix) <- colnames(SeuratObject)
    colnames(harmony_matrix) <- paste0('LSI_', 1:ncol(harmony_matrix))

    SeuratObject@reductions$harmony <- Seurat::CreateDimReducObject(
        embeddings=harmony_matrix,
        assay="peaks"
    )
  } else if(reducedDims1 == 'IterativeLSI'){
    LSI_matrix <- ArchRProject@reducedDims$IterativeLSI$matSVD
    rownames(LSI_matrix) <- colnames(SeuratObject)
    colnames(LSI_matrix) <- paste0('LSI_', 1:ncol(LSI_matrix))

    SeuratObject@reductions$IterativeLSI <- Seurat::CreateDimReducObject(
        embeddings=LSI_matrix,
        assay="peaks"
    )
  } else if(reducedDims1 == 'IterativeLSI2'){
    LSI_matrix2 <- ArchRProject@reducedDims$IterativeLSI2$matSVD
    rownames(LSI_matrix2) <- colnames(SeuratObject)
    colnames(LSI_matrix2) <- paste0('LSI_', 1:ncol(LSI_matrix2))

    SeuratObject@reductions$IterativeLSI2 <- Seurat::CreateDimReducObject(
        embeddings=LSI_matrix2,
        assay="peaks"
    )
  }

  print("In Progress:")
  print("add reducedDims2 From ArchRProject to SeuratObject")
  if(reducedDims2 == 'Harmony'){
    harmony_matrix <- ArchRProject@reducedDims$Harmony$matDR
    rownames(harmony_matrix) <- colnames(SeuratObject)
    colnames(harmony_matrix) <- paste0('LSI_', 1:ncol(harmony_matrix))

    SeuratObject@reductions$harmony <- Seurat::CreateDimReducObject(
        embeddings=harmony_matrix,
        assay="peaks"
    )
  } else if(reducedDims2 == 'IterativeLSI'){
    LSI_matrix <- ArchRProject@reducedDims$IterativeLSI$matSVD
    rownames(LSI_matrix) <- colnames(SeuratObject)
    colnames(LSI_matrix) <- paste0('LSI_', 1:ncol(LSI_matrix))

    SeuratObject@reductions$IterativeLSI <- Seurat::CreateDimReducObject(
        embeddings=LSI_matrix,
        assay="peaks"
    )
  }  else if(reducedDims2 == 'IterativeLSI2'){
    LSI_matrix2 <- ArchRProject@reducedDims$IterativeLSI2$matSVD
    rownames(LSI_matrix2) <- colnames(SeuratObject)
    colnames(LSI_matrix2) <- paste0('LSI_', 1:ncol(LSI_matrix2))

    SeuratObject@reductions$IterativeLSI2 <- Seurat::CreateDimReducObject(
        embeddings=LSI_matrix2,
        assay="peaks"
    )
  }

  print("Return SeuratObject")
  SeuratObject

}
