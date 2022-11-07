
ALLOWED_SPECIES <- c("mm", "hs", "rn")

H5_DATASET_FILE_NAME <- "data.h5"
DATASET_FILE_NAME <- "dataset.json"
EXP_DATA_FILE_NAME <- "exp_data.json"
MARKERS_FILE_NAME <- "markers.json"
PLOT_DATA_FILE_NAME <- "plot_data.json"
FILES_FOLDER_NAME <- "files"
GMT_PREFIX <- "%s.modules.gmt"
GMT_ANNOTATION_FILE <- "modules.annotation.json"

#' Migrate Seurat object to JSON files of scNavigator
#'
#' @param object Seurat object to convert
#' @param species Species of the dataset. Supported values are: "mm", "hs", and "rn"
#' @param assay Which assay to use to generate expression table and/or markers
#' @param slot Which slot to use from the given assay to generate expression table and/or markers
#' @param outdir Output folder to put the results. Will be created if doesn't exist
#' @param userAnnotation List of data.frame objects. Annotation to be added to the cells. Rownames must be
#' the same as the barcodes
#' @param token Unique token for scNavigator - if not given will be randomly generated.
#' @param name Name of the dataset
#' @param description Desription of the dataset
#' @param link External link for the dataset
#' @param public Should be dataset list on the main page in public. T/F
#' @param curated Should be dataset list on the main page in curated T/F
#' @param debug Enable test features for this dataset
#' @param markers markers - data.frame or list of data.frames with generated markers
#' @param generateGMTS if markers are not given, should we calculate DE and generate GMTS? T/F
#'
#' @return
#'
#' @import Matrix
#' @import stringi
#' @import Seurat
#' @import jsonlite
#'
#' @examples
#' @export
migrateSeuratObject <- function(
  object,
  species,
  assay="RNA",
  slot="counts",
  outdir=".",
  userAnnotations=list(),

  token=NULL,
  name="",
  description="",
  link="",
  public=F,
  curated=F,
  debug=F,


  generateCenters=T,
  generateMasks=T,

  markers=NULL,
  generateMarkers=F,
  generateGMTS=public && (!is.null(markers) || generateMarkers),
  ## TODO: implement this maxReductionDims, PLEASE
  maxReductionDims=5,
  maxCellsPerIdent=200,
  compressionLevel=9
) {

  ## Argument checkers

  if (typeof(object) != "S4" | class(object) != "Seurat") {
    stop("Object doesn't seem to be Seurat object")
  }
  if (!(species %in% ALLOWED_SPECIES)) {
    stop(sprintf("Unsupported species - %s. \n Supported species are: %s",
             species, paste0(ALLOWED_SPECIES, collapse = ", ")))
  }

  presentAssays <- Assays(object)
  if (!(assay %in% presentAssays)) {
    stop(sprintf("Assay %s is not present in Seurat object"), assay)
  }

  dir.create(outdir, recursive = T)

  plotDataForJson <- generatePlotData(object,
                                      userAnnotations,
                                      maxReductionDims=maxReductionDims,
                                      generateCenters=generateCenters,
                                      generateMasks=generateMasks)
  write(toJSON(plotDataForJson), file.path(outdir, PLOT_DATA_FILE_NAME))
  message(sprintf("%s - generated", file.path(outdir, PLOT_DATA_FILE_NAME)))

  counts <- GetAssayData(object, slot=slot, assay=assay)
  expDataForJson <- writeH5ExpressionData(counts, file.path(outdir, H5_DATASET_FILE_NAME), compressionLevel=compressionLevel)
  write(toJSON(expDataForJson), file.path(outdir, EXP_DATA_FILE_NAME))
  message(sprintf("%s - generated", file.path(outdir, EXP_DATA_FILE_NAME)))

  if (is.null(token)) {
    token <- stri_rand_strings(1, 20)
  }

  datasetDescrptor <- list(
    "token"=unbox(token),
    "name"=unbox(name),
    "description"=unbox(description),
    "link"=unbox(link),
    "species"=unbox(species),
    "cells"=unbox(ncol(object)),
    "public"=unbox(public),
    "curated"=unbox(curated),
    "debug"=unbox(debug)
  )
  write(toJSON(datasetDescrptor, pretty=T),
        file.path(outdir, DATASET_FILE_NAME))
  message(sprintf("%s - generated", file.path(outdir, DATASET_FILE_NAME)))


  mmarkers <- NULL
  if (!is.null(markers)) {
    if (class(markers) == "list") {
      keys <- names(markers)
      mmarkers <- markers
      totalMarkers <- generateMarkers(mmarkers, keys)
      write(toJSON(totalMarkers), file.path(outdir, MARKERS_FILE_NAME))
      message(sprintf("%s - generated", file.path(outdir, MARKERS_FILE_NAME)))

    } else if (class(markers) == "data.frame") {
      mmarkers <- list("markers"=markers)
      totalMarkers <- generateMarkers(mmarkers, "markers")
      write(toJSON(totalMarkers), file.path(outdir, MARKERS_FILE_NAME))
      message(sprintf("%s - generated", file.path(outdir, MARKERS_FILE_NAME)))
    } else {
      stop("Provided markers are of the wrong class")
    }
  } else {
    if (generateMarkers) {
      markers <- FindAllMarkers(
        object, assay=assay, slot=slot, only.pos=T, max.cells.per.ident=maxCellsPerIdent
      )
      mmarkers <- list("markers"=markers)
      totalMarkers <- generateMarkers(mmarkers, "markers")
      write(toJSON(totalMarkers), file.path(outdir, MARKERS_FILE_NAME))
      message(sprintf("%s - generated", file.path(outdir, MARKERS_FILE_NAME)))
    }
  }

  if (generateGMTS) {
    if (is.null(mmarkers)) {
      stop("GenerateGMTs is true, but markers were never provided")
    }
    outfile <- file.path(outdir, sprintf(GMT_PREFIX, species))
    gmts <- generateGMTs(token, rownames(counts), species, totalMarkers)
    writeLines(gmts, outfile)
    message(sprintf("%s - generated", outfile))

    totalJson <- list()
    for (gmt in names(gmts)) {
      totalJson[[gmt]] <- list(
        "token"=unbox(token),
        "name"=unbox(name),
        "title"=unbox(description),
        "species"=unbox(species),
        "link"=unbox(link)
      )
    }
    annotationsFile <- file.path(outdir, GMT_ANNOTATION_FILE)
    write(toJSON(totalJson, pretty=T), annotationsFile)
    message(sprintf("%s - generated", annotationsFile))
  }

  filesFolder <- file.path(outdir, FILES_FOLDER_NAME)
  dir.create(filesFolder, showWarnings = F, recursive = T)
  message(sprintf("Folder %s - generated", filesFolder))

}
