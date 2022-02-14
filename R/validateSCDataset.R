#' validateSCDataset
#'
#' @param folder
#'
#' @return
#' @export
#'
#' @examples
validateSCDataset <- function(folder) {

  descriptorFile <- file.path(folder, DATASET_FILE_NAME)
  plotDataFile <- file.path(folder, PLOT_DATA_FILE_NAME)

  if (!file.exists(descriptorFile)) {
    stop(sprintf("Error: descriptor file `%s` does not exist", DATASET_FILE_NAME))
  }

  if (!file.exists(plotDataFile)) {
    stop(sprintf("Error: plot data file `%s` does not exist", PLOT_DATA_FILE_NAME))
  }

  descriptor <- fromJSON(descriptorFile)
  plotData <- fromJSON(plotDataFile, simplifyVector = F)

  validateDatasetDescriptor(descriptor)
  validatePlotData(plotData)


  h5path <- file.path(folder, H5_DATASET_FILE_NAME)
  expDataFile <- file.path(folder, EXP_DATA_FILE_NAME)

  if (!file.exists(h5path) && file.exists(expDataFile)) {
    stop(sprintf("Error: files `%s` and `%s` should be both present. Now `%s` is missing",
                 h5path, expDataFile, h5path))
  }

  if (file.exists(h5path) && !file.exists(expDataFile)) {
    stop(sprintf("Error: files `%s` and `%s` should be both present. Now `%s` is missing",
                 h5path, expDataFile, expDataFile))
  }

  if (file.exists(h5path) && file.exists(expDataFile)) {
    expDataJson <- fromJSON(expDataFile)

    if (descriptor$cells == length(plotData$data)) {
      message("Descriptor `cells` field is consistent with plot data")
    } else {
      stop(sprintf("Error: descriptor `cells` field is inconsistent with plot data: %d and %d respectively",
                   descriptor$cells, length(plotData$data)))
    }

    validateExpressionData(expDataJson, h5path)

    if (descriptor$cells == length(expDataJson$barcodes)) {
      message("Descriptor `cells` field is consistent with expression data")
    } else {
      stop(sprintf("Error: descriptor `cells` field is inconsistent with expression data: %d and %d respectively",
                   descriptor$cells, length(expDataJson$barcodes)))
    }
  }

  if (!file.exists(h5path) && !file.exists(expDataFile)) {
    message("Expression data is not present")
  }

  markersFile <- file.path(folder, MARKERS_FILE_NAME)

  if (file.exists(markersFile)) {
    markers <- fromJSON(file.path(folder, MARKERS_FILE_NAME), simplifyVector = F)
    validateMarkersJson(markers)
  } else {
    message("Markers data is not present")
  }


  gmtFile <- file.path(folder, sprintf(GMT_PREFIX, descriptor$species))
  gmtAnnotationFile <- file.path(folder, GMT_ANNOTATION_FILE)

  if (!file.exists(gmtFile) && file.exists(gmtAnnotationFile)) {
    stop(sprintf("Error: files `%s` and `%s` should be both present. Now `%s` is missing",
                 gmtFile, gmtAnnotationFile, gmtFile))
  }

  if (file.exists(gmtFile) && !file.exists(gmtAnnotationFile)) {
    stop(sprintf("Error: files `%s` and `%s` should be both present. Now `%s` is missing",
                 gmtFile, gmtAnnotationFile, gmtAnnotationFile))
  }

  if (file.exists(gmtFile) && file.exists(gmtAnnotationFile)) {
    gmtAnnotation <- fromJSON(gmtAnnotationFile, simplifyVector = F)
    gmtInfo <- parseGMTFile(gmtFile)
    validateGMTFiles(gmtInfo, gmtAnnotation)
  }

  if (!file.exists(gmtFile) && !file.exists(gmtAnnotationFile)) {
    message("Gene Signature data is not present")
  }


}
