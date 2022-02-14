#' convert old dataset to newer version
#'
#' @param Folder
#'
#' @return
#'
#' @import jsonlite
#' @import Matrix
#' @import rhdf5
#'
#' @examples
#' @export
folderIsConsistent <- function(folder) {

  descriptor <- fromJSON(file.path(folder, DATASET_FILE_NAME))
  plotData <- fromJSON(file.path(folder, PLOT_DATA_FILE_NAME), simplifyVector = F)

  h5dataset <- file.path(folder, "data.h5")
  x <- as.numeric(h5read(h5dataset, "X/data"))
  p <- as.numeric(h5read(h5dataset, "X/indptr"))
  j <- as.numeric(h5read(h5dataset, "X/indices"))

  id <- H5Fopen(h5dataset)
  group <- H5Gopen(id, "X")
  shapeId <- H5Aopen(group, "shape")
  shape <- H5Aread(shapeId)
  h5closeAll()

  counts <- sparseMatrix(
    j = j, p = p, x = x, repr="R", dims=c(shape[2], shape[1]), index1=F
  )

  expDataJson <- fromJSON(file.path(folder, "exp_data.json"))
  genes <- expDataJson$features
  barcodes <- expDataJson$barcodes
  shapeExp <- c(length(barcodes), length(genes))

  return(all(shape==shapeExp))
}
