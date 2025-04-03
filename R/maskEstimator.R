#' Estimating masks for plots
#'
#' For a grid of dimColumns estimated borders for points grouped by maskField.
#'
#' @return
#' @import data.table
#' @importFrom mascarade generateMask
#' @export
#'
#' @examples
#' \dontrun{
#' clusterMasking <- maskEstimator(dataForPlot, "Cluster", dimColumns=c("UMAP_1", "UMAP_2"))
#' }
maskEstimator <- function(data, maskField,
                          gridSize=200,
                          threshold=0.0001,
                          dimColumns=c("tSNE_1", "tSNE_2")) {

  dims <- cbind(data[, dimColumns[1]], data[, dimColumns[2]])
  colnames(dims) <- dimColumns
  rownames(dims) <- rownames(data)
  
  clusters <- get(maskField, data)
  
  maskTable <- mascarade::generateMask(dims=dims, clusters=clusters)
  maskTable <- maskTable[, c(dimColumns, "cluster", "group"), with=FALSE]
  data.table::setnames(maskTable, "cluster", maskField)
  return(maskTable)
}
