#' write H5 expression data to a file
#'
#' @param counts counts table (will be converted to RsparseMatrix)
#' @param newH5File file name for new H5 file
#'
#' @return
#'
#' @import rhdf5
#' @import Matrix
#' @import jsonlite
#'
#' @examples
writeH5ExpressionData <- function(counts, newH5File) {
  asIs <- F
  if (sum(dim(counts)) == 0) {
    stop("Provided slot in provided assay does not exists, cannot generate data.h5")
  }

  if ("dgCMatrix" %in% class(counts)) {
    # converting to row-based sparse matrix
    counts <- as(counts, "RsparseMatrix")
  } else if ("matrix" %in% class(counts)) {
    ## converting to sparse anyway
    counts <- as(counts, "RsparseMatrix")
    asIs <- T
  }

  h5createFile(newH5File)
  h5createGroup(newH5File, "X")

  genes <- nrow(counts)
  barcodes <- ncol(counts)
  nonZero <- length(counts@x)

  h5createDataset(newH5File, "X/indptr", c(genes + 1),
                  storage.mode = "integer", level=9)
  h5write(counts@p, newH5File, "X/indptr")

  h5createDataset(newH5File, "X/indices", c(nonZero),
                  storage.mode = "integer", level=9)
  h5write(counts@j, newH5File, "X/indices")

  h5createDataset(newH5File, "X/data", c(nonZero),
                  storage.mode = "double", level=9)
  h5write(as.double(counts@x), newH5File, "X/data")
  h5closeAll()

  id <- H5Fopen(newH5File)
  group <- H5Gopen(id, "X")
  h5writeAttribute(c(barcodes, genes), group, "shape")
  h5closeAll()

  message(sprintf("%s - generated", newH5File))

  return(list(
    "expType"=unbox(asIs),
    "features"=rownames(counts),
    "barcodes"=colnames(counts),
    "totalCounts"=colSums(counts)
  ))
}
