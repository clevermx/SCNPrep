#' convert old dataset to newer version
#'
#' @param path
#'
#' @return
#'
#' @import jsonlite
#' @import Matrix
#' @import rhdf5
#'
#' @examples
convertOldToNew <- function(path) {

  dataFileName <- "data.h5"
  expDataFileName <- "exp_data.json"

  filePath <- file.path(path, dataFileName)
  x <- as.numeric(h5read(filePath, "X/data"))
  p <- as.numeric(h5read(filePath, "X/indptr"))
  j <- as.numeric(h5read(filePath, "X/indices"))

  expDataJson <- fromJSON(file.path(path, expDataFileName))
  genes <- expDataJson$features
  barcodes <- expDataJson$barcodes

  counts <- sparseMatrix(
    j = j, p = p, x = x, repr="R", dims=c(length(genes), length(barcodes)), index1=F
  )

  rownames(counts) <- genes
  colnames(counts) <- barcodes

  countsCopy <- counts
  countsCopy@x <- rep(1, length(counts@x))
  featureCounts <- as.list(rowSums(countsCopy))
  featureCounts <- lapply(featureCounts, unbox)
  names(featureCounts) <- rownames(counts)

  jsonDataToOverwrite <- list(
    "expType"=unbox(expDataJson$expType),
    "features"=rownames(counts),
    "featureCounts"=featureCounts,
    "barcodes"=colnames(counts),
    "totalCounts"=colSums(counts)
  )

  write(toJSON(jsonDataToOverwrite), file.path(path, expDataFileName))
  message(sprintf("%s was overwriten", file.path(path, expDataFileName)))
}
