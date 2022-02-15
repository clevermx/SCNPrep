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

  ## rewrite markers if needed

  markersFile <- file.path(path, MARKERS_FILE_NAME)
  if (file.exists(markersFile)) {
    markers <- fromJSON(markersFile)
    newMarkers <- list()

    for (i in 1:length(markers)) {
      table <- markers[[i]]
      table$cluster <- as.factor(table$cluster)

      if ("avg_log2FC" %in% names(table)) {
        j <- which(names(table) == "avg_log2FC")
        names(table)[j] <- "avg_logFC"
      }

      newMarkers[[i]] <- table
    }
    names(newMarkers) <- names(markers)

    write(toJSON(newMarkers), markersFile)
    message(sprintf("%s was overwriten", markersFile))
  }


  write(toJSON(jsonDataToOverwrite), file.path(path, expDataFileName))
  message(sprintf("%s was overwriten", file.path(path, expDataFileName)))
}
