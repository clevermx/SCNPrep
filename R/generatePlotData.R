ALL_REDUCTIONS <- c("FItSNE", "tsne", "umap", "pca")

#' Generate data.frame with plot data
#'
#' @param object SeuratObject
#' @param userAnnotations list of data.frames with user-defined annotation for cells
#' @param maxReductionDims maximum number of dimensions for each dim-reduction technique to take.
#' Mostly used for reducing number of PCs reported in the object.
#'
#' @return
#'
#' @importFrom dplyr group_by
#' @importFrom dplyr summarize_at
#' @importFrom dplyr sym
#' @importFrom Seurat Idents
#' @import jsonlite
#'
#' @examples
generatePlotData <- function(object,
                             userAnnotations,
                             maxReductionDims,
                             generateCenters=T,
                             generateMasks=T) {
  presentAssays <- Assays(object)
  reductions <- ALL_REDUCTIONS
  reductions <- reductions[reductions %in% names(object@reductions)]

  if (length(reductions) == 0) {
    stop("No TSNE, UMAP or PCA found, please make sure to run some dimensionality reduction first")
  }

  embeddings <- lapply(reductions, function(red) {
    emb <- object@reductions[[red]]@cell.embeddings
    dimMax <- min(ncol(emb), maxReductionDims)
    reduced <- emb[, 1:dimMax]
    reduced
  })

  dataForPlot <- object@meta.data
  dataForPlot <- cbind(dataForPlot, embeddings)

  if ("orig.ident" %in% names(object@meta.data)) {
    if (length(levels(object$orig.ident)) > 1 || length(unique(object$orig.ident)) > 1) {
      dataForPlot$Sample <- object$orig.ident
      dataForPlot$Sample <- as.factor(dataForPlot$Sample)
    } else {
      samples <- gsub("(\\w+)(_|:)[NATGCx]+", "\\1", colnames(object))
      if (length(unique(samples)) > 1 && length(unique(samples)) < 100) {
        dataForPlot$Sample <- as.factor(samples)
      }
    }
  }

  if (!("Cluster" %in% colnames(dataForPlot))) {
    dataForPlot$Cluster <-  Seurat::Idents(object = object)
  }

  for (userAnnotation in userAnnotations) {
    dataForPlot <- cbind(dataForPlot, userAnnotation[rownames(dataForPlot), ])
  }


  colnames(dataForPlot) <- sapply(colnames(dataForPlot), function(colName){
    gsub("\\.", "_", colName)
  })

  clusterColnames <- grep("^(C|c)luster|_res", colnames(dataForPlot), value = T)
  fields <- list()

  for (column in colnames(dataForPlot)) {
    thisColumn <- dataForPlot[, column]
    if (is.numeric(thisColumn)) {
      fields[[column]] = list(
        type=unbox("numeric"),
        range=c(min(thisColumn), max(thisColumn))
      )
    }
    if (is.factor(thisColumn)) {
      fields[[column]] = list(
        type=unbox("factor"),
        levels=levels(thisColumn)
      )
    }
  }


  annotations <- list()

  for (reduction in reductions) {
    for (clusterCol in clusterColnames) {
      dimColnames <- colnames(object@reductions[[reduction]]@cell.embeddings)[1:2]


      if (generateCenters) {
        centers <- dataForPlot %>%
          dplyr::group_by(!!dplyr::sym(clusterCol)) %>%
          dplyr::summarize_at(dimColnames, median)
        centers <- as.data.frame(centers)
        centers$Text <- centers[, clusterCol]

        centersName = sprintf("%s_%s_centers", reduction, clusterCol)

        annotations[[centersName]] <- list(
          type=unbox("text"),
          value=unbox(clusterCol),
          coords=dimColnames,
          data=centers
        )
      }

      if (generateMasks) {
        masks <- maskEstimator(dataForPlot, clusterCol, dimColumns = dimColnames)
        bordersName = sprintf("%s_%s_borders", reduction, clusterCol)

        annotations[[bordersName]] <- list(
          type=unbox("polygon"),
          value=unbox("group"),
          coords=dimColnames,
          data=masks
        )
      }

    }
  }

  plotDataForJson = list(
    fields = fields,
    data = dataForPlot,
    annotations = annotations
  )
  return(plotDataForJson)
}
