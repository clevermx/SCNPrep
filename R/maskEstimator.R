#' Estimating masks for plots
#'
#' For a grid of dimColumns estimated borders for points grouped by maskField.
#'
#' @return
#' @import data.table
#' @importFrom MASS kde2d
#' @import igraph
#' @import alphahull
#' @importFrom reshape2 melt
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

  classes <- get(maskField, data)
  occurrences <- table(classes)
  goodClasses <- names(occurrences[occurrences > 2])


  z <- MASS::kde2d(data[, dimColumns[1]], data[, dimColumns[2]], n=gridSize,
                   lims = c(range(data[, dimColumns[1]]) * 1.05,
                            range(data[, dimColumns[2]]) * 1.05))
  gridX <- z$x
  gridY <- z$y
  fullGrid <- expand.grid(gridX, gridY)

  # lets estimate two-dimensional kernel density for each cluster
  allDensities <- lapply(sort(goodClasses), function(i) {
    subset <- data %>% filter(.data[[maskField]] == i)
    z <- MASS::kde2d(subset[, dimColumns[1]], subset[, dimColumns[2]], n=gridSize, lims = c(range(data[, dimColumns[1]]) * 1.05, range(data[, dimColumns[2]]) * 1.05))
    as.numeric(z$z)
  })
  allDensities <- do.call(cbind, allDensities)
  colnames(allDensities) <- sort(goodClasses)

  # for each point in the grid find the cluster with higher probability
  # higher than the threshold
  # cluster -1 is the outside
  allDensitiesMax <- apply(allDensities, 1, function(x) {
    tt <- which.max(x)
    ifelse(x[tt] > threshold, tt, -1)
  })

  ## object for ggplot plotting
  allDensitiesT <- cbind(fullGrid, allDensitiesMax)
  colnames(allDensitiesT) <- c(dimColumns, maskField)
  allDensitiesT <- data.table(allDensitiesT)

  # getting grid graph
  gridM <- matrix(get(maskField, allDensitiesT), ncol=gridSize)
  graph_t <- as.data.table(reshape2::melt(gridM))
  setnames(graph_t, c("x", "y", "id"))
  graph_t <- graph_t[id != -1]
  graph_t[, name := paste0(x, "_", y)]
  graph_t <- graph_t[, list(name, id, x, y)]
  node_colors <- graph_t[, setNames(id, name)]

  edges_t <- rbind(graph_t[, list(from=name, to=paste0(x+1, "_", y))],
                   graph_t[, list(from=name, to=paste0(x-1, "_", y))],
                   graph_t[, list(from=name, to=paste0(x, "_", y+1))],
                   graph_t[, list(from=name, to=paste0(x, "_", y-1))]
  )
  edges_t <- edges_t[to %in% graph_t$name]

  edges_t <- edges_t[node_colors[from] == node_colors[to]]

  g <- graph_from_data_frame(edges_t, directed=FALSE, vertices = graph_t)

  gcomp <- components(g)
  gcomp$csize

  shapes <- rbindlist(lapply(seq_len(gcomp$no), function(i) {
    if (gcomp$csize[i] < 20) {
      return(NULL)
    }
    piece <- V(g)[gcomp$membership == i]
    zz <- graph_t[name %in% piece$name]
    zzz <- data.frame(x=gridX[zz$x]+rnorm(nrow(zz), sd = 0.001),
                      y=gridY[zz$y]+rnorm(nrow(zz), sd = 0.001))
    shape <- ashape(zzz, alpha=2)
    res <- data.table(shape$edges)[, list(x1, y1, x2, y2)]
    ggt <- data.table(from=paste0("v", shape$edges[, "ind1"]),
                      to=paste0("v", shape$edges[, "ind2"]))

    gg <- graph_from_data_frame(shape$edges[, c("ind1", "ind2")],directed=F)
    components(gg)$csize
    V(gg)$membership <- components(gg)$membership
    new_order <- igraph::dfs(gg, 1)$order
    new_order1 <- do.call("c", lapply(split(new_order, new_order$membership), function(x) c(x, x[1])))
    new_order1i <- as.integer(new_order1$name)
    new_order1m <- new_order1$membership

    shape_x <- shape$x[new_order1i, 1:2]
    colnames(shape_x) <- dimColumns

    res <- data.table(shape_x)
    res[, c(maskField) := list(sort(goodClasses)[piece$id[1]])]
    res[, group := paste0("gr", i, "_", new_order1m)]
    res[]
  }))

  return(shapes)
}
