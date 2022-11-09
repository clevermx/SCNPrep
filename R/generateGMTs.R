#' Generate GMT files
#'
#' @param token token string used as prefix in all gmt titles
#' @param universe list of all expressed genes
#' @param species species
#' @param dflist list of data.frames - results of FindAllMarkers
#' @param pval pval threshold
#'
#' @return
#' @importFrom dplyr filter
#' @importFrom dplyr pull
#'
#' @examples
generateGMTs <- function(token, universe, species,
                         dflist, pval=0.05) {
  ## universe is all expressed genes in the dataset
  idTypeFrom <- guessIdType(universe)
  universeEntrez <- convertGeneIdsToEntrez(universe, idTypeFrom, species)


  universeName <- sprintf("%s#UNIVERSE", token)
  universeGMT <- sprintf("%s\t%s\t%s", universeName, universeName,
                         paste0(universeEntrez, collapse=","))
  gmts <- c(universeGMT)
  names <- c(universeName)

  for (i in 1:length(dflist)) {
    key <- names(dflist)[i]
    markers <- dflist[[i]]
    clusters <- levels(markers$cluster)
    for (clusterId in clusters) {

      goodMarkers <- markers %>%
        dplyr::filter(p_val_adj < pval & cluster==clusterId) %>%
        dplyr::pull(gene)
      goodMarkersEntrez <- convertGeneIdsToEntrez(goodMarkers, idTypeFrom, species)

      if (length(goodMarkersEntrez) > 0) {
        moduleName <- sprintf("%s#%s#%s", token, key, clusterId)
        gmtLine <- sprintf("%s\t%s\t%s", moduleName, universeName,
                           paste0(goodMarkersEntrez, collapse=","))
        gmts <- c(gmts, gmtLine)
        names <- c(names, moduleName)
      }
    }
  }

  names(gmts) <- names
  return(gmts)
}
