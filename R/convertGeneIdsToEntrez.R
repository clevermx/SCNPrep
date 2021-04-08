#' Convert Gene Ids to Entrez
#'
#' @param ids ids to convert
#' @param from from which type to convert (symbol, ensembl)
#' @param species species
#'
#' @return
#' @export
#'
#' @examples
convertGeneIdsToEntrez <- function(ids, from, species) {
  converter <- converters[[from]][[species]]
  entrezIds <- converter[ids, "V2", drop=T]
  entrezIds <- na.omit(entrezIds)
  return(entrezIds)
}
