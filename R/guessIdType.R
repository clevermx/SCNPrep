#' Guess Gene ID type
#'
#' @param ids
#'
#' @return
#' @export
#'
#' @examples
guessIdType <- function(ids) {
  overlaps <- sapply(allIds, function(typedIds) {
    length(intersect(typedIds, ids))
  })
  return(names(overlaps)[which.max(overlaps)])
}
