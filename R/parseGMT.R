
## parseGMT


parseGMTFile <- function(path) {
  lines <- readLines(path)
  lines <- strsplit(lines, '\t')

  sapply(lines, function(x) {
    if (length(x) != 3) {
      stop(sprintf("Error: gmt line doesn't have 3 elements The line: %s", x))
    }
  })

  gmtNames <- sapply(lines, function(x) x[1])
  universeNames <- sapply(lines, function(x) x[2])
  genes <- sapply(lines, function(x) strsplit(x[3], ",")[[1]])

  return(list(
    geneSetNames=gmtNames,
    universeNames=universeNames,
    genes=genes
  ))

}
