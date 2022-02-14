
# required field GMT annotation
REQUIRED_FIELDS_GMT_ANNOTATION <- c("token", "name", "title", "species", "link")

validateGMTFiles <- function(gmtInfo, gmtAnnotation) {

  for (i in 1:length(gmtInfo$geneSetNames)) {
    geneSet <- gmtInfo$geneSetNames[i]
    genes <- gmtInfo$genes[i]
    if (!(geneSet %in% names(gmtAnnotation))) {
      stop(sprintf("Error: GMT gene set `%s` is described in GMT file but missing in JSON annotation", geneSet))
    }

    if (length(genes) == 0) {
      stop(sprintf("Error: GMT gene set `%s` is empty.", geneSet))
    }

  }

  for (geneSet in names(gmtAnnotation)) {
    if (!(geneSet %in% gmtInfo$geneSetNames)) {
      stop(sprintf("Error: GMT gene set `%s` is described in JSON annotation but missing in GMT file", geneSet))
    }
  }

  for (gmtDescriptor in gmtAnnotation) {
    sapply(REQUIRED_FIELDS_GMT_ANNOTATION, function(subField) {
      if (!(subField %in% names(gmtDescriptor)))
        stop(sprintf("
Error: GMT Json Descriptor is missing required field `%s`
GMT Json Descriptor:
%s", subField, prettify(toJSON(gmtDescriptor, auto_unbox = T))))
    })
  }

  message(sprintf('GMT files look OK. Total of %d gene sets are present.', length(gmtAnnotation)))
}
