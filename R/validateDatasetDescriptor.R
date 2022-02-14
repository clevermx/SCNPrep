# validate dataset.descriptor
#
# data class SCdescriptor (
#   val token: String,
#   val name: String? = null,
#   val description: String? = null,
#   val link: String? = null,
#   val species: Species,
#   val cells: Int = 0,
#   val public: Boolean = false,
#   val curated: Boolean = false,
#   val debug: Boolean = false
# )

ALLOWED_FIELDS_DATASET_DESCRIPTOR = c(
  "token", "name", "description", "link",
  "species", "cells", "public", "curated", "debug"
)


#' validateDatasetDescriptor
#'
#' @param descriptor
#'
#' @return Nothing if all checks are correct, stops otherwise
#'
#' @examples
validateDatasetDescriptor <- function(descriptor) {
  if (!("token" %in% names(descriptor)))
    stop("Error: dataset descriptor doesn't have `token`")
  if (!("species" %in% names(descriptor)))
    stop("Error: dataset descriptor doesn't have `species`")
  if (!("cells" %in% names(descriptor)))
    stop("Error: dataset descriptor doesn't have `cells`")
  if (!(descriptor$species %in% ALLOWED_SPECIES))
    stop(sprintf("Error: dataset descriptor species is not from the list: %s", paste(ALLOWED_SPECIES, collapse = ", ")))

  sapply(names(descriptor), function(x) {
    if (!(x %in% ALLOWED_FIELDS_DATASET_DESCRIPTOR))
      stop(sprintf("Error: dataset descriptor has a field `%s`. This field is not from allowed list: %s", x, paste(ALLOWED_FIELDS_DATASET_DESCRIPTOR, collapse = ", ")))
  })


  message("Dataset descriptor looks OK")
  message(prettify(toJSON(descriptor, auto_unbox = T)))

  return(invisible(NULL))
}
