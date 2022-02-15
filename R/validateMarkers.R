
#' @Serializable
#' data class MarkerEntry(
#'   @SerialName("p_val")
#'   val pValue: Double,
#'
#'   @SerialName("p_val_adj")
#'   val pValueAdjusted: Double,
#'
#'   @SerialName("avg_logFC")
#'   val averageLogFoldChange: Double,
#'
#'   @SerialName("pct.1")
#'   val pct1: Double,
#'
#'   @SerialName("pct.2")
#'   val pct2: Double,
#'
#'   @SerialName("cluster")
#'   val cluster: String,
#'
#'   @SerialName("gene")
#'   val gene: String
#' )

REQUIRED_FIELDS_MARKERS <- c("gene", "cluster", "avg_logFC", "p_val", "p_val_adj", "pct.1", "pct.2")
REQUIRED_MARKER_TYPE_CHECKS <- list(
  "gene"=is.character,
  "cluster"=is.character,
  "avg_logFC"=is.numeric,
  "p_val"=is.numeric,
  "p_val_adj"=is.numeric,
  "pct.1"=is.numeric,
  "pct.2"=is.numeric
)

REQUIRED_MARKER_TYPES <- list(
  "gene"="character",
  "cluster"="character",
  "avg_logFC"="numeric",
  "p_val"="numeric",
  "p_val_adj"="numeric",
  "pct.1"="numeric",
  "pct.2"="numeric"
)

validateMarkersJson <- function(markersJson) {
  tableNames <- names(markersJson)

  for (table in tableNames) {

    markers <- markersJson[[table]]

    for (marker in markers) {
      sapply(REQUIRED_FIELDS_MARKERS, function(subField) {
        if (!(subField %in% names(marker)))
          stop(sprintf("
Error: A marker in table `%s` don't have required field `%s`
Marker data:
%s", table, subField, prettify(toJSON(marker, auto_unbox = T))))


        if(!(REQUIRED_MARKER_TYPE_CHECKS[[subField]](marker[[subField]]))) {
          stop(sprintf("
Error: In table %s, a marker has `%s` field of unsupported type. Expected type is `%s`
Marker data:
%s", table, subField, REQUIRED_MARKER_TYPES[[subField]],
                       prettify(toJSON(marker, auto_unbox = T))))
        }

      })
    }

    message(sprintf("Marker table `%s` looks OK and contains %d markers", table, length(markers)))

  }

}
