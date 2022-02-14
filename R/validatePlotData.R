ALLOWED_FIELDS_PLOT_DATA = c("fields", "data", "annotations")
ALLOWED_FIELD_TYPES = c("numeric", "factor")
ALLOWED_ANNOTATION_TYPES = c("text", "arrow", "polygon")


validatePlotData <- function(plotData) {

  if (!("fields" %in% names(plotData)))
    stop("Error: Plot data file doesn't have `fields`")
  if (!("data" %in% names(plotData)))
    stop("Error: Plot data file doesn't have `data`")

  sapply(names(plotData), function(x) {
    if (!(x %in% ALLOWED_FIELDS_PLOT_DATA))
      stop(sprintf("Error: Plot data file has a field `%s`. This field is not from allowed list: %s", x, paste(ALLOWED_FIELDS_PLOT_DATA, collapse = ", ")))
  })

  allFields <- names(plotData$fields)
  sapply(allFields, function(field) {

    if (!("type" %in% names(plotData$fields[[field]])))
      stop(sprintf("Error: Plot data field `%s` doesn't have `type", field))
    if (!(plotData$fields[[field]]$type %in% ALLOWED_FIELD_TYPES))
      stop(sprintf("Error: Plot data field `%s` is of unsupported type `%s`. Supported types: %s",
           field, plotData$fields[[field]]$type, paste(ALLOWED_FIELD_TYPES, collapse = ", ")))
  })

  numericFields <- Filter(function(x) plotData$fields[[x]]$type == "numeric", allFields)
  factorFields <- Filter(function(x) plotData$fields[[x]]$type == "factor", allFields)
  isNumericField <- function(field) field %in% numericFields
  isFactorField <- function(field) field %in% factorFields


  checkNumericField <- function(dataPoint, field) {

    if ((dataPoint[[field]] < plotData$fields[[field]]$range[1]) ||
        (plotData$fields[[field]]$range[2] < dataPoint[[field]]))
      stop(sprintf("
Error: numeric value is out of range
Expected range for `%s` is [%s]
Actual value: %s
Whole data point:
%s", field, paste(plotData$fields[[field]]$range, collapse=", "),
     dataPoint[[field]], prettify(toJSON(dataPoint, auto_unbox = T))))
  }

  checkFactorField <- function(dataPoint, field) {

    if (!(dataPoint[[field]] %in% plotData$fields[[field]]$levels)) {
      stop(sprintf("
Error: factor value is not from the levels
Levels for `%s` are [%s]
Actual value: %s
Whole data point:
%s", field, paste(plotData$fields[[field]]$fields, collapse=", "),
     dataPoint[[field]], prettify(toJSON(dataPoint, auto_unbox = T))))
    }
  }

  checkField <- function(dataPoint, field) {
    if (isNumericField(field)) checkNumericField(dataPoint, field)
    if (isFactorField(field)) checkFactorField(dataPoint, field)
  }

  sapply(numericFields, function(field) {
    if (!("range" %in% names(plotData$fields[[field]])))
      stop(sprintf("Error: Numeric field `field` doesn't have `range`", field))
  })

  sapply(factorFields, function(field) {
    if (!("levels" %in% names(plotData$fields[[field]])))
      stop(sprintf("Error: Factors field `field` doesn't have `levels`", field))
  })



  sapply(plotData$data, function(point) {
    sapply(allFields, function(field) {

      if (!(field %in% names(point)))
        stop(sprintf("Error: field `%s` is not present in the data point
Field `%s` is described in plot_data.json
Data point:
%s", field, field, prettify(toJSON(point, auto_unbox = T))))

      checkField(point, field)
    })
  })

  if ("annotations" %in% names(plotData)) {

    annotations <- names(plotData$annotations)
    for (annotation in annotations) {

      annObject <- plotData$annotations[[annotation]]
      if (!(annObject$type %in% ALLOWED_ANNOTATION_TYPES))
        stop(sprintf("Error: Plot data file has an annotation of unsupported type `%s`.
Supported annotation type are: %s", annObject$type, paste(ALLOWED_ANNOTATION_TYPES, collapse = ", ")))

      # checking coordinates

      if (!all(annObject$coords %in% numericFields))
        stop(sprintf("Error: Annotation coordinates must be from numeric fields.
For annotation `%s`, coordinates found: %s", annotation, paste(annObject$coords, collapse = ", ")))

      if (annObject$type == "text") {
        if (!("data" %in% names(annObject)))
          stop(sprintf("Error: Annotation `%s` missing data field.", annotation))
        if (!("value" %in% names(annObject)))
          stop(sprintf("Error: Annotation `%s` missing value field.", annotation))

        subFields <- c(annObject$coords[[1]], annObject$coords[[2]], annObject$value)

        for (dataPoint in annObject$data) {
          sapply(subFields, function(subField) {
            if (!(subField %in% names(dataPoint)))
              stop(sprintf("
Error: Text annotation doesn't have `%s`
Coordinates are: [%s]
Value is %s:
Whole data point:
%s", subField, paste(annObject$coords, collapse=", "),
                           annObject$value, prettify(toJSON(dataPoint, auto_unbox = T))))
          })
        }


      }

      if (annObject$type == "arrows") {
        if (!("data_start" %in% names(annObject)))
          stop(sprintf("Error: Annotation `%s` missing data_start field.", annotation))
        if (!("data_end" %in% names(annObject)))
          stop(sprintf("Error: Annotation `%s` missing data_end field.", annotation))


        data_start <- annObject$data_start
        data_end <- annObject$data_end

        if (!(length(data_start) == length(data_end))) {
          stop(sprintf("Error: Arrow annotation `%s` has different lengths of data_start and data_end:
                       data_start: %d, data_end: %d",
                       annotation, length(data_start), length(data_end)))
        }

        subFields <- c(annObject$coords[[1]], annObject$coords[[2]])

        for (dataPoint in annObject$data_start) {
          sapply(subFields, function(subField) {
            if (!(subField %in% names(dataPoint)))
              stop(sprintf("
Error: Arrow annotation (data_start) doesn't have `%s`
Coordinates are: [%s]
Whole data point:
%s", subField, paste(annObject$coords, collapse=", "), prettify(toJSON(dataPoint, auto_unbox = T))))
          })
        }

        for (dataPoint in annObject$data_end) {
          sapply(subFields, function(subField) {
            if (!(subField %in% names(dataPoint)))
              stop(sprintf("
Error: Arrow annotation (data_end) doesn't have `%s`
Coordinates are: [%s]
Whole data point:
%s", subField, paste(annObject$coords, collapse=", "), prettify(toJSON(dataPoint, auto_unbox = T))))
          })
        }

      }

      if (annObject$type == "polygon") {
        if (!("data" %in% names(annObject)))
          stop(sprintf("Error: Annotation `%s` missing data field.", annotation))
        if (!("value" %in% names(annObject)))
          stop(sprintf("Error: Annotation `%s` missing value field.", annotation))
      }

      subFields <- c(annObject$coords[[1]], annObject$coords[[2]], annObject$value)

      for (dataPoint in annObject$data) {
        sapply(subFields, function(subField) {
          if (!(subField %in% names(dataPoint)))
            stop(sprintf("
Error: Polygon annotation doesn't have `%s`
Coordinates are: [%s]
Value is %s:
Whole data point:
%s", subField, paste(annObject$coords, collapse=", "),
                         annObject$value, prettify(toJSON(dataPoint, auto_unbox = T))))
        })
      }

    }
  }

  message("Plot data looks OK")
  message(sprintf("Plot data consists of %d data points", length(plotData$data)))
  message("Fields:")
  message(sprintf("Numeric fields: [%s]", paste(numericFields, collapse=", ")))
  message(sprintf("Factor fields: [%s]", paste(factorFields, collapse=", ")))

  if ("annotations" %in% names(plotData)) {
    message("Available annotations:")
    annotations <- names(plotData$annotations)
    # for (annotation in annotations) {
    #   annObject <- plotData$annotations[[annotation]]
    #   message(sprintf("Annotation `%s`, type: `%s`, coordinates: [%s]",
    #                   annotation, annObject$type, paste(annObject$coords, collapse=", ")))
    # }
    message(paste(annotations, collapse=", "))
  }

  return(invisible(NULL))
}
