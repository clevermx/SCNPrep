
ALLOWED_EXP_TYPE_VALUES <- c("counts", "as_is")

validateExpressionData <- function(expDataJson, h5path) {

  x <- as.numeric(h5read(h5path, "X/data"))
  p <- as.numeric(h5read(h5path, "X/indptr"))
  j <- as.numeric(h5read(h5path, "X/indices"))

  id <- H5Fopen(h5path)
  group <- H5Gopen(id, "X")
  shapeId <- H5Aopen(group, "shape")
  shape <- H5Aread(shapeId)
  h5closeAll()

  counts <- sparseMatrix(j = j, p = p, x = x, repr="R", dims=c(shape[2], shape[1]), index1=F)


  if (!("features" %in% names(expDataJson)))
    stop("Error: dataset expression data descriptor doesn't have `features`")

  if (!("barcodes" %in% names(expDataJson)))
    stop("Error: dataset expression data descriptor doesn't have `barcodes`")

  if (!("featureCounts" %in% names(expDataJson)))
    stop("Error: dataset expression data descriptor doesn't have `featureCounts`")

  if (!("totalCounts" %in% names(expDataJson)))
    stop("Error: dataset expression data descriptor doesn't have `totalCounts`")

  if (!("expType" %in% names(expDataJson)))
    stop("Error: dataset expression data descriptor doesn't have `expType`")

  if (!(expDataJson$expType %in% ALLOWED_EXP_TYPE_VALUES))
    stop(sprintf("Error: dataset expression data descriptor has unsupported `expType`: %s.
Supported values are: %s", expDataJson$expType, paste(ALLOWED_EXP_TYPE_VALUES, collapse = ", ")))


  genes <- expDataJson$features
  barcodes <- expDataJson$barcodes
  shapeExp <- c(length(barcodes), length(genes))

  if (length(barcodes) != length(expDataJson$totalCounts)) {
    stop(sprintf("Error: dataset expression data descriptor has barcodes and totalCounts of different lengths: %d and %d respectively",
                 length(barcodes), length(expDataJson$totalCounts)))
  }

  if (length(genes) != length(expDataJson$featureCounts)) {
    stop(sprintf("Error: dataset expression data descriptor has `features` and `featureCounts` of different lengths: %d and %d respectively",
                 length(genes), length(expDataJson$featureCounts)))
  }


  shapeExp <- c(length(barcodes), length(genes))

  if (!all(shape==shapeExp)) {
    stop(sprintf("Error: dataset exp_data.json and data.h5 have different sizes.
exp_data.json: [%s]
data.h5: [%s]", paste(shapeExp, collapse = ", "), paste(shape, collapse = ", ")))
  }

  message("Expression data looks OK")
  message(sprintf("Expression data consists of %d features and %d barcodes", length(genes), length(barcodes)))

  return(invisible(NULL))
}
