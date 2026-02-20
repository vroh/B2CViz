#' @import(png)
#' @import(jpeg)
#' @import(tiff)
#' @import(dplyr)
#' @import(Seurat)
#' @import(sf)
#' @import(RBioFormats)
NULL

#' Load image file
#'
#' @param path Path to the image file used for bin2cell segmentation
#' @return A data frame with image data
#' @export
load_img <- function(path) {

  if(grepl("\\.png", path)) {
    im <- png::readPNG(path)
  } else if(grepl("\\.jpg", path)) {
    im <- jpeg::readJPEG(path)
  } else if(grepl("\\.tiff", path)) {
    im <- tiff::readTIFF(path)
  }
  imrgb <- expand.grid(x = 1:dim(im)[1], y = 1:dim(im)[2])
  imrgb$r <- as.vector(im[,,1])
  imrgb$g <- as.vector(im[,,2])
  imrgb$b <- as.vector(im[,,3])
  imrgb$color <- rgb(imrgb$r, imrgb$g, imrgb$b)
  imrgb
}

#' Helper function to convert spaceranger to bin2cell format
#'
#' @param obj Spaceranger object to convert
#' @param slice Which slice to use
#' @return A bin2cell-like object
#' @export
space2b2c <- function(obj = NULL, slice = NULL) {

  # Populate the data slot
  obj <- SetAssayData(
    object = obj,
    assay = "Spatial.Polygons",
    layer = "data",
    new.data = GetAssayData(obj[["Spatial.Polygons"]], layer = "counts")
  )

  # Create the dim reduction with spatial coordinates
  coords <- obj@images[[paste0(slice, ".polygons")]]$centroids@coords
  rownames(coords) <- obj@images[[paste0(slice, ".polygons")]]$centroids@cells
  colnames(coords) <- c("SPATIAL_1", "SPATIAL_2")

  spatial_dr <- Seurat::CreateDimReducObject(
    embeddings = as.matrix(coords),
    key = "SPATIAL_",
    assay = "Spatial.Polygons"
  )

  obj[["spatial"]] <- spatial_dr
  obj
}

#' Load B2C object
#'
#' @param pre Pre-aggregated Seurat object
#' @param post Post-aggregated Seurat object
#' @param path Path to the image file used for bin2cell segmentation
#' @param data Type of input data, either bin2cell (b2c) or spaceranger
#' @param slice Which slice to use in the object (only if using spaceranger data)
#' @param scale.factor The scale factor of the image used (if using spaceranger data or if using a different image than the original bin2cell input image)
#' @return B2C object (A list containing pre, post, image data and other object information)
#' @export
load_b2c <- function(pre = NULL, post = NULL, path = NULL, data = "b2c", slice = NULL, scale.factor = NULL)
{
  if(data == "spaceranger") {
    if(is.null(slice)) stop("Provide slice name when using spaceranger data")
    if(is.null(scale.factor)) stop("Provide scale factor when using spaceranger data")
    pre <- post
  }
  library(png)
  library(jpeg)
  library(tiff)
  library(dplyr)
  library(shiny)
  library(imager)
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(ggnewscale)
  library(tidyr)
  b2c <- list(pre = pre, post = post, path = path)
  b2c$img <- load_img(path)
  b2c$data <- "b2c"
  if(!is.null(scale.factor)) {
    # Rescale the image coordinates
    b2c$img <- b2c$img %>%
      dplyr::mutate(
        x = x / scale.factor,
        y = y / scale.factor
      )
  }
  if(data == "spaceranger") {
    b2c$data <- "spaceranger"
    b2c$slice <- slice

    # convert spaceranger to b2c
    b2c$pre <- space2b2c(post, slice)
    b2c$post <- space2b2c(post, slice)
  }
  if(is.null(scale.factor)) {
    b2c$scale.factor <- 1
  } else {
    b2c$scale.factor <- scale.factor
  }
  b2c
}

#' Scale down original image
#'
#' @param b2c B2C object
#' @param grid.size grid size defining the resolution of the scaled-down image
#' @return B2C object
#' @export
scaledown_img <- function(b2c, grid.size = 10) {
  # Define the grid size
  grid.size <- grid.size

  # Aggregate the data
  output <- b2c$img %>%
    dplyr::mutate(
      x_bin = floor(x / grid.size) * grid.size,
      y_bin = floor(y / grid.size) * grid.size
    ) %>%
    dplyr::group_by(x_bin, y_bin) %>%
    dplyr::summarise(
      r = mean(r),
      g = mean(g),
      b = mean(b)
    ) %>%
    dplyr::select(x_bin, y_bin, r, g, b) %>%
    dplyr::rename(x = "x_bin", y = "y_bin") %>%
    dplyr::mutate(color = rgb(r, g, b))
  b2c$img_sd <- output
  b2c
}

#' Upscale ROI image using original OME-TIFF image
#'
#' @param b2c B2C (cropped) object
#' @param path Path to the original OME-TIFF image
#' @param series Index of the selected OME-TIFF image series
#' @param resolution Index of the selected OME-TIFF image resolution
#' @return B2C (cropped) object
#' @export
upscale_roi <- function(b2c, path, series = 1, resolution = 1) {
  if(!requireNamespace("RBioFormats", quietly = TRUE)) stop("RBioFormats is required but not installed.")
  library(RBioFormats)
  print(read.metadata(path))
  roi <- read.image(
    path,
    series = series,
    resolution = resolution,
    subset = list(X = seq(from = min(b2c$img$y), to = max(b2c$img$y)),
                  Y = seq(from = min(b2c$img$x), to = max(b2c$img$x)))
  )
  roirgb <- data.frame(x = rep(seq(from = min(b2c$img$x), to = max(b2c$img$x)), each = diff(range(b2c$img$y))+1),
                       y = rep(seq(from = min(b2c$img$y), to = max(b2c$img$y)), diff(range(b2c$img$x))+1))
  roirgb$r <- as.vector(roi@.Data[, , 1])
  roirgb$g <- as.vector(roi@.Data[, , 2])
  roirgb$b <- as.vector(roi@.Data[, , 3])
  roirgb$color <- rgb(roirgb$r, roirgb$g, roirgb$b)
  b2c$img <- roirgb
  b2c
}
