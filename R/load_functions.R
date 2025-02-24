#' Load image file
#'
#' @param path Path to the image file used for bin2cell segmentation
#' @return A data frame with image data
#' @export
load_img <- function(path) {

  # load all required packages
  library(Seurat)
  library(png)
  library(jpeg)
  library(tiff)
  library(shiny)
  library(imager)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(FNN)

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

#' Load B2C object
#'
#' @param pre Pre-aggregated Seurat object
#' @param post Post-aggregated Seurat object
#' @param path Path to the image file used for bin2cell segmentation
#' @return B2C object (A list containing pre, post, and image data)
#' @export
load_b2c <- function(pre = NULL, post = NULL, path = NULL) {
  b2c <- list(pre = pre, post = post, path = path)
  b2c$img <- load_img(path)
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
    mutate(
      x_bin = floor(x / grid.size) * grid.size,
      y_bin = floor(y / grid.size) * grid.size
    ) %>%
    group_by(x_bin, y_bin) %>%
    summarise(
      r = mean(r),
      g = mean(g),
      b = mean(b)
    ) %>%
    select(x_bin, y_bin, r, g, b) %>%
    rename(x = "x_bin", y = "y_bin") %>%
    mutate(color = rgb(r, g, b))
  b2c$img_sd <- output
  b2c
}

#' Load ENACT object
#'
#' @param poly ENACT table (csv) containing polygon data
#' @param post Post-aggregated (with ENACT) Seurat object
#' @param path Path to the image file used for ENACT segmentation
#' @return B2C object (A list containing pre, post, and image data, adapted from ENACT objects)
#' @export
load_enact <- function(poly = NULL, post = NULL, path = NULL) {



  b2c <- list(pre = pre, post = post, path = path)
  b2c$img <- load_img(path)
  b2c
}
