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
  library(shiny)
  library(imager)
  library(ggplot2)
  library(dplyr)
  library(FNN)

  if(grepl("\\.png", path)) {
    im <- png::readPNG(path)
  } else if(grepl("\\.jpg", path)) {
    im <- jpeg::readJPEG(path)
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
