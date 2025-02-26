#' @import(png)
#' @import(jpeg)
#' @import(tiff)
#' @import(dplyr)
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

#' Crop B2C object to the dimension of selected ROI
#'
#' @param b2c B2C object
#' @param roi ROI id to crop from
#' @param label.id Label ID of the bin2cell segmentation feature
#' @return B2C object
#' @export
crop_b2c <- function(b2c, roi = 1, label.id = "labels_he_expanded") {
  cells_to_keep <- FetchData(b2c$post, c("SPATIAL_1", "SPATIAL_2"))
  cells_to_keep <- filter(cells_to_keep,
                          SPATIAL_1 >= b2c$coord[[roi]]$xmin,
                          SPATIAL_1 <= b2c$coord[[roi]]$xmax,
                          SPATIAL_2 >= b2c$coord[[roi]]$ymin,
                          SPATIAL_2 <= b2c$coord[[roi]]$ymax)
  b2c$post <- b2c$post[,colnames(b2c$post) %in% rownames(cells_to_keep)]
  b2c$pre <- b2c$pre[,as.character(b2c$pre@meta.data[,label.id]) %in% rownames(cells_to_keep)]
  b2c$img <- filter(b2c$img,
                    y >= b2c$coord[[roi]]$xmin,
                    y <= b2c$coord[[roi]]$xmax,
                    x >= b2c$coord[[roi]]$ymin,
                    x <= b2c$coord[[roi]]$ymax)
  b2c
}
