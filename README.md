# B2CViz

B2CViz is an R package for bin2cell or Spaceranger processed VisiumHD spatial single-cell data, including ROI selection, multi-feature plotting and quantification. Figures presented in this page were generated with Bin2cell and B2CViz using data from [10x genomics](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-human-breast-cancer-ff-ultima)

![Original 2 microns bins, aggregated 8 microns bins, bin2cell centroids and B2CViz cells](man/figures/B2CViz.jpg)

## Installation

You can install the development version of B2CViz from Gitlab with:

``` r
install.packages("devtools")
devtools::install_gitlab("vroh/B2CViz")
```

B2CViz depends on the following packages: 'jpeg', 'png', 'tiff', 'Seurat', 'dplyr', 'ggplot2', 'ggrepel', 'imager', 'shiny', 'ggnewscale', 'tidyr', 'sf'

## Preprocessing

Initial versions of B2CViz required objects that have been preprocessed by [Bin2cell](https://github.com/Teichlab/bin2cell). Two objects are required, the first one corresponds to the pre-aggregated object (generated after the `b2c.expand_labels` step), while the second is the final aggregated object (after `b2c.bin_to_cell`). Both anndata objects are converted using [sceasy](https://github.com/cellgeni/sceasy) and the following function:

``` r
#devtools::install_github("cellgeni/sceasy")
library(reticulate)
library(sceasy)
obj <- convertFormat(obj = "/path/to/obj.h5ad",
                     from="anndata",
                     to="seurat",
                     outFile='/path/to/obj.rds')
```

Since January 2026, objects that have been segmented directly with Spaceranger can also be used (see how to use them below).

## Visualization

To visualize data, load the objects, select a region of interest, crop it, and call the plotting functions

``` r
library(B2CViz)
library(Seurat)
object_pre <- readRDS("/path/to/seurat_object.rds")
object_post <- readRDS("/path/to/seurat_object.rds")
image_path <- "/path/to/image.jpg"

# Normalize bin2cell data
object_post <- NormalizeData(object_post)

# create bin2cell object (works with jpg or png, provide path to image used for bin2cell)
b2c <- load_b2c(pre = object_pre, post = object_post, path = image_path)
```

If you are using data directly segmented by Spaceranger, you can load them as follows

``` r
# Specify the path to your Spaceranger outputs
data.dir <- "/path/to/outs/"
# Choose which data slice to use
slice <- "slice1"

# Load the Visium HD dataset with segmentations
obj <- Load10X_Spatial(
  data.dir = data.dir,
  slice = slice,
  bin.size = c("polygons")
)

# The scaling factor is required. For example, when using spaceranger's hires scale it is reported here:
sf <- obj@images[[paste0(slice, ".polygons")]]@scale.factors$hires

# Create a bin2cell (B2C) object
image_path <- "/path/to/outs/segmented_outputs/spatial/tissue_hires_image.png"

b2c <- load_b2c(post = obj, path = image_path, data = "spaceranger", slice = slice, scale.factor = sf)

# normalize with NormalizeData or SCTransform (though latter probably better)
b2c$post <- NormalizeData(b2c$post)
```

You can use Spaceranger's hires image as input, though I would recommend to generate a higher definition image. If you have massive OME-TIFF files, you can prepare jpgs of manageable size using [OT-LIP](https://gitlab.com/vroh/ot-lip). Record the scale factor used, it is required by the `load_b2c` function.

### Overview

``` r
# plot overview
b2c <- scaledown_img(b2c = b2c)
overview_b2c(b2c = b2c, feat = "CDH1")
```

![Overview plot](man/figures/overview.jpg)

### ROI selection

Select your region of interest using `set_roi()` function. Draw rectangles or polygons and add the ROIs to the b2c object

``` r
# set region of interest
b2c <- set_roi(b2c = b2c)

# crop b2c object using ROI limits
b2c_1 <- crop_b2c(b2c, roi = 1)
b2c_2 <- crop_b2c(b2c, roi = 2)
```

![ROI selection](man/figures/roi_selector.jpg)

### Segmentation plot

You can inspect the segmentation results by running

``` r
plot_segmentation(b2c = b2c_1)
```

### Feature plot

The default feature plot function is `plot_b2c`

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1")
plot_b2c(b2c = b2c_2, feat = "CDH1")
```

![Default plot](man/figures/Slide1.JPG)

### Adjust cells transparency

Set alpha.low to a value above 0 to show all cells with positive feature counts

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1", alpha.low = 0.1)
```

### Adjust color gradient scale

Provide a list of vectors to adjust the minimum and maximum values of each color scale

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1", scale.min.max = list(c(0,3)))
```

### Threshold for cells displayed

Adjust min.visible to only show cells that have feature counts above the desired threshold

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1", alpha.low = 0.1, min.visible = 3)
```

![Higher minimum threshold for visible cells](man/figures/Slide2.JPG)

### H&E visibility adjustments

Adjust the visibility of the displayed H&E picture with he_alpha

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1", alpha.low = 0.1, min.visible = 2, he_alpha = 0.1)
```

![Higher transparency for the H&E image](man/figures/Slide3.JPG)

### Cells display format

Choose between points (bin2cell centroids), hulls (cells) or a combination of the two

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1", plot.type = "points")
```

![Displaying centroids only](man/figures/Slide4.JPG)

### Gradient representation

Differentiate level of counts with a gradient of 2 colors instead of transparency

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1", col.low = "blue", col.mid = "white", col.high = "red", alpha.low = 1, alpha.mid = 1)
```

![2-color scale](man/figures/Slide5.JPG)

### Labels

Show labels if you need to identify cells of interest (slow and overcrowded if too many cells are displayed, so keep the number of cells low!)

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1", min.visible = 3.5, show.labels = T)
```

![Cell ID identification](man/figures/Slide6.JPG)

### Highlight cells of interest

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1", outline.hulls = c(154900, 157962, 157507))
```

![Highlighting cells](man/figures/Slide7.JPG)

### Pre-filter to only keep cells of interest

You can select features and threshold expression levels to be used for data pre-filtering. If using multiple features, add them in lists of vectors in the same order than the feat parameter. Filters for multiple features are combined according to `filter.type` which can be either in `or` mode (default) or `ànd` mode.

``` r
plot_b2c(b2c = b2c_1, feat = "CDH1", filter.feat = "CLU", filter.threshold = 5)
plot_b2c(b2c = b2c_1, feat = "CDH1", filter.feat = list(c("CLU", "EPCAM")), filter.threshold = list(c(5, 2)))
```

![Pre-filtering cells for high CLU expression](man/figures/Slide8.JPG)

### Multiple features

You can plot multiple features, in that case provide a set of colors (matching the number of features)

``` r
plot_b2c(b2c_1, c("CDH1", "CD19", "COL1A1", "MYH11"), plot.type = "hulls", col.high = c("dodgerblue", "orangered", "gold", "seagreen"), he_alpha = 0.3)
```

![Multi-feature plot](man/figures/Slide9.JPG)

If multiple features overlap a lot, it can be better to use centroids representation only to identify multi-positive cells more easily

``` r
plot_b2c(b2c_1, c("CD3E", "CD4"), plot.type = "points", col.high = c("dodgerblue", "orangered"))
```

![Overlapping features](man/figures/multi_overlap.jpg)

min.visible, alpha.low, alpha.mid and alpha.high can be provided as vector when plotting multiple features to adjust parameters for each feature independently

``` r
plot_b2c(b2c = b2c_1, feat = c("CDH1", "CD4"), min.visible = c(2, 2.5), col.high = c("orangered", "seagreen2"))
```

## Quantification

Save the plot in a variable to compute cell-cell distances

![Multi-feature plot](man/figures/Slide10.JPG)

``` r
p <- plot_b2c(b2c = b2c_1, feat = c("CDH1", "CD4"), min.visible = c(2, 2.5), plot = F)
output_data <- get_dist(p, radius = 1000)
```

This generates a data frame containing all 2 by 2 distances, ids, and feature expression levels for the cells displayed in the plot.\
You can then plot the distribution of the distances with `plot_dist()`

``` r
plot_dist(output_data, binwidth = 20)
```

![Distribution of distances](man/figures/Slide11.JPG)
