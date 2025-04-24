# B2CViz

B2CViz is an R package for bin2cell processed VisiumHD spatial single-cell data, including ROI selection and multi-feature plotting.

## Installation

You can install the development version of B2CViz from Gitlab with:

``` r
install.packages("devtools")
devtools::install_gitlab("vroh/B2CViz")
```

B2CViz depends on the following packages: 'jpeg', 'png', 'tiff', 'Seurat', 'dplyr', 'ggplot2', 'ggrepel', 'imager', 'shiny', 'ggnewscale', 'tidyr'

## Preprocessing

B2CViz requires objects that have been preprocessed by [Bin2cell](https://github.com/Teichlab/bin2cell). Two objects are required, the first one corresponds to the pre-aggregated object (generated after the `b2c.expand_labels` step), while the second is the final aggregated object (after `b2c.bin_to_cell`). Both anndata objects are converted using [sceasy](https://github.com/cellgeni/sceasy) and the following function:

``` r
#devtools::install_github("cellgeni/sceasy")
library(reticulate)
library(sceasy)
obj <- convertFormat(obj = "/path/to/obj.h5ad",
                     from="anndata",
                     to="seurat",
                     outFile='/path/to/obj.rds')

```

## Visualization

To visualize data, load the objects, select a region of interest, crop it, and call the plotting functions

``` r
object_pre <- readRDS("/path/to/your/seurat_object.rds")
object_post <- readRDS("/path/to/your/seurat_object.rds")
image_path <- "/path/to/your/image.jpg"

# create bin2cell object (works with jpg or png, provide path to image used for bin2cell)
b2c <- load_b2c(pre = object_pre, post = object_post, path = image_path)

# plot overview
b2c <- scaledown_img(b2c = b2c)
overview_b2c(b2c = b2c, feat = "Cdh1")

# set region of interest
b2c <- set_roi(b2c = b2c)

# crop b2c object using ROI limits
b2c_1 <- crop_b2c(b2c, roi = 1)
b2c_2 <- crop_b2c(b2c, roi = 2)
```

### Default plot

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1")
plot_b2c(b2c = b2c_2, feat = "Cdh1")
```

### Adjust cells transparency

Set alpha.low to a value above 0 to show all cells with positive feature counts

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1", alpha.low = 0.1)
```

### Adjust color gradient scale

Provide a list of vectors to adjust the minimum and maximum values of each color scale

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1", scale.min.max = list(c(0,3)))
```

### Threshold for cells displayed

Adjust min.visible to only show cells that have feature counts above the desired threshold

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1", alpha.low = 0.1, min.visible = 2)
```

### H&E visibility adjustments

Adjust the visibility of the displayed H&E picture with he_alpha

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1", alpha.low = 0.1, min.visible = 2, he_alpha = 0.1)
```

### Cells display format

Choose between points, hulls or a combination of the two

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1", plot.type = "points")
```

### Gradient representation

Differentiate level of counts with a gradient of 2 colors instead of transparency

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1", col.low = "lightblue", col.mid = "white", col.high = "orangered", alpha.low = 1, alpha.mid = 1)
```

### Labels

Show labels if you need to identify cells of interest (slow if too many cells are displayed, so keep the number of cells low!)

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1", min.visible = 6, show.labels = T)
```

### Highlight cells of interest

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1", outline.hulls = c(149249, 186039))
```

### Prefilter to only keep cells of interest

You can select one feature and a threshold expression level to be used for data pre-filtering

``` r
plot_b2c(b2c = b2c_1, feat = "Cdh1", filter.feat = "Clu", filter.threshold = 2)
```

### Multiple features

You can plot multiple features, in that case provide a set of colors (matching the number of features)

``` r
plot_b2c(b2c = b2c_1, feat = c("Cdh1", "Ros1"), col.high = c("orangered", "seagreen2"))
```

min.visible, alpha.low, alpha.mid and alpha.high can be provided as vector when plotting multiple features to adjust parameters for each feature independently

``` r
plot_b2c(b2c = b2c_1, feat = c("Cdh1", "Ros1"), min.visible = c(4, 0), col.high = c("orangered", "seagreen2"))
```

## Quantification

Save the plot in a variable to compute cell-cell distances

``` r
p <- plot_b2c(b2c = b2c_1, feat = c("Cdh1", "Cd4"), plot = F)
output_data <- get_dist(p)
```

This generates a data frame containing all 2 by 2 distances, ids, and feature expression levels for the cells represented in the plot.  
You can then plot the distribution of the distances with plot_dist()

``` r
plot_dist(output_data)
```
