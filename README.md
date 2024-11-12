# B2CViz

B2CViz is an R package for bin2cell processed VisiumHD spatial single-cell data, including ROI selection and multi-feature plotting.

## Installation

You can install the development version of B2CViz from Gitlab with:

``` r
# install.packages("devtools")
devtools::install_gitlab("vroh/B2CViz")
```

B2CViz depends on the following libraries: 'circlize', 'jpeg', 'png', 'FNN', 'Seurat', 'dplyr', 'ggplot2', 'ggrepel', 'imager', 'shiny'

## Preprocessing

B2CViz requires objects that have been processed by [Bin2cell](https://github.com/Teichlab/bin2cell). Two objects are required, the first one corresponds to the pre-aggregated object (generated after the `b2c.expand_labels` step), while the second is the final aggregated object (after `b2c.bin_to_cell`). Both anndata objects are converted using [sceasy](https://github.com/cellgeni/sceasy) and the following function:

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

To visualize data, load the objects, select a region of interest and call the plotting function:

``` r
library(B2CViz)

# create bin2cell object (works with jpg or png, provide path to image used for bin2cell segmentation)
b2c <- load_b2c(pre = object_pre, post = object_post, path = image_path)

# set region of interest (choose a small region)
set_roi(b2c = b2c)

# plot features
# points is faster and require "post" object only
# hulls require more computation time and both "pre" and "post" Seurat object
plot_b2c(b2c = b2c,
         features = c("Feature1", "Feature2"),
         intensity = 10,
         he_alpha = 0.4,
         pt_size = 0.5,
         plot.type = c("points", "hulls"),
         label.id = "labels_he_expanded"), # column name of the bin2cell segmentation feature
         outline.hulls = c(7966))
```

More informations are available in the vignette
