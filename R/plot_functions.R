#' @import(Seurat)
#' @import(ggplot2)
#' @import(ggrepel)
#' @import(dplyr)
#' @import(ggnewscale)
NULL

#' Plot spatial feature from a B2C object (requires scaled-down image)
#'
#' @param b2c B2C object
#' @param feat Gene feature to plot
#' @param pt.size Size of the points
#' @param he_alpha Alpha value for H&E image
#' @param col.low Low color for colorscale
#' @param col.high High color for colorscale
#' @param scalebar scalebar size in microns, FALSE for no scalebar
#' @export
overview_b2c <- function(b2c, feat, pt.size = 0.001, he_alpha = 0.4, col.low = "gray", col.high = "seagreen2", scalebar = 500) {

  if (is.null(b2c$img_sd)) {
    stop("run scaledown_img() on b2c object first")
  }

  df <- FeaturePlot(b2c$post, feat, reduction = "spatial")[[1]]$data
  colnames(df)[4] <- "feat"

  p <-
    ggplot(filter(df, feat > 0), aes(y, x)) +
    geom_raster(data = b2c$img_sd, aes(fill = color), alpha = he_alpha) +
    geom_point(aes(x = SPATIAL_1, y = SPATIAL_2, col = feat),
               alpha = scales::rescale(filter(df, feat > 0)$feat),
               size = pt.size) +
    coord_fixed(ratio = 1) +
    scale_fill_identity() +
    scale_color_continuous(name = feat, low = col.low, high = col.high) +
    xlim(min(df$SPATIAL_1), max(df$SPATIAL_1)) +
    ylim(max(df$SPATIAL_2), min(df$SPATIAL_2)) +
    theme_void() +
    ggtitle(feat) +
    theme(plot.title = element_text(hjust = 0.5))

  # scalebar
  if(scalebar) {
    meta <- cbind(b2c$pre@meta.data, b2c$pre@reductions$spatial@cell.embeddings)

    # Calculate micron-to-plot conversion factor
    adjacent_spots <- meta %>%
      filter(array_col == min(array_col)) %>% # Same column
      arrange(array_row) %>%
      slice(1:2) # First two rows in same column

    if(nrow(adjacent_spots) < 2) stop("Insufficient adjacent spots for scale calculation, can't add scalebar, set scalebar_micron = FALSE")

    y_distance <- abs(adjacent_spots$SPATIAL_2[2] - adjacent_spots$SPATIAL_2[1])
    microns_per_bin <- 2 # Your known bin spacing
    plot_units_per_micron <- y_distance / microns_per_bin

    # Convert desired microns to plot units
    scalebar_adj <- scalebar * plot_units_per_micron

    # Calculate scalebar position (bottom-left corner)
    x_range <- range(df$SPATIAL_1)
    y_range <- range(df$SPATIAL_2)
    x_pos <- x_range[1] + 0.02 * diff(x_range)  # 2% from left edge
    y_pos <- y_range[2] - 0.03 * diff(y_range)  # 3% from visual bottom (after flip)

    # Add scalebar to plot
    p <-
      p +
      annotate("rect",
               xmin = x_pos, xmax = x_pos + scalebar_adj,
               ymin = y_pos, ymax = y_pos + diff(y_range) * 0.01,
               fill = "black",
               color = NA
      ) +
      annotate("text",
               x = x_pos + scalebar_adj / 2,
               y = y_pos - diff(y_range) * 0.015, # Adjusted for better visibility
               label = paste(scalebar, "μm"),
               size = 3,
               color = "black"
      )
  }

  print(p)

}

#' Plot B2C object
#'
#' @param b2c B2C object
#' @param feat Features to plot
#' @param label.id Label ID of the bin2cell segmentation feature
#' @param min.visible Feature value threshold for visualization (single value or vector for each feature)
#' @param col.low Low color for colorscale (single value or vector for each feature)
#' @param col.mid Mid color for colorscale (single value or vector for each feature)
#' @param col.high High color for colorscale (single value or vector for each feature)
#' @param alpha.low Low alpha value for colorscale (single value or vector for each feature)
#' @param alpha.mid Mid alpha value for colorscale (single value or vector for each feature)
#' @param alpha.high High alpha value for colorscale (single value or vector for each feature)
#' @param scale.min.max List of vectors (of length 2) indicating the min and max value for the color gradient scale
#' @param pt_size Point size
#' @param shape Shape of the points. Defaults to solid circle (16) for 1 feature, empty circle (21) for multi-feature plot to better see multi-positive cells
#' @param he_alpha Alpha value for H&E image
#' @param title Plot title
#' @param plot.type Type of plot (points or hulls or both)
#' @param show.bins Whether to show the 2um bins data. Defaults to "no", use "yes" to show bins corresponding to the feature and "all" to show all bins data (in white)
#' @param outline.hulls Character vector of hulls to outline (by expanded labels ID)
#' @param show.labels Whether or not to plot the hulls labels
#' @param plot Whether to display or return the plot and the data
#' @param scalebar Scalebar size in microns, FALSE for no scalebar
#' @param scalebar.width width of the scalebar
#' @param translate Whether or not to translate the plot to the (0, 0) origin (can be useful to adjust plot size when comparing multiple ROIs)
#' @param filter.feat Features (list of vectors) to pre-filter the data (e.g. keep only cells expressing these features, order has to match order used in feat, use "" for no filtering)
#' @param filter.threshold Threshold (list of vectors) levels for filter.feat
#' @export
plot_b2c <- function(b2c, feat, label.id = "labels_he_expanded", min.visible = 0,
                     col.low = NULL, col.mid = NULL, col.high = "orangered", alpha.low = 0, alpha.mid = 0.5, alpha.high = 1, scale.min.max = NULL,
                     pt.size = 1, shape = NULL, he_alpha = 0.3, title = NULL, plot.type = c("points", "hulls"), show.bins = "no",
                     outline.hulls = NULL, show.labels = F, plot = T, scalebar = 200,
                     scalebar.width = 10, translate = T, filter.feat = NULL, filter.threshold = 0) {

  # prefilter data?
  if(!is.null(filter.feat)) {
    if(length(filter.threshold) == 1 & is.numeric(filter.threshold)) {
      filter.threshold <- sapply(lapply(filter.feat, length), function(u) list(rep(filter.threshold, u)))
    }
  }

  # adjust parameters
  if(length(min.visible) == 1) {
    min.visible <- rep(min.visible, length(feat))
  }
  if(length(alpha.low) == 1) {
    alpha.low <- rep(alpha.low, length(feat))
  }
  if(length(alpha.mid) == 1) {
    alpha.mid <- rep(alpha.mid, length(feat))
  }
  if(length(alpha.high) == 1) {
    alpha.high <- rep(alpha.high, length(feat))
  }
  if(is.null(col.mid)) {
    col.mid <- col.high
  }
  if(is.null(col.low)) {
    col.low <- col.high
  }
  if(is.null(shape)) {
    shape <- ifelse(length(feat) > 1, 21, 16)
  }

  # fetch data
  df_post <- FetchData(b2c$post, vars = c("SPATIAL_1", "SPATIAL_2", feat, unlist(filter.feat)))
  df_post[label.id] <- row.names(df_post)
  if("hulls" %in% plot.type) {
    df_pre <- FetchData(b2c$pre, vars = c("SPATIAL_1", "SPATIAL_2", label.id)) %>%
      group_by(across(label.id)) %>%
      slice(chull(SPATIAL_1, SPATIAL_2))
    df <- merge(df_post, df_pre, by = label.id)
  }

  # translate to origin
  if(translate) {
    if("hulls" %in% plot.type) {

      translate_sp1 <- min(df$SPATIAL_1.y)
      translate_sp2 <- min(df$SPATIAL_2.y)

      df_post$SPATIAL_1 <- df_post$SPATIAL_1 - translate_sp1
      df_post$SPATIAL_2 <- df_post$SPATIAL_2 - translate_sp2
      df$SPATIAL_1.y <- df$SPATIAL_1.y - translate_sp1
      df$SPATIAL_2.y <- df$SPATIAL_2.y - translate_sp2
      b2c$img$x <- b2c$img$x - translate_sp2
      b2c$img$y <- b2c$img$y - translate_sp1

    } else {

      translate_sp1 <- min(df_post$SPATIAL_1)
      translate_sp2 <- min(df_post$SPATIAL_2)

      df_post$SPATIAL_1 <- df_post$SPATIAL_1 - translate_sp1
      df_post$SPATIAL_2 <- df_post$SPATIAL_2 - translate_sp2
      b2c$img$x <- b2c$img$x - translate_sp2
      b2c$img$y <- b2c$img$y - translate_sp1

    }
  }

  # plot H&E
  p <-
    ggplot(df_post) +
    geom_raster(data = b2c$img, aes(x = y, y = x, fill = color), alpha = he_alpha) +
    scale_fill_identity() +
    ggnewscale::new_scale_fill()

  # show bins
  if(show.bins == "yes" | show.bins == "all") {
    bins <- FetchData(b2c$pre, vars = c("SPATIAL_1", "SPATIAL_2", feat))
    colnames(bins) <- paste0(colnames(bins), "_bins")
    bins$SPATIAL_1 <- bins$SPATIAL_1 - ifelse(translate, translate_sp1, 0)
    bins$SPATIAL_2 <- bins$SPATIAL_2 - ifelse(translate, translate_sp2, 0)
    for(i in 1:length(feat)) {
      p <-
        p +
        geom_point(data = bins, aes(x = SPATIAL_1, y = SPATIAL_2, col = !!sym(paste0(feat[i], "_bins"))), size = 0.3) +
        scale_color_gradient2(mid = ifelse(show.bins == "all", alpha("white", alpha = 0.65/length(feat)), alpha(col.mid[i], alpha = 0)),
                              high = col.high[i],
                              na.value = "transparent") +
        ggnewscale::new_scale_color()
    }
  }

  # plot cells
  plotted <- NULL
  if("points" %in% plot.type & !("hulls" %in% plot.type)) {
    for(i in 1:length(feat)) {
      data.points <- df_post[df_post[[feat[i]]] > min.visible[i], ]
      if(!is.null(filter.feat)) {
        for(j in 1:length(filter.feat[[i]])) {
          if(filter.feat[[i]][j] == "") {
            next
          } else {
            data.points <- data.points[data.points[[filter.feat[[i]][j]]] > filter.threshold[[i]][j], ]
          }
        }
      }
      p <-
        p +
        geom_point(data = data.points,
                   aes_string(x = "SPATIAL_1", y = "SPATIAL_2", col = feat[i]), shape = shape, fill = NA, size = pt.size*i, stroke = 2*pt.size/3) +
        scale_color_gradient2(low = alpha(col.low[i], alpha = alpha.low[i]),
                              mid = alpha(col.mid[i], alpha = alpha.mid[i]),
                              high = alpha(col.high[i], alpha = alpha.high[i]),
                              midpoint = ifelse(is.null(scale.min.max), max(df_post[feat[i]])/2, mean(scale.min.max[[i]])),
                              na.value = "transparent",
                              limits = c(ifelse(is.null(scale.min.max), min(df_post[feat[i]]), scale.min.max[[i]][1]),
                                         ifelse(is.null(scale.min.max), max(df_post[feat[i]]), scale.min.max[[i]][2]))) +
        ggnewscale::new_scale_color()
      plotted <- c(plotted, data.points[,length(data.points)])
    }
  } else if("hulls" %in% plot.type & !("points" %in% plot.type)) {
    for(i in 1:length(feat)) {
      data.points <- df_post[df_post[[feat[i]]] > min.visible[i], ]
      data.hulls <-  df[df[[feat[i]]] > min.visible[i], ]
      if(!is.null(filter.feat)) {
        for(j in 1:length(filter.feat[[i]])) {
          if(filter.feat[[i]][j] == "") {
            next
          } else {
            data.points <- data.points[data.points[[filter.feat[[i]][j]]] > filter.threshold[[i]][j], ]
            data.hulls <- data.hulls[data.hulls[[filter.feat[[i]][j]]] > filter.threshold[[i]][j], ]
          }
        }
      }
      p <-
        p +
        geom_polygon(data = data.hulls,
                     aes_string(x = "SPATIAL_1.y", y= "SPATIAL_2.y", group = label.id, fill = feat[i]), color = NA) +
        scale_fill_gradient2(low = alpha(col.low[i], alpha = alpha.low[i]),
                             mid = alpha(col.mid[i], alpha = alpha.mid[i]),
                             high = alpha(col.high[i], alpha = alpha.high[i]),
                             midpoint = ifelse(is.null(scale.min.max), max(df_post[feat[i]])/2, mean(scale.min.max[[i]])),
                             na.value = "transparent",
                             limits = c(ifelse(is.null(scale.min.max), min(df_post[feat[i]]), scale.min.max[[i]][1]),
                                        ifelse(is.null(scale.min.max), max(df_post[feat[i]]), scale.min.max[[i]][2]))) +
        ggnewscale::new_scale_fill()
      plotted <- c(plotted, data.points[,length(data.points)])
    }
  } else if("points" %in% plot.type & "hulls" %in% plot.type) {
    for(i in 1:length(feat)) {
      data.points <- df_post[df_post[[feat[i]]] > min.visible[i], ]
      data.hulls <-  df[df[[feat[i]]] > min.visible[i], ]
      if(!is.null(filter.feat)) {
        for(j in 1:length(filter.feat[[i]])) {
          if(filter.feat[[i]][j] == "") {
            next
          } else {
            data.points <- data.points[data.points[[filter.feat[[i]][j]]] > filter.threshold[[i]][j], ]
            data.hulls <- data.hulls[data.hulls[[filter.feat[[i]][j]]] > filter.threshold[[i]][j], ]
          }
        }
      }
      p <-
        p +
        geom_polygon(data = data.hulls,
                     aes_string(x = "SPATIAL_1.y", y= "SPATIAL_2.y", group = label.id, fill = feat[i]), color = NA, show.legend = F) +
        scale_fill_gradient2(low = alpha(col.low[i], alpha = alpha.low[i]),
                             mid = alpha(col.mid[i], alpha = alpha.mid[i]),
                             high = alpha(col.high[i], alpha = alpha.high[i]),
                             midpoint = ifelse(is.null(scale.min.max), max(df_post[feat[i]])/2, mean(scale.min.max[[i]])),
                             na.value = "transparent",
                             limits = c(ifelse(is.null(scale.min.max), min(df_post[feat[i]]), scale.min.max[[i]][1]),
                                        ifelse(is.null(scale.min.max), max(df_post[feat[i]]), scale.min.max[[i]][2]))) +
        ggnewscale::new_scale_fill() +
        geom_point(data = data.points,
                   aes_string(x = "SPATIAL_1", y = "SPATIAL_2", col = feat[i]), shape = shape, fill = NA, size = pt.size*i, stroke = 2*pt.size/3) +
        scale_color_gradient2(low = alpha(col.low[i], alpha = alpha.low[i]),
                              mid = alpha(col.mid[i], alpha = alpha.mid[i]),
                              high = alpha(col.high[i], alpha = alpha.high[i]),
                              midpoint = ifelse(is.null(scale.min.max), max(df_post[feat[i]])/2, mean(scale.min.max[[i]])),
                              na.value = "transparent",
                              limits = c(ifelse(is.null(scale.min.max), min(df_post[feat[i]]), scale.min.max[[i]][1]),
                                         ifelse(is.null(scale.min.max), max(df_post[feat[i]]), scale.min.max[[i]][2]))) +
        ggnewscale::new_scale_color()
      plotted <- c(plotted, data.points[,length(data.points)])
    }
  }

  data.points <- df_post[df_post[, length(df_post)] %in% unique(plotted),]

  # plot labels
  if(show.labels) {
    p <-
      p +
      ggrepel::geom_text_repel(data = data.points, aes_string(x = "SPATIAL_1", y= "SPATIAL_2", label = label.id), color = "black", min.segment.length = 0, max.overlaps = Inf)
  }

  # outline hulls
  if(!is.null(outline.hulls)) {
    p <-
      p +
      geom_polygon(data = dplyr::filter(df, !!sym(label.id) %in% outline.hulls),
                   aes_string(x = "SPATIAL_1.y", y= "SPATIAL_2.y", group = label.id), fill = NA, color = "black")
  }

  # plot scalebar
  if(scalebar) {
    # Calculate scalebar position (bottom-left corner)
    x_range <- range(df_post$SPATIAL_1)
    y_range <- range(df_post$SPATIAL_2)
    x_pos <- x_range[1]# + 0.02 * diff(x_range)  # 2% from left edge
    y_pos <- y_range[2]# - 0.03 * diff(y_range)  # 3% from visual bottom (after flip)

    # Add scalebar to plot
    p <-
      p +
      annotate("rect",
               xmin = x_pos, xmax = x_pos + scalebar / 2,
               ymin = y_pos, ymax = y_pos - scalebar.width,
               fill = "black",
               color = NA
      ) +
      annotate("text",
               x = x_pos + scalebar / 4,
               y = y_pos - scalebar.width * 3, # Adjusted for better visibility
               label = paste(scalebar, "μm"),
               size = scalebar/100,
               color = "black"
      )
  }

  # wrap plot
  p <-
    p +
    coord_fixed(ratio = 1) +
    theme_void() +
    scale_x_continuous(expand = c(0, 5)) +
    scale_y_reverse(expand = c(0, 5)) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(fill = NA),
          plot.margin = margin(5, 5, 5, 5))

  # display plot or return object
  if(plot) {
    print(p)
  } else {
    list(plot = p, cells = data.points, threshold = min.visible)
  }
}

#' Plot cell-cell distances distribution from a ROI plot (previously quantified with get_dist)
#'
#' @param df Data frame of distances (output from get_dist)
#' @param binwidth Binwidth of the histograms
#' @export
plot_dist <- function(df, binwidth = 5) {
  print(
    ggplot(df, aes(x = distance, fill = neighbor_marker)) +
      geom_histogram(binwidth = binwidth, color = NA, alpha = 0.3, position = "identity") +
      stat_bin(binwidth = binwidth, aes(y = after_stat(count), group = neighbor_marker, color = neighbor_marker), geom = "smooth", se = FALSE, linewidth = 0.5, position = "identity") +
      geom_vline(data = dplyr::group_by(df, neighbor_marker, origin_marker) %>% dplyr::summarise(meandist = mean(distance)), aes(xintercept = meandist, col = neighbor_marker), linetype = 2) +
      facet_grid(origin_marker ~ ., scales = "free_y") +
      xlab("distance (microns)") +
      theme_light() +
      theme(strip.background = element_blank(),
            strip.text = element_text(color = "black"))
  )
}
