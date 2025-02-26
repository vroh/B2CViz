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
#' @export
overview_b2c <- function(b2c, feat, pt.size = 0.001, he_alpha = 0.4, col.low = "gray", col.high = "seagreen2") {

  if(is.null(b2c$img_sd)) {
    stop("run scaledown_img() on b2c object first")
  }

  df <- FeaturePlot(b2c$post, feat, reduction = "spatial")[[1]]$data
  colnames(df)[4] <- "feat"
  print(
    ggplot(filter(df, feat > 0), aes(y, x)) +
      geom_raster(data = b2c$img_sd, aes(fill = color), alpha = he_alpha) +
      geom_point(aes(x = SPATIAL_1, y = SPATIAL_2, col = feat), alpha = scale(filter(df, feat > 0)$feat), size = pt.size) +
      coord_fixed(ratio = 1) +
      scale_fill_identity() +
      scale_color_continuous(name = feat, low = col.low, col.high = "seagreen2") +
      xlim(min(df$SPATIAL_1), max(df$SPATIAL_1)) +
      ylim(max(df$SPATIAL_2), min(df$SPATIAL_2)) +
      theme_void() +
      ggtitle(feat) +
      theme(plot.title = element_text(hjust = 0.5)))
}

#' Plot B2C object
#'
#' @param b2c B2C object
#' @param feat Features to plot
#' @param label.id Label ID of the bin2cell segmentation feature
#' @param min.visible Feature value threshold for visualization (single value or vector for each feature)
#' @param col.mid Mid color for colorscale (single value or vector for each feature)
#' @param col.high High color for colorscale (single value or vector for each feature)
#' @param alpha.mid Mid alpha value for colorscale (single value or vector for each feature)
#' @param alpha.high High alpha value for colorscale (single value or vector for each feature)
#' @param pt_size Point size
#' @param he_alpha Alpha value for H&E image
#' @param title Plot title
#' @param plot.type Type of plot (points or hulls or both)
#' @param outline.hulls Character vector of hulls to outline (by expanded labels ID)
#' @param show.labels Whether or not to plot the hulls labels
#' @param plot Whether to display or return the plot
#' @export
plot_b2c <- function(b2c, feat, label.id = "labels_he_expanded", min.visible = 0,
                     col.mid = "orangered", col.high = "orangered", alpha.mid = 0, alpha.high = 1,
                     pt.size = 1, he_alpha = 0.3, title = NULL,
                     plot.type = c("points", "hulls"), outline.hulls = NULL, show.labels = F, plot = T) {

  # adjust parameters length
  if(length(min.visible) == 1) {
    min.visible <- rep(min.visible, length(feat))
  }
  if(length(alpha.mid) == 1) {
    alpha.mid <- rep(alpha.mid, length(feat))
  }
  if(length(alpha.high) == 1) {
    alpha.high <- rep(alpha.high, length(feat))
  }

  # fetch data
  df_post <- FetchData(b2c$post, vars = c("SPATIAL_1", "SPATIAL_2", feat))
  df_post[label.id] <- row.names(df_post)
  if("hulls" %in% plot.type) {
    df_pre <- FetchData(b2c$pre, vars = c("SPATIAL_1", "SPATIAL_2", label.id)) %>%
      group_by(across(label.id)) %>%
      slice(chull(SPATIAL_1, SPATIAL_2))
    df <- merge(df_post, df_pre, by = label.id)
  }

  # plot H&E
  p <-
    ggplot(df_post) +
    geom_raster(data = b2c$img, aes(x = y, y = x, fill = color), alpha = he_alpha) +
    scale_fill_identity() +
    ggnewscale::new_scale_fill()

  # plot cells
  if("points" %in% plot.type & !("hulls" %in% plot.type)) {
    for(i in 1:length(feat)) {
      p <-
        p +
        geom_point(data = dplyr::filter(df_post, !!sym(feat[i]) > min.visible[i]),
                   aes_string(x = "SPATIAL_1", y = "SPATIAL_2", col = feat[i]), shape = 21, fill = NA, size = pt.size*i, stroke = 2*pt.size/3) +
        scale_color_gradient2(mid = alpha(col.mid[i], alpha = alpha.mid[i]), high = alpha(col.high[i], alpha = alpha.high[i])) +
        ggnewscale::new_scale_color()
    }
  } else if("hulls" %in% plot.type & !("points" %in% plot.type)) {
    for(i in 1:length(feat)) {
      p <-
        p +
        geom_polygon(data = dplyr::filter(df, !!sym(feat[i]) > min.visible[i]),
                     aes_string(x = "SPATIAL_1.y", y= "SPATIAL_2.y", group = label.id, fill = feat[i]), color = NA) +
        scale_fill_gradient2(mid = alpha(col.mid[i], alpha = alpha.mid[i]), high = alpha(col.high[i], alpha = alpha.high[i])) +
        ggnewscale::new_scale_fill()
    }
  } else if("points" %in% plot.type & "hulls" %in% plot.type) {
    for(i in 1:length(feat)) {
      p <-
        p +
        geom_polygon(data = dplyr::filter(df, !!sym(feat[i]) > min.visible[i]),
                     aes_string(x = "SPATIAL_1.y", y= "SPATIAL_2.y", group = label.id, fill = feat[i]), color = NA, show.legend = F) +
        scale_fill_gradient2(mid = alpha(col.mid[i], alpha = alpha.mid[i]), high = alpha(col.high[i], alpha = alpha.high[i])) +
        ggnewscale::new_scale_fill() +
        geom_point(data = dplyr::filter(df_post, !!sym(feat[i]) > min.visible[i]),
                   aes_string(x = "SPATIAL_1", y = "SPATIAL_2", col = feat[i]), shape = 21, fill = NA, size = pt.size*i, stroke = 2*pt.size/3) +
        scale_color_gradient2(mid = alpha(col.mid[i], alpha = alpha.mid[i]), high = alpha(col.high[i], alpha = alpha.high[i])) +
        ggnewscale::new_scale_color()
    }
  }

  # plot labels
  if(show.labels) {
    for(i in 1:length(feat)) {
      p <-
        p +
        ggrepel::geom_text_repel(data = dplyr::filter(df_post, !!sym(feat[i]) > min.visible[i]), aes_string(x = "SPATIAL_1", y= "SPATIAL_2", label = label.id), color = "black", min.segment.length = 0, max.overlaps = Inf)
    }
  }

  # outline hulls
  if(!is.null(outline.hulls)) {
    p <-
      p +
      geom_polygon(data = dplyr::filter(df, !!sym(label.id) %in% outline.hulls),
                   aes_string(x = "SPATIAL_1.y", y= "SPATIAL_2.y", group = label.id), fill = NA, color = "black")
  }

  # wrap plot
  p <-
    p +
    coord_fixed(ratio = 1) +
    theme_void() +
    scale_y_reverse() +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5))

  # display plot or return object
  if(plot) {
    print(p)
  } else {
    p
  }
}

