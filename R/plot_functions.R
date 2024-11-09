#' Plot MSD data
#'
#' @param object Post-aggregated Seurat object
#' @param features Features to plot
#' @param img Image data
#' @param msd Mulit Spatial Data
#' @param coord Coordinates
#' @param he_alpha Alpha value for H&E image
#' @param pt_size Point size
#' @param plot Whether to plot or return ggplot object
#' @return ggplot object or plot
#' @export
plot_msd <- function(object, features = NULL, img, msd, coord, he_alpha = 0.4, pt_size = 1, plot = TRUE) {
  img <- img[img$y >= coord$xmin & img$y <= coord$xmax & img$x >= coord$ymin & img$x <= coord$ymax,]
  msd[[1]] <- msd[[1]][msd[[1]]$y >= coord$xmin & msd[[1]]$y <= coord$xmax & msd[[1]]$x >= coord$ymin & msd[[1]]$x <= coord$ymax,]


  to_plot <- ggplot(msd[[1]], aes(x = y, y = x)) +
    geom_raster(data = img, aes(fill = color), alpha = he_alpha) +
    geom_point( # create fake legend
      data =
        data.frame(
          features = features,
          x = 0,
          y = 0),
      aes(col = features)
    )

  for(i in 1:length(features)) {
    to_plot <- to_plot +
      geom_point(color = msd[[1]][,features[i]], size = pt_size*i, shape = 21, fill = NA, stroke = pt_size)
  }

  if(plot) {
    print(
      to_plot +
        xlim(coord$xmin, coord$xmax) +
        ylim(coord$ymax, coord$ymin) +
        scale_fill_identity() +
        scale_color_manual(breaks = features, values = msd[[2]], name = "") +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        coord_fixed(ratio = 1) +
        theme_void() +
        theme(legend.position = "right")
    )
  } else {

    to_plot +
      xlim(coord$xmin, coord$xmax) +
      ylim(coord$ymax, coord$ymin) +
      scale_fill_identity() +
      scale_color_manual(breaks = features, values = msd[[2]], name = "") +
      guides(color = guide_legend(override.aes = list(alpha = 1))) +
      coord_fixed(ratio = 1) +
      theme_void() +
      theme(legend.position = "right")
  }
}


#' Wrapper for MSD plotting
#'
#' @param object Post-aggregated Seurat object
#' @param img Image data
#' @param features Features to plot
#' @param coord Coordinates
#' @param intensity Color intensity
#' @param he_alpha Alpha value for H&E image
#' @param pt_size Point size
#' @param plot Whether to plot or return ggplot object
#' @return ggplot object or plot
#' @export
wrap_msd <- function(object, img, features, coord, intensity = 10, he_alpha = 0.4, pt_size = 0.5, plot = TRUE) {
  msd <- get_msd(object = object, features = features, reduction = "spatial")
  msd <- color_msd(msd = msd, intensity = intensity)
  plot_msd(object, features, img = img, msd = msd, coord = coord, he_alpha = he_alpha, pt_size = pt_size, plot = plot)
}



#' Plot B2C object
#'
#' @param b2c B2C object
#' @param features Features to plot
#' @param colors Color palette
#' @param intensity Color intensity
#' @param he_alpha Alpha value for H&E image
#' @param pt_size Point size
#' @param plot.type Type of plot
#' @param label.id Label ID of the bin2cell segmentation feature
#' @export
plot_b2c <- function(b2c = NULL, features = NULL, colors = NULL, intensity = 10, he_alpha = 0.4, pt_size = 0.5, plot.type = c("points", "hulls"), label.id = "labels_he_expanded") {
  if(is.null(b2c$coord)) {
    stop("Run set_roi() on b2c object to set the desired coordinates")
  }

  # load features
  message("loading features")
  if(is.null(b2c[["msd"]])) {
    msd <- get_msd(object = b2c$post, features = features, reduction = "spatial")
    b2c$msd <<- msd
  } else if(all(features %in% colnames(b2c$msd))){

  } else {
    feat_to_add <- features[!features %in% colnames(b2c$msd)]
    feat_to_add_msd <- get_msd(object = b2c$post, features = feat_to_add, reduction = "spatial")
    msd <- cbind(b2c$msd, feat_to_add_msd[, (length(feat_to_add_msd)-(length(feat_to_add)-1)):length(feat_to_add_msd)])
    if(length(feat_to_add) == 1) {
      colnames(msd)[length(msd)] <- feat_to_add
    }
    b2c$msd <<- msd
  }
  if(is.null(b2c$bins) & "hulls" %in% plot.type) {
    message("loading hull data")
    bins <- cbind(b2c$pre@meta.data, FetchData(b2c$pre, vars = c("SPATIAL_1", "SPATIAL_2")))
    bins <- bins[bins[, label.id] > 0,]
    b2c$bins <<- bins
  }

  # color msd
  message("adding color")
  if(!is.null(b2c[["msd"]])) {
    msd_col <- color_msd(msd = b2c$msd, colors = colors, intensity = intensity)
  } else {
    msd_col <- color_msd(msd = msd, colors = colors, intensity = intensity)
  }
  b2c$msd_col <<- msd_col

  # plot
  if(length(plot.type) == 1 & plot.type[1] == "points") {
    message("plotting points")
    plot_msd(object = b2c$post, features = features, img = b2c$img, msd = msd_col, coord = b2c$coord, he_alpha = he_alpha, pt_size = pt_size)
  } else if(length(plot.type) == 1 & plot.type[1] == "hulls") {
    message("generatig hulls")
    if(is.null(b2c$bins)) {
      df <- bins[bins$SPATIAL_1 >= b2c$coord$xmin & bins$SPATIAL_1 <= b2c$coord$xmax & bins$SPATIAL_2 >= b2c$coord$ymin & bins$SPATIAL_2 <= b2c$coord$ymax,]
    } else {
      df <- b2c$bins[b2c$bins$SPATIAL_1 >= b2c$coord$xmin & b2c$bins$SPATIAL_1 <= b2c$coord$xmax & b2c$bins$SPATIAL_2 >= b2c$coord$ymin & b2c$bins$SPATIAL_2 <= b2c$coord$ymax,]
    }
    df[[label.id]] <- factor(df[[label.id]], levels = unique(df[[label.id]]))
    hull <- df %>%
      group_by(across(label.id)) %>%
      slice(chull(SPATIAL_1, SPATIAL_2))


    img <- b2c$img[b2c$img$y >= b2c$coord$xmin & b2c$img$y <= b2c$coord$xmax & b2c$img$x >= b2c$coord$ymin & b2c$img$x <= b2c$coord$ymax,]

    # compute nearest neighbour
    message("computing nearest neighbour (NN)")
    # points
    p1 <- wrap_msd(object = b2c$post, features = features, img = b2c$img, coord = b2c$coord, he_alpha = he_alpha, pt_size = pt_size, plot = FALSE)
    p1 <- p1$data[p1$data$y >= b2c$coord$xmin & p1$data$y <= b2c$coord$xmax & p1$data$x >= b2c$coord$ymin & p1$data$x <= b2c$coord$ymax,]
    if(length(p1) == 4) {
      p1 <- p1[p1[,4] != "#00000000",]
    } else {
      p1 <- p1[apply(p1[,4:length(p1)], 1, function(u) !all(grepl("#00000000", u))),]
    }
    # centroid
    centroids <- hull %>% group_by(across(label.id)) %>% summarise(y = mean(SPATIAL_1), x = mean(SPATIAL_2))
    # NN
    # Sample data frames
    df1 <- select(p1, x, y)
    df2 <- centroids
    # Combine x and y into a matrix for both data frames
    points_df1 <- as.matrix(df1)
    points_df2 <- as.matrix(df2[, c("x", "y")])  # Only use x and y from df2
    # Find the nearest neighbors
    # k = 1 means we want the nearest neighbor
    nearest_indices <- get.knnx(points_df2, points_df1, k = 1)$nn.index
    p1$NN <- nearest_indices

    # generate hull colors
    hull_col <- list()
    for(j in 4:(length(p1)-1)) {
      message(paste("feature:", colnames(p1)[j]))
      current_hull_col <- data.frame()
      current <- p1[,c(j, length(p1))]
      current <- current[current[,1] != "#00000000",]
      for(i in 1:nrow(current)) {
        message(paste0("Computing ", colnames(p1)[j], " NN: ", round(100*i/nrow(current), 1), "% Done"))
        color <- current[i,1]
        cell <- centroids[current[i,]$NN,][[label.id]]
        output <- hull[hull[[label.id]] == cell,] %>% mutate(color = color)
        current_hull_col <- rbind(current_hull_col, output)
      }
      hull_col[[j-3]] <- current_hull_col
    }

    message("plotting hulls")
    to_plot <- ggplot(img, aes(y, x)) +
      geom_raster(aes(fill = color), alpha = 0.4) +
      geom_point( # create fake legend
        data =
          data.frame(
            features = features,
            x = 0,
            y = 0),
        aes(col = features))

    for(i in 1:length(features)) {
      to_plot <- to_plot +
        geom_polygon(data = hull_col[[i]], aes_string(x = "SPATIAL_1", y= "SPATIAL_2", group = label.id), fill = hull_col[[i]]$color, alpha = 1/length(features), color = NA)
    }

    if(is.null(b2c[["msd_col"]])) {
      print(
        to_plot +
          xlim(b2c$coord$xmin, b2c$coord$xmax) +
          ylim(b2c$coord$ymax, b2c$coord$ymin) +
          scale_fill_identity() +
          scale_color_manual(breaks = features, values = msd_col[[2]], name = "") +
          guides(color = guide_legend(override.aes = list(alpha = 1))) +
          coord_fixed(ratio = 1) +
          theme_void() +
          theme(legend.position = "right")
      )
    } else {
      print(
        to_plot +
          xlim(b2c$coord$xmin, b2c$coord$xmax) +
          ylim(b2c$coord$ymax, b2c$coord$ymin) +
          scale_fill_identity() +
          scale_color_manual(breaks = features, values = b2c$msd_col[[2]], name = "") +
          guides(color = guide_legend(override.aes = list(alpha = 1))) +
          coord_fixed(ratio = 1) +
          theme_void() +
          theme(legend.position = "right")
      )
    }

  } else if(all(c("points", "hulls") %in% plot.type)) {
    message("generatig hulls")
    if(is.null(b2c$bins)) {
      df <- bins[bins$SPATIAL_1 >= b2c$coord$xmin & bins$SPATIAL_1 <= b2c$coord$xmax & bins$SPATIAL_2 >= b2c$coord$ymin & bins$SPATIAL_2 <= b2c$coord$ymax,]
    } else {
      df <- b2c$bins[b2c$bins$SPATIAL_1 >= b2c$coord$xmin & b2c$bins$SPATIAL_1 <= b2c$coord$xmax & b2c$bins$SPATIAL_2 >= b2c$coord$ymin & b2c$bins$SPATIAL_2 <= b2c$coord$ymax,]
    }
    df[[label.id]] <- factor(df[[label.id]], levels = unique(df[[label.id]]))
    hull <- df %>%
      group_by(across(label.id)) %>%
      slice(chull(SPATIAL_1, SPATIAL_2))


    img <- b2c$img[b2c$img$y >= b2c$coord$xmin & b2c$img$y <= b2c$coord$xmax & b2c$img$x >= b2c$coord$ymin & b2c$img$x <= b2c$coord$ymax,]

    # compute nearest neighbour
    message("computing nearest neighbour (NN)")
    # points
    p1 <- wrap_msd(object = b2c$post, features = features, img = b2c$img, coord = b2c$coord, he_alpha = he_alpha, pt_size = pt_size, plot = FALSE)
    p1 <- p1$data[p1$data$y >= b2c$coord$xmin & p1$data$y <= b2c$coord$xmax & p1$data$x >= b2c$coord$ymin & p1$data$x <= b2c$coord$ymax,]
    points_to_plot <- p1
    if(length(p1) == 4) {
      p1 <- p1[p1[,4] != "#00000000",]
    } else {
      p1 <- p1[apply(p1[,4:length(p1)], 1, function(u) !all(grepl("#00000000", u))),]
    }
    # centroid
    centroids <- hull %>% group_by(across(label.id)) %>% summarise(y = mean(SPATIAL_1), x = mean(SPATIAL_2))
    # NN
    # Sample data frames
    df1 <- select(p1, x, y)
    df2 <- centroids
    # Combine x and y into a matrix for both data frames
    points_df1 <- as.matrix(df1)
    points_df2 <- as.matrix(df2[, c("x", "y")])  # Only use x and y from df2
    # Find the nearest neighbors
    # k = 1 means we want the nearest neighbor
    nearest_indices <- get.knnx(points_df2, points_df1, k = 1)$nn.index
    p1$NN <- nearest_indices

    # generate hull colors
    hull_col <- list()
    for(j in 4:(length(p1)-1)) {
      message(paste("feature:", colnames(p1)[j]))
      current_hull_col <- data.frame()
      current <- p1[,c(j, length(p1))]
      current <- current[current[,1] != "#00000000",]
      for(i in 1:nrow(current)) {
        message(paste0("Computing ", colnames(p1)[j], " NN: ", round(100*i/nrow(current), 1), "% Done"))
        color <- current[i,1]
        cell <- centroids[current[i,]$NN,][[label.id]]
        output <- hull[hull[[label.id]] == cell,] %>% mutate(color = color)
        current_hull_col <- rbind(current_hull_col, output)
      }
      hull_col[[j-3]] <- current_hull_col
    }

    message("plotting points and hulls")
    to_plot <- ggplot(img, aes(y, x)) +
      geom_raster(aes(fill = color), alpha = 0.4) +
      geom_point( # create fake legend
        data =
          data.frame(
            features = features,
            x = 0,
            y = 0),
        aes(col = features))

    for(i in 1:length(features)) {
      to_plot <- to_plot +
        geom_polygon(data = hull_col[[i]], aes_string(x = "SPATIAL_1", y= "SPATIAL_2", group = label.id), fill = hull_col[[i]]$color, alpha = 1/length(features), color = NA) +
        geom_point(data = points_to_plot, color = points_to_plot[,features[i]], size = pt_size*i, shape = 21, fill = NA, stroke = pt_size)
    }

    if(is.null(b2c[["msd_col"]])) {
      print(
        to_plot +
          xlim(b2c$coord$xmin, b2c$coord$xmax) +
          ylim(b2c$coord$ymax, b2c$coord$ymin) +
          scale_fill_identity() +
          scale_color_manual(breaks = features, values = msd_col[[2]], name = "") +
          guides(color = guide_legend(override.aes = list(alpha = 1))) +
          coord_fixed(ratio = 1) +
          theme_void() +
          theme(legend.position = "right")
      )
    } else {
      print(
        to_plot +
          xlim(b2c$coord$xmin, b2c$coord$xmax) +
          ylim(b2c$coord$ymax, b2c$coord$ymin) +
          scale_fill_identity() +
          scale_color_manual(breaks = features, values = b2c$msd_col[[2]], name = "") +
          guides(color = guide_legend(override.aes = list(alpha = 1))) +
          coord_fixed(ratio = 1) +
          theme_void() +
          theme(legend.position = "right")
      )
    }
  }

}
