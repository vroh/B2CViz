#' Get MultiSpatialData (MSD) data
#'
#' @param ... Parameters for FeaturePlot
#' @return MSD data frame
#' @export
get_msd <- function(...) {
  pw <- tryCatch(FeaturePlot(...), warning=function(w) w)
  if(is(pw,"warning")) {
    stop(pw)
  }
  df <- list()
  for(i in 1:length(pw)) {
    df[[i]] <- pw[[i]]$data
  }
  df <- do.call("cbind", df)
  df <- df[,!duplicated(colnames(df))]
  colnames(df)[1:2] <- c("y", "x")
  df
}

#' Color MultiSpatialData (MSD) data
#'
#' @param msd MSD data frame
#' @param colors Color palette
#' @param intensity Color intensity
#' @return List of colored MSD data and legend colors
#' @export
color_msd <- function(msd, colors = NULL, intensity = 10) {
  if(is.null(colors)) {
    colors <-
      data.frame(
        mid = c( "orange", "indianred4",  "dodgerblue4", "seagreen"),#, "magenta4"), # a bit too close to red
        high = c("gold",   "indianred1", "dodgerblue",  "seagreen1")#,  "magenta2")
      )
  }
  legend_colors <- NULL
  for(i in 4:length(msd)) {
    mid <- colors$mid[ifelse(i%%nrow(colors) == 0, nrow(colors), i%%nrow(colors))]
    high <- colors$high[ifelse(i%%nrow(colors) == 0, nrow(colors), i%%nrow(colors))]
    msd[,i] <- alpha(circlize::colorRamp2(c(min(msd[,i]), mean(msd[,i]), max(msd[,i])),
                                          c("#00000000",
                                            mid,
                                            high))
                     (msd[,i]), intensity*msd[,i]/max(msd[,i]))
    legend_colors <- c(legend_colors, mid)
  }
  list(msd, legend_colors)
}
