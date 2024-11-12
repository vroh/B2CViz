#' ROI Selector
#'
#' @param path Path to the image file used for bin2cell segmentation
#' @param post Post-aggregated Seurat object
#' @return list of ROI coordinates
#' @export
roi_selector <- function(path, post) {
  roi_coords_list <- reactiveVal(list())
  redraw_trigger <- reactiveVal(0)

  # Load the original image
  original_img <- load.image(path)

  # Create a lower resolution version for the main display
  max_display_size <- 800  # Maximum width or height for display
  scale_factor <- min(max_display_size / width(original_img),
                      max_display_size / height(original_img),
                      1)  # Ensure we don't upscale small images

  # Use imresize function from imager package
  display_img <- imresize(original_img, scale = scale_factor)

  # Calculate exact scaling factors
  x_scale <- width(original_img) / width(display_img)
  y_scale <- height(original_img) / height(display_img)

  ui <- fluidPage(
    titlePanel("ROI Selector"),
    sidebarLayout(
      sidebarPanel(
        actionButton("add", "Add ROI"),
        actionButton("finish", "Finish and Close"),
        verbatimTextOutput("roi_list"),
        width = 3
      ),
      mainPanel(
        fluidRow(
          column(6, plotOutput("main_plot", height = "800px",
                               click = "plot_click",
                               brush = brushOpts(id = "plot_brush"))),
          column(6, plotOutput("roi_plot", height = "400px"))
        )
      )
    )
  )

  server <- function(input, output, session) {
    output$main_plot <- renderPlot({
      redraw_trigger()  # Depend on this to trigger redraws
      plot(display_img, axes = FALSE)
      rect(min(post@reductions$spatial[[1:dim(post@reductions$spatial)[1]]][,1])*scale_factor,
           min(post@reductions$spatial[[1:dim(post@reductions$spatial)[1]]][,2])*scale_factor,
           max(post@reductions$spatial[[1:dim(post@reductions$spatial)[1]]][,1])*scale_factor,
           max(post@reductions$spatial[[1:dim(post@reductions$spatial)[1]]][,2])*scale_factor,
           border = "green", lwd = 2)

      # Draw existing ROIs in blue
      for (roi in roi_coords_list()) {
        rect(roi$xmin / x_scale, roi$ymin / y_scale,
             roi$xmax / x_scale, roi$ymax / y_scale,
             border = "blue", lwd = 2)
      }

      # Draw current selection in red
      if (!is.null(input$plot_brush)) {
        rect(input$plot_brush$xmin, input$plot_brush$ymin,
             input$plot_brush$xmax, input$plot_brush$ymax,
             border = "red", lwd = 2)
      }
    })

    output$roi_plot <- renderPlot({
      if (!is.null(input$plot_brush)) {
        # Calculate the ROI coordinates for the original image
        orig_xmin <- max(1, round(input$plot_brush$xmin * x_scale))
        orig_xmax <- min(width(original_img), round(input$plot_brush$xmax * x_scale))
        orig_ymin <- max(1, round(input$plot_brush$ymin * y_scale))
        orig_ymax <- min(height(original_img), round(input$plot_brush$ymax * y_scale))

        # Ensure minimum size of 1x1 pixel
        if (orig_xmax <= orig_xmin) orig_xmax <- orig_xmin + 1
        if (orig_ymax <= orig_ymin) orig_ymax <- orig_ymin + 1

        roi <- imsub(original_img,
                     x %inr% c(orig_xmin, orig_xmax),
                     y %inr% c(orig_ymin, orig_ymax))
        plot(roi, axes = FALSE)
      }
    })

    output$roi_list <- renderPrint({
      cat("Selected ROIs:\n")
      for (i in seq_along(roi_coords_list())) {
        cat(sprintf("ROI %d: xmin=%d, xmax=%d, ymin=%d, ymax=%d\n",
                    i,
                    roi_coords_list()[[i]]$xmin,
                    roi_coords_list()[[i]]$xmax,
                    roi_coords_list()[[i]]$ymin,
                    roi_coords_list()[[i]]$ymax))
      }
    })

    observeEvent(input$add, {
      if (!is.null(input$plot_brush)) {
        new_roi <- list(
          xmin = max(1, round(input$plot_brush$xmin * x_scale)),
          xmax = min(width(original_img), round(input$plot_brush$xmax * x_scale)),
          ymin = max(1, round(input$plot_brush$ymin * y_scale)),
          ymax = min(height(original_img), round(input$plot_brush$ymax * y_scale))
        )
        # Ensure minimum size of 1x1 pixel
        if (new_roi$xmax <= new_roi$xmin) new_roi$xmax <- new_roi$xmin + 1
        if (new_roi$ymax <= new_roi$ymin) new_roi$ymax <- new_roi$ymin + 1

        roi_coords_list(c(roi_coords_list(), list(new_roi)))
        showNotification("ROI added successfully.", type = "message")

        # Trigger a redraw
        redraw_trigger(redraw_trigger() + 1)
      } else {
        showNotification("Please select a region before adding.", type = "warning")
      }
    })

    observeEvent(input$finish, {
      stopApp(roi_coords_list())
    })
  }

  result <- runApp(list(ui = ui, server = server))

  # Convert the reactive value to a regular list before returning
  if (is.reactive(result)) {
    return(result())
  } else {
    return(result)
  }
}

#' Set ROI for B2C object
#'
#' @param b2c B2C object
#' @export
set_roi <- function(b2c = NULL) {
  b2c$coord <<- roi_selector(path = b2c$path, post = b2c$post)
}
