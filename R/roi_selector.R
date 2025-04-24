#' @import(shiny)
#' @import(imager)
#' @import(dplyr)
NULL

#' ROI Selector
#'
#' @param path Path to the image file used for bin2cell segmentation
#' @param post Post-aggregated Seurat object
#' @return list of ROI coordinates
#' @export
roi_selector <- function(path, post = NULL) {
  # --- State Management ---
  roi_coords_list <- reactiveVal(list()) # Stores final ROIs
  redraw_trigger <- reactiveVal(0) # Triggers plot redraws (mainly for polygon updates)

  # Polygon drawing state
  drawing_mode <- reactiveVal(FALSE) # Are we currently drawing a polygon?
  current_polygon_points <- reactiveVal(list()) # Points added to the current polygon (original coords)
  selected_rectangle_orig <- reactiveVal(NULL) # Stores the *original* image coords of the brushed rectangle
  # selected_rectangle_disp <- reactiveVal(NULL) # Store display coords if needed elsewhere

  # --- Add caching for the ROI image snippet ---
  cached_roi_img <- reactiveVal(NULL)

  # --- Image Loading and Scaling ---
  if (!file.exists(path)) {
    stop("Image file not found at path: ", path)
  }
  original_img <- tryCatch(
    load.image(path),
    error = function(e) {
      stop("Failed to load image. Error: ", e$message)
    }
  )

  max_display_size <- 800
  scale_factor <- min(max_display_size / width(original_img),
                      max_display_size / height(original_img),
                      1)

  display_img <- imresize(original_img, scale = scale_factor)
  x_scale <- width(original_img) / width(display_img)
  y_scale <- height(original_img) / height(display_img)

  # --- UI Definition (Identical to previous version) ---
  ui <- fluidPage(
    titlePanel("ROI Selector"),
    sidebarLayout(
      sidebarPanel(
        h4("Instructions:"),
        p("1. Drag a rectangle on the main image."),
        p("2. Click 'Add Rectangle ROI' OR"),
        p("3. Click 'Start/Reset Polygon' to draw within the zoom view."),
        p("4. Click points in the zoom view. Click near start to finish."),
        p("5. Click 'Add Polygon ROI'."),
        hr(),
        actionButton("add_rect", "Add Rectangle ROI", icon = icon("vector-square")),
        actionButton("start_poly", "Start/Reset Polygon", icon = icon("draw-polygon")),
        actionButton("add_poly", "Add Polygon ROI", icon = icon("check")),
        actionButton("clear_poly", "Clear Polygon Drawing", icon = icon("eraser"), disabled = TRUE),
        hr(),
        actionButton("finish", "Finish and Close", icon = icon("sign-out-alt")),
        hr(),
        h4("Selected ROIs:"),
        verbatimTextOutput("roi_list"),
        width = 3
      ),
      mainPanel(
        fluidRow(
          column(6,
                 h4("Full Image"),
                 plotOutput("main_plot", height = "auto",
                            click = "plot_click",
                            brush = brushOpts(id = "plot_brush", resetOnNew = FALSE))
          ),
          column(6,
                 h4("Zoomed Region / Polygon Drawing"),
                 # Make height slightly larger maybe? or keep auto
                 plotOutput("roi_plot", height = "auto",
                            click = "roi_plot_click")
          )
        )
      )
    )
  )

  # --- Server Logic ---
  server <- function(input, output, session) {

    # Get spatial extent if 'post' object is provided
    spatial_extent_orig <- reactive({
      if (!is.null(post) && "spatial" %in% names(post@reductions) && dim(post@reductions$spatial)[1] > 0) {
        coords <- post@reductions$spatial[[1:dim(post@reductions$spatial)[1]]]
        list(
          xmin = min(coords[,1]), ymin = min(coords[,2]),
          xmax = max(coords[,1]), ymax = max(coords[,2])
        )
      } else { NULL }
    })

    # -- Main Plot Rendering --
    output$main_plot <- renderPlot({
      # Use redraw_trigger() OR selected_rectangle_orig() change to redraw
      # Dependency on roi_coords_list() is implicit via redraw_trigger increment on add
      redraw_trigger()

      par(mar = c(0, 0, 0, 0))
      plot(display_img, axes = FALSE, interp = FALSE)

      ext <- spatial_extent_orig()
      if (!is.null(ext)) {
        rect(ext$xmin / x_scale, ext$ymin / y_scale,
             ext$xmax / x_scale, ext$ymax / y_scale,
             border = "green", lwd = 2, lty = "dashed")
      }

      current_rois <- roi_coords_list()
      if (length(current_rois) > 0) {
        for (roi in current_rois) {
          if (roi$type == "rectangle") {
            coords <- roi$coords
            rect(coords$xmin / x_scale, coords$ymin / y_scale,
                 coords$xmax / x_scale, coords$ymax / y_scale,
                 border = "green", lwd = 2)
          } else if (roi$type == "polygon") {
            poly_points_scaled <- lapply(roi$points, function(p) list(x = p$x / x_scale, y = p$y / y_scale))
            polygon(sapply(poly_points_scaled, `[[`, "x"), sapply(poly_points_scaled, `[[`, "y"),
                    border = "green", lwd = 2)
          }
        }
      }

      # Use selected_rectangle_disp() if available for consistency
      current_brush <- input$plot_brush # Read reactive once
      if (!is.null(current_brush)) {
        rect(current_brush$xmin, current_brush$ymin,
             current_brush$xmax, current_brush$ymax,
             border = "red", lwd = 2)
      }
    }, height = function() {
      max(400, min(max_display_size, round(height(display_img) * (session$clientData$output_main_plot_width / width(display_img)))))
    })


    # -- Observe Brush Selection -> Update Cache --
    observeEvent(input$plot_brush, {
      brush <- input$plot_brush # Read reactive once
      if (!is.null(brush)) {
        # Calculate original coordinates (ensure validity)
        orig_xmin <- max(1, round(brush$xmin * x_scale))
        orig_xmax <- min(width(original_img), round(brush$xmax * x_scale))
        orig_ymin <- max(1, round(brush$ymin * y_scale))
        orig_ymax <- min(height(original_img), round(brush$ymax * y_scale))

        if (orig_xmax <= orig_xmin) orig_xmax <- orig_xmin + 1
        if (orig_ymax <= orig_ymin) orig_ymax <- orig_ymin + 1
        orig_xmax <- min(width(original_img), orig_xmax)
        orig_ymax <- min(height(original_img), orig_ymax)

        sel_orig <- list(xmin = orig_xmin, xmax = orig_xmax, ymin = orig_ymin, ymax = orig_ymax)
        selected_rectangle_orig(sel_orig) # Store original coords

        # --- Cache the image subset ---
        # This is the key change for caching
        if (sel_orig$xmax > sel_orig$xmin && sel_orig$ymax > sel_orig$ymin) {
          roi_subset <- tryCatch(
            imsub(original_img,
                  x %inr% c(sel_orig$xmin, sel_orig$xmax),
                  y %inr% c(sel_orig$ymin, sel_orig$ymax)),
            error = function(e) {
              showNotification(paste("Error subsetting image:", e$message), type="error")
              return(NULL)
            }
          )
          cached_roi_img(roi_subset) # Update the cache
        } else {
          cached_roi_img(NULL) # Invalid selection dimensions
          showNotification("Selected region is too small.", type = "warning")
        }
        # --- End caching update ---

        # Reset polygon state for the new selection
        drawing_mode(FALSE)
        current_polygon_points(list())
        updateActionButton(session, "clear_poly", disabled = TRUE)
        updateActionButton(session, "add_poly", disabled = TRUE)
        updateActionButton(session, "add_rect", disabled = FALSE)
        updateActionButton(session, "start_poly", disabled = FALSE)

      } else {
        # Brush cleared
        selected_rectangle_orig(NULL)
        cached_roi_img(NULL) # Clear cache
        drawing_mode(FALSE)
        current_polygon_points(list())
        updateActionButton(session, "clear_poly", disabled = TRUE)
        updateActionButton(session, "add_poly", disabled = TRUE)
        updateActionButton(session, "add_rect", disabled = TRUE)
        updateActionButton(session, "start_poly", disabled = TRUE)
      }
      # Trigger roi_plot redraw explicitly when brush changes (uses cache)
      # No need for redraw_trigger here, plot depends on cached_roi_img & selected_rectangle_orig
      # redraw_trigger(redraw_trigger() + 1) # Not strictly needed for cache update
    }, ignoreNULL = FALSE) # Trigger even when brush becomes NULL


    # -- ROI Plot Rendering (Uses Cache) --
    output$roi_plot <- renderPlot({
      # Depend on redraw_trigger for polygon updates, and cached_roi_img for image data
      redraw_trigger()
      roi <- cached_roi_img() # Get from cache
      rect_orig <- selected_rectangle_orig() # Get original coords for offsets

      if (!is.null(roi) && !is.null(rect_orig) && width(roi) > 0 && height(roi) > 0) {
        par(mar = c(0, 0, 0, 0))
        # --- Plot the cached image ---
        plot(roi, axes = FALSE, interp = FALSE)

        # --- Overlay polygon drawing (same logic as before) ---
        poly_points <- current_polygon_points()
        if (length(poly_points) > 0) {
          # Convert original polygon coords to ROI plot coords (relative to subset)
          roi_plot_points <- lapply(poly_points, function(p) {
            list(x = p$x - rect_orig$xmin + 1,
                 y = p$y - rect_orig$ymin + 1)
          })

          # Draw points
          points(sapply(roi_plot_points, `[[`, "x"), sapply(roi_plot_points, `[[`, "y"),
                 col = "red", pch = 19, cex = 1.5)

          # Draw lines connecting points
          if (length(roi_plot_points) > 1) {
            lines(sapply(roi_plot_points, `[[`, "x"), sapply(roi_plot_points, `[[`, "y"),
                  col = "red", lwd = 2)
          }

          # Indicate closure possibility / status
          if (drawing_mode() && length(roi_plot_points) >= 3) {
            first_pt <- roi_plot_points[[1]]
            points(first_pt$x, first_pt$y, col = "orange", pch = 1, cex = 3, lwd=2)
            title(main = "Click points. Click near orange circle to close.", line = -1.1, cex.main = 1.5, col.main = "black")
          } else if (drawing_mode()) {
            title(main = "Click points to draw polygon.", line = -1.1, cex.main = 1.5, col.main = "black")
          } else if (!drawing_mode() && length(poly_points) > 0) {
            lines(c(sapply(roi_plot_points, `[[`, "x"), roi_plot_points[[1]]$x),
                  c(sapply(roi_plot_points, `[[`, "y"), roi_plot_points[[1]]$y),
                  col = "green", lwd = 2, lty="dashed")
            title(main = "Polygon finalized. Click 'Add Polygon ROI'.", line = -1.1, cex.main = 1.5, col.main = "black")
          }
        } else if (drawing_mode()) {
          title(main = "Click points to start drawing polygon.", line = -1.1, cex.main = 1.5, col.main = "black")
        }
        # --- End polygon overlay ---
      } else {
        # No selection or invalid cache
        plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "Select a region in the main plot")
      }
    }, height = function(){
      roi <- cached_roi_img() # Check cache for dimensions
      if (!is.null(roi)){
        # Calculate height based on cached roi aspect ratio
        max(200, min(400, round(height(roi) * (session$clientData$output_roi_plot_width / width(roi)))))
      } else {
        400 # Default height
      }
    })


    # -- Handle Clicks in ROI Plot for Polygon Drawing --
    observeEvent(input$roi_plot_click, {
      # Important: Use isolate() to read values that should not trigger this observer
      # Only trigger based on the click event itself.
      if (isolate(drawing_mode()) && !is.null(input$roi_plot_click)) {

        click <- input$roi_plot_click
        rect_orig <- isolate(selected_rectangle_orig()) # Read without dependency

        if(is.null(rect_orig)) return() # Should not happen if drawing_mode is TRUE, but safe check

        # Convert click coordinates (plot data coords) to original image coordinates.
        # Assuming plot data coords map directly to original subset indices.
        orig_x <- round(rect_orig$xmin + click$x - 1)
        orig_y <- round(rect_orig$ymin + click$y - 1)

        # Clamp coordinates within the original image bounds (optional but safer)
        orig_x <- max(1, min(width(isolate(original_img)), orig_x))
        orig_y <- max(1, min(height(isolate(original_img)), orig_y))

        new_point <- list(x = orig_x, y = orig_y)

        # Use functional update for safer modification of reactiveVal list
        current_points <- current_polygon_points() # Get current list

        close_threshold_orig <- 10 # Distance threshold in original coordinates

        closed_polygon <- FALSE
        if (length(current_points) >= 3) {
          first_point <- current_points[[1]]
          distance <- sqrt((new_point$x - first_point$x)^2 + (new_point$y - first_point$y)^2)

          if (distance < close_threshold_orig) {
            drawing_mode(FALSE) # Stop drawing mode
            closed_polygon <- TRUE
            showNotification("Polygon closed. Click 'Add Polygon ROI'.", type = "message")
            updateActionButton(session, "add_poly", disabled = FALSE)
            updateActionButton(session, "clear_poly", disabled = FALSE)
          }
        }

        # Add the point ONLY if not closing the polygon
        if (!closed_polygon) {
          current_polygon_points(c(current_points, list(new_point)))
        }

        # --- Trigger redraw for polygon update ---
        # This is essential for seeing points/lines appear
        redraw_trigger(redraw_trigger() + 1)
      }
    })

    # -- Button Actions (Mostly Unchanged, ensure redraw_trigger on Add/Clear) --

    # Add Rectangle ROI
    observeEvent(input$add_rect, {
      rect_orig <- selected_rectangle_orig()
      if (!is.null(rect_orig)) {
        new_roi <- list(type = "rectangle", coords = rect_orig)
        roi_coords_list(c(roi_coords_list(), list(new_roi))) # Update list
        showNotification("Rectangle ROI added.", type = "message")

        # Reset state
        selected_rectangle_orig(NULL)
        cached_roi_img(NULL) # Clear cache
        drawing_mode(FALSE)
        current_polygon_points(list())
        updateActionButton(session, "clear_poly", disabled = TRUE)
        updateActionButton(session, "add_poly", disabled = TRUE)
        updateActionButton(session, "add_rect", disabled = TRUE)
        updateActionButton(session, "start_poly", disabled = TRUE)
        session$resetBrush("plot_brush") # Attempt to clear brush visually

        redraw_trigger(redraw_trigger() + 1) # Redraw main plot
      } else {
        showNotification("Please select a rectangular region first.", type = "warning")
      }
    })

    # Start Polygon Drawing
    observeEvent(input$start_poly, {
      # Check if a region is selected and cached
      if (!is.null(selected_rectangle_orig()) && !is.null(cached_roi_img())) {
        drawing_mode(TRUE)
        current_polygon_points(list()) # Reset points
        updateActionButton(session, "clear_poly", disabled = FALSE)
        updateActionButton(session, "add_poly", disabled = TRUE)
        updateActionButton(session, "add_rect", disabled = TRUE)
        showNotification("Polygon drawing started. Click points in the zoom view.", type = "message")
        redraw_trigger(redraw_trigger() + 1) # Update roi_plot title/state
      } else {
        showNotification("Please select a valid rectangular region first.", type = "warning")
      }
    })

    # Clear Polygon Drawing
    observeEvent(input$clear_poly, {
      drawing_mode(FALSE)
      current_polygon_points(list())
      updateActionButton(session, "clear_poly", disabled = TRUE)
      updateActionButton(session, "add_poly", disabled = TRUE)
      updateActionButton(session, "add_rect", disabled = is.null(selected_rectangle_orig())) # Re-enable if selection exists
      showNotification("Polygon drawing cleared.", type = "message")
      redraw_trigger(redraw_trigger() + 1) # Update roi_plot display
    })

    # Add Polygon ROI
    observeEvent(input$add_poly, {
      # Use isolate to prevent dependency loops if poly_points changes elsewhere
      poly_points <- isolate(current_polygon_points())
      # Check if drawing is finished and polygon is valid
      if (!isolate(drawing_mode()) && length(poly_points) >= 3) {
        new_roi <- list(type = "polygon", points = poly_points)
        roi_coords_list(c(roi_coords_list(), list(new_roi))) # Update list
        showNotification("Polygon ROI added.", type = "message")

        # Reset state
        selected_rectangle_orig(NULL)
        cached_roi_img(NULL) # Clear cache
        drawing_mode(FALSE)
        current_polygon_points(list())
        updateActionButton(session, "clear_poly", disabled = TRUE)
        updateActionButton(session, "add_poly", disabled = TRUE)
        updateActionButton(session, "add_rect", disabled = TRUE)
        updateActionButton(session, "start_poly", disabled = TRUE)
        session$resetBrush("plot_brush") # Attempt to clear brush

        redraw_trigger(redraw_trigger() + 1) # Redraw main plot

      } else if (isolate(drawing_mode())) {
        showNotification("Please finish drawing the polygon (click near start point).", type = "warning")
      } else {
        showNotification("No valid polygon drawn or selected (need >= 3 points).", type = "warning")
      }
    })

    # -- ROI List Output (Unchanged) --
    output$roi_list <- renderPrint({
      rois <- roi_coords_list()
      # ... (same as before) ...
      if (length(rois) == 0) {
        cat("No ROIs selected yet.\n")
      } else {
        cat("Selected ROIs:\n")
        for (i in seq_along(rois)) {
          roi <- rois[[i]]
          if (roi$type == "rectangle") {
            cat(sprintf("ROI %d (Rect): xmin=%d, xmax=%d, ymin=%d, ymax=%d\n",
                        i, roi$coords$xmin, roi$coords$xmax, roi$coords$ymin, roi$coords$ymax))
          } else if (roi$type == "polygon") {
            cat(sprintf("ROI %d (Poly): %d points\n", i, length(roi$points)))
          }
        }
      }
    })

    # -- Finish Button (Unchanged) --
    observeEvent(input$finish, {
      stopApp(roi_coords_list())
    })

    # -- Initial Button States (Unchanged) --
    observe({
      updateActionButton(session, "add_rect", disabled = TRUE)
      updateActionButton(session, "start_poly", disabled = TRUE)
      updateActionButton(session, "add_poly", disabled = TRUE)
      updateActionButton(session, "clear_poly", disabled = TRUE)
    })

  } # End server function

  # --- Run the App (Unchanged) ---
  app_port <- httpuv::randomPort()
  print(paste("Starting Shiny app on port", app_port))
  result <- runApp(list(ui = ui, server = server), port = app_port, launch.browser = getOption("shiny.launch.browser", interactive()))

  if (is.reactive(result)) { return(result()) } else { return(result) }
}

#' Set ROI for B2C object
#'
#' @param b2c B2C object
#' @return B2C object
#' @export
set_roi <- function(b2c = NULL) {
  b2c$coord <- roi_selector(path = b2c$path, post = b2c$post)
  b2c
}

# Define the point-in-polygon function outside the main function for clarity
point_in_polygon <- function(point_x, point_y, polygon_df) {
  polygon_x <- polygon_df$x
  polygon_y <- polygon_df$y
  n <- nrow(polygon_df)
  if (n < 3) return(FALSE) # Need at least 3 vertices for a polygon

  inside <- FALSE
  p1x <- polygon_x[1]
  p1y <- polygon_y[1]
  for (i in 1:n) {
    p2x <- polygon_x[ifelse(i == n, 1, i + 1)] # Wrap around to the first point
    p2y <- polygon_y[ifelse(i == n, 1, i + 1)]

    # Check if the point is on the same horizontal line as a segment
    # This handles horizontal edges correctly in the ray casting algorithm
    if (point_y == p1y && point_y == p2y && point_x >= min(p1x, p2x) && point_x <= max(p1x, p2x)) {
      return(TRUE) # Point is on a horizontal boundary segment
    }
    # Check if the point is on the same vertical line as a segment
    if (point_x == p1x && point_x == p2x && point_y >= min(p1y, p2y) && point_y <= max(p1y, p2y)) {
      return(TRUE) # Point is on a vertical boundary segment
    }

    # Ray casting algorithm part
    if (((p1y <= point_y && point_y < p2y) || (p2y <= point_y && point_y < p1y)) &&
        (point_x < (p2x - p1x) * (point_y - p1y) / (p2y - p1y) + p1x)) {
      inside <- !inside
    }
    p1x <- p2x
    p1y <- p2y
  }
  return(inside)
}

#' Crop B2C object to the dimension of selected ROI
#'
#' @param b2c B2C object
#' @param roi ROI id to crop from
#' @param label.id Label ID of the bin2cell segmentation feature
#' @return B2C object
#' @export
crop_b2c <- function(b2c, roi = 1, label.id = "labels_he_expanded") {
  # --- 1. Fetch initial data ---
  coord_type <- b2c$coord[[roi]]$type
  cells_data <- FetchData(b2c$post, c("SPATIAL_1", "SPATIAL_2"))
  # Ensure we have rownames to link back to the original object
  cells_data$original_rownames <- rownames(cells_data)

  # --- 2. Determine Coordinates and Apply Filtering ---
  if (coord_type == "rectangle") {
    xmin <- b2c$coord[[roi]]$coords$xmin
    xmax <- b2c$coord[[roi]]$coords$xmax
    ymin <- b2c$coord[[roi]]$coords$ymin
    ymax <- b2c$coord[[roi]]$coords$ymax

    # Filter cells based on rectangle
    cells_to_keep_df <- filter(cells_data,
                               SPATIAL_1 >= xmin,
                               SPATIAL_1 <= xmax,
                               SPATIAL_2 >= ymin,
                               SPATIAL_2 <= ymax)

    # Filter image data based on rectangle
    img_to_keep <- filter(b2c$img,
                          y >= xmin,
                          y <= xmax,
                          x >= ymin,
                          x <= ymax)

  } else if (coord_type == "polygon") {
    # Extract polygon points
    points_list <- b2c$coord[[roi]]$points
    points_df <- data.frame(
      x = sapply(points_list, function(p) p$x),
      y = sapply(points_list, function(p) p$y)
    )

    # Calculate bounding box of the polygon
    xmin_poly <- min(points_df$x)
    xmax_poly <- max(points_df$x)
    ymin_poly <- min(points_df$y)
    ymax_poly <- max(points_df$y)

    # --- Step 1: Broad Phase - Filter by Bounding Box ---
    # Pre-filter cells using the polygon's bounding box
    # Note: Using the variable names as in FetchData (SPATIAL_1, SPATIAL_2)
    pre_filtered_cells <- filter(cells_data,
                                 SPATIAL_1 >= xmin_poly,
                                 SPATIAL_1 <= xmax_poly,
                                 SPATIAL_2 >= ymin_poly,
                                 SPATIAL_2 <= ymax_poly)

    # Pre-filter image data using the polygon's bounding box
    # Note: Using the variable names as in b2c$img (y, x)
    # Assuming b2c$img$y corresponds to SPATIAL_1 and b2c$img$x to SPATIAL_2
    pre_filtered_img <- filter(b2c$img,
                               y >= xmin_poly,
                               y <= xmax_poly,
                               x >= ymin_poly,
                               x <= ymax_poly)

    # --- Step 2: Narrow Phase - Apply Point-in-Polygon Test ---
    if (nrow(pre_filtered_cells) > 0) {
      # Apply the precise polygon test only to the pre-filtered cells
      inside_flags_cells <- apply(pre_filtered_cells, 1, function(row) {
        point_in_polygon(point_x = as.numeric(row["SPATIAL_1"]),
                         point_y = as.numeric(row["SPATIAL_2"]),
                         polygon_df = points_df)
      })
      cells_to_keep_df <- pre_filtered_cells[inside_flags_cells, ]
    } else {
      # Handle case where bounding box is empty or outside all points
      cells_to_keep_df <- pre_filtered_cells[FALSE, ] # Return empty dataframe structure
    }

    if (nrow(pre_filtered_img) > 0) {
      # Apply the precise polygon test only to the pre-filtered image points
      inside_flags_img <- apply(pre_filtered_img, 1, function(row) {
        point_in_polygon(point_x = as.numeric(row["y"]), # Corresponds to SPATIAL_1
                         point_y = as.numeric(row["x"]), # Corresponds to SPATIAL_2
                         polygon_df = points_df)
      })
      img_to_keep <- pre_filtered_img[inside_flags_img, ]
    } else {
      img_to_keep <- pre_filtered_img[FALSE, ] # Return empty dataframe structure
    }


  } else {
    stop("Invalid coordinate type. Must be 'rectangle' or 'polygon'.")
  }

  # --- 3. Update the b2c object ---
  cells_to_keep_rownames <- cells_to_keep_df$original_rownames

  # Filter Seurat objects (post and pre)
  b2c$post <- b2c$post[, colnames(b2c$post) %in% cells_to_keep_rownames]
  # Ensure the metadata column exists before filtering pre
  if (label.id %in% colnames(b2c$pre@meta.data)) {
    pre_labels_to_keep <- as.character(b2c$pre@meta.data[, label.id]) %in% cells_to_keep_rownames
    b2c$pre <- b2c$pre[, pre_labels_to_keep]
  } else {
    warning(paste("Label ID '", label.id, "' not found in b2c$pre@meta.data. Skipping filtering for b2c$pre.", sep=""))
  }


  # Update the image data frame
  b2c$img <- img_to_keep

  return(b2c)
}
