#' Spatial Visualization Functions
#'
#' @description A collection of functions for visualizing spatial data.
#'
#' @importFrom leaflet leaflet addTiles addPolygons addCircleMarkers addLegend addLayersControl
#' @importFrom ggplot2 ggplot aes geom_sf scale_fill_viridis_c theme_minimal labs
#' @importFrom plotly ggplotly
#' @importFrom RColorBrewer brewer.pal
#' @importFrom sf st_centroid st_coordinates
#' @importFrom dplyr %>%
#'
#' @name visualization
NULL

#' Create an interactive map using leaflet
#'
#' @param data An sf object or list of sf objects
#' @param variable Name of the variable to visualize (for choropleth maps)
#' @param palette Color palette to use
#' @param reverse Whether to reverse the color palette
#' @param basemap Basemap to use ("osm", "satellite", "toner", "terrain")
#' @param title Map title
#' @param labels Labels for features (column name or function)
#'
#' @return A leaflet map object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create choropleth map of population density
#' map <- create_leaflet_map(counties, "pop_density", title = "Population Density")
#'
#' # Create map with multiple layers
#' map <- create_leaflet_map(list(counties = counties, cities = cities))
#' }
create_leaflet_map <- function(data, variable = NULL, palette = "viridis", reverse = FALSE,
                              basemap = "osm", title = NULL, labels = NULL) {
  # Initialize leaflet map
  map <- leaflet::leaflet() %>%
    add_basemap(basemap)
  
  # Add title if provided
  if (!is.null(title)) {
    map <- map %>%
      leaflet::addControl(html = paste0("<h4>", title, "</h4>"), position = "topright")
  }
  
  # Process data based on type
  if (inherits(data, "sf")) {
    # Single sf object
    map <- add_sf_to_leaflet(map, data, variable, palette, reverse, labels)
  } else if (is.list(data)) {
    # List of sf objects
    for (i in seq_along(data)) {
      layer_name <- names(data)[i]
      if (is.null(layer_name)) layer_name <- paste("Layer", i)
      
      layer_data <- data[[i]]
      if (!inherits(layer_data, "sf")) {
        warning("Skipping non-sf object in list: ", layer_name)
        next
      }
      
      # Get variable for this layer
      layer_variable <- NULL
      if (!is.null(variable)) {
        if (is.list(variable) && length(variable) >= i) {
          layer_variable <- variable[[i]]
        } else if (is.character(variable) && variable %in% names(layer_data)) {
          layer_variable <- variable
        }
      }
      
      # Add layer to map
      map <- add_sf_to_leaflet(map, layer_data, layer_variable, palette, reverse, labels, group = layer_name)
    }
    
    # Add layers control
    map <- map %>%
      leaflet::addLayersControl(
        overlayGroups = names(data),
        options = leaflet::layersControlOptions(collapsed = FALSE)
      )
  } else {
    stop("data must be an sf object or a list of sf objects")
  }
  
  return(map)
}

#' Add an sf object to a leaflet map
#'
#' @param map A leaflet map object
#' @param data An sf object
#' @param variable Name of the variable to visualize
#' @param palette Color palette to use
#' @param reverse Whether to reverse the color palette
#' @param labels Labels for features
#' @param group Group name for layer control
#'
#' @return A leaflet map object with the sf object added
#' @keywords internal
add_sf_to_leaflet <- function(map, data, variable = NULL, palette = "viridis", reverse = FALSE,
                             labels = NULL, group = NULL) {
  # Determine geometry type
  geom_type <- sf::st_geometry_type(data, by_geometry = FALSE)
  
  # Process labels
  if (is.null(labels)) {
    feature_labels <- NULL
  } else if (is.function(labels)) {
    feature_labels <- labels(data)
  } else if (is.character(labels) && labels %in% names(data)) {
    feature_labels <- data[[labels]]
  } else {
    feature_labels <- NULL
  }
  
  # Add data based on geometry type and whether a variable is specified
  if (grepl("POLYGON", geom_type)) {
    if (!is.null(variable) && variable %in% names(data)) {
      # Choropleth map
      values <- data[[variable]]
      
      # Determine color palette
      if (is.factor(values) || is.character(values)) {
        # Categorical variable
        n_colors <- length(unique(values))
        if (palette == "viridis") {
          pal <- leaflet::colorFactor("viridis", domain = values, reverse = reverse)
        } else {
          pal <- leaflet::colorFactor(palette, domain = values, reverse = reverse)
        }
        legend_title <- variable
      } else {
        # Continuous variable
        if (palette == "viridis") {
          pal <- leaflet::colorNumeric("viridis", domain = values, reverse = reverse)
        } else {
          pal <- leaflet::colorNumeric(palette, domain = values, reverse = reverse)
        }
        legend_title <- variable
      }
      
      # Add polygons with colors
      map <- map %>%
        leaflet::addPolygons(
          data = data,
          fillColor = ~pal(values),
          fillOpacity = 0.7,
          weight = 1,
          color = "#333333",
          opacity = 1,
          label = feature_labels,
          group = group,
          highlightOptions = leaflet::highlightOptions(
            weight = 2,
            color = "#666666",
            fillOpacity = 0.9,
            bringToFront = TRUE
          )
        ) %>%
        leaflet::addLegend(
          position = "bottomright",
          pal = pal,
          values = values,
          title = legend_title,
          opacity = 0.7
        )
    } else {
      # Simple polygon map
      map <- map %>%
        leaflet::addPolygons(
          data = data,
          fillColor = "#3388ff",
          fillOpacity = 0.5,
          weight = 1,
          color = "#333333",
          opacity = 1,
          label = feature_labels,
          group = group,
          highlightOptions = leaflet::highlightOptions(
            weight = 2,
            color = "#666666",
            fillOpacity = 0.9,
            bringToFront = TRUE
          )
        )
    }
  } else if (grepl("POINT", geom_type)) {
    if (!is.null(variable) && variable %in% names(data)) {
      # Proportional symbol map or colored points
      values <- data[[variable]]
      
      # Determine if values should be used for size or color
      if (is.numeric(values) && !is.factor(values)) {
        # Use for both size and color
        if (palette == "viridis") {
          pal <- leaflet::colorNumeric("viridis", domain = values, reverse = reverse)
        } else {
          pal <- leaflet::colorNumeric(palette, domain = values, reverse = reverse)
        }
        
        # Scale radius between 3 and 15 based on values
        radius <- scales::rescale(values, to = c(3, 15))
        
        map <- map %>%
          leaflet::addCircleMarkers(
            data = data,
            radius = radius,
            color = "#333333",
            weight = 1,
            fillColor = ~pal(values),
            fillOpacity = 0.7,
            label = feature_labels,
            group = group,
            popup = ~paste0("<strong>", variable, ":</strong> ", values)
          ) %>%
          leaflet::addLegend(
            position = "bottomright",
            pal = pal,
            values = values,
            title = variable,
            opacity = 0.7
          )
      } else {
        # Categorical variable, use for color only
        if (is.factor(values) || is.character(values)) {
          if (palette == "viridis") {
            pal <- leaflet::colorFactor("viridis", domain = values, reverse = reverse)
          } else {
            pal <- leaflet::colorFactor(palette, domain = values, reverse = reverse)
          }
          
          map <- map %>%
            leaflet::addCircleMarkers(
              data = data,
              radius = 6,
              color = "#333333",
              weight = 1,
              fillColor = ~pal(values),
              fillOpacity = 0.7,
              label = feature_labels,
              group = group,
              popup = ~paste0("<strong>", variable, ":</strong> ", values)
            ) %>%
            leaflet::addLegend(
              position = "bottomright",
              pal = pal,
              values = values,
              title = variable,
              opacity = 0.7
            )
        }
      }
    } else {
      # Simple point map
      map <- map %>%
        leaflet::addCircleMarkers(
          data = data,
          radius = 6,
          color = "#333333",
          weight = 1,
          fillColor = "#3388ff",
          fillOpacity = 0.7,
          label = feature_labels,
          group = group
        )
    }
  } else if (grepl("LINE", geom_type)) {
    if (!is.null(variable) && variable %in% names(data)) {
      # Colored lines
      values <- data[[variable]]
      
      if (is.factor(values) || is.character(values)) {
        # Categorical variable
        if (palette == "viridis") {
          pal <- leaflet::colorFactor("viridis", domain = values, reverse = reverse)
        } else {
          pal <- leaflet::colorFactor(palette, domain = values, reverse = reverse)
        }
      } else {
        # Continuous variable
        if (palette == "viridis") {
          pal <- leaflet::colorNumeric("viridis", domain = values, reverse = reverse)
        } else {
          pal <- leaflet::colorNumeric(palette, domain = values, reverse = reverse)
        }
      }
      
      map <- map %>%
        leaflet::addPolylines(
          data = data,
          color = ~pal(values),
          weight = 3,
          opacity = 0.7,
          label = feature_labels,
          group = group,
          highlightOptions = leaflet::highlightOptions(
            weight = 5,
            color = "#666666",
            opacity = 1,
            bringToFront = TRUE
          )
        ) %>%
        leaflet::addLegend(
          position = "bottomright",
          pal = pal,
          values = values,
          title = variable,
          opacity = 0.7
        )
    } else {
      # Simple line map
      map <- map %>%
        leaflet::addPolylines(
          data = data,
          color = "#3388ff",
          weight = 3,
          opacity = 0.7,
          label = feature_labels,
          group = group,
          highlightOptions = leaflet::highlightOptions(
            weight = 5,
            color = "#666666",
            opacity = 1,
            bringToFront = TRUE
          )
        )
    }
  }
  
  return(map)
}

#' Add a basemap to a leaflet map
#'
#' @param map A leaflet map object
#' @param basemap Basemap to use
#'
#' @return A leaflet map object with the basemap added
#' @keywords internal
add_basemap <- function(map, basemap = "osm") {
  if (basemap == "osm") {
    map <- map %>% leaflet::addTiles()
  } else if (basemap == "satellite") {
    map <- map %>% leaflet::addProviderTiles(leaflet::providers$Esri.WorldImagery)
  } else if (basemap == "toner") {
    map <- map %>% leaflet::addProviderTiles(leaflet::providers$Stamen.Toner)
  } else if (basemap == "terrain") {
    map <- map %>% leaflet::addProviderTiles(leaflet::providers$Stamen.TerrainBackground)
  } else {
    warning("Unsupported basemap: ", basemap, ". Using OpenStreetMap.")
    map <- map %>% leaflet::addTiles()
  }
  
  return(map)
}

#' Create a static map using ggplot2
#'
#' @param data An sf object
#' @param variable Name of the variable to visualize
#' @param palette Color palette to use
#' @param title Map title
#' @param subtitle Map subtitle
#' @param legend_title Legend title
#'
#' @return A ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create choropleth map of population density
#' map <- create_static_map(counties, "pop_density", title = "Population Density")
#' }
create_static_map <- function(data, variable = NULL, palette = "viridis",
                             title = NULL, subtitle = NULL, legend_title = NULL) {
  # Check if input is an sf object
  if (!inherits(data, "sf")) {
    stop("data must be an sf object")
  }
  
  # Create base map
  p <- ggplot2::ggplot()
  
  # Add data based on whether a variable is specified
  if (!is.null(variable) && variable %in% names(data)) {
    # Choropleth map
    p <- p + ggplot2::geom_sf(data = data, ggplot2::aes(fill = .data[[variable]]), color = "#333333", size = 0.2)
    
    # Add color scale based on variable type
    if (is.factor(data[[variable]]) || is.character(data[[variable]])) {
      # Categorical variable
      if (palette == "viridis") {
        p <- p + ggplot2::scale_fill_viridis_d(option = "viridis")
      } else {
        p <- p + ggplot2::scale_fill_brewer(palette = palette)
      }
    } else {
      # Continuous variable
      if (palette == "viridis") {
        p <- p + ggplot2::scale_fill_viridis_c(option = "viridis")
      } else {
        p <- p + ggplot2::scale_fill_distiller(palette = palette)
      }
    }
    
    # Set legend title
    if (is.null(legend_title)) {
      legend_title <- variable
    }
    p <- p + ggplot2::labs(fill = legend_title)
  } else {
    # Simple map
    p <- p + ggplot2::geom_sf(data = data, fill = "#3388ff", color = "#333333", size = 0.2)
  }
  
  # Add titles
  p <- p + ggplot2::labs(
    title = title,
    subtitle = subtitle
  )
  
  # Apply theme
  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  
  return(p)
}

#' Create an interactive plot using plotly
#'
#' @param data A data frame or sf object
#' @param x Name of the x variable
#' @param y Name of the y variable
#' @param color Name of the color variable (optional)
#' @param size Name of the size variable (optional)
#' @param type Plot type ("scatter", "line", "bar", "box")
#' @param title Plot title
#' @param xlab X-axis label
#' @param ylab Y-axis label
#'
#' @return A plotly object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create scatter plot of population vs. area
#' plot <- create_interactive_plot(counties, "area", "population", color = "region")
#' }
create_interactive_plot <- function(data, x, y, color = NULL, size = NULL,
                                   type = "scatter", title = NULL, xlab = NULL, ylab = NULL) {
  # Convert sf object to data frame if necessary
  if (inherits(data, "sf")) {
    data <- sf::st_drop_geometry(data)
  }
  
  # Check if variables exist
  if (!x %in% names(data)) {
    stop("x variable not found in data: ", x)
  }
  if (!y %in% names(data)) {
    stop("y variable not found in data: ", y)
  }
  if (!is.null(color) && !color %in% names(data)) {
    stop("color variable not found in data: ", color)
  }
  if (!is.null(size) && !size %in% names(data)) {
    stop("size variable not found in data: ", size)
  }
  
  # Create base ggplot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x]], y = .data[[y]]))
  
  # Add geometry based on plot type
  if (type == "scatter") {
    if (!is.null(color) && !is.null(size)) {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = .data[[color]], size = .data[[size]]))
    } else if (!is.null(color)) {
      p <- p + ggplot2::geom_point(ggplot2::aes(color = .data[[color]]))
    } else if (!is.null(size)) {
      p <- p + ggplot2::geom_point(ggplot2::aes(size = .data[[size]]))
    } else {
      p <- p + ggplot2::geom_point()
    }
  } else if (type == "line") {
    if (!is.null(color)) {
      p <- p + ggplot2::geom_line(ggplot2::aes(color = .data[[color]]))
    } else {
      p <- p + ggplot2::geom_line()
    }
  } else if (type == "bar") {
    if (!is.null(color)) {
      p <- p + ggplot2::geom_col(ggplot2::aes(fill = .data[[color]]))
    } else {
      p <- p + ggplot2::geom_col()
    }
  } else if (type == "box") {
    if (!is.null(color)) {
      p <- p + ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[color]]))
    } else {
      p <- p + ggplot2::geom_boxplot()
    }
  } else {
    stop("Unsupported plot type: ", type)
  }
  
  # Add labels
  p <- p + ggplot2::labs(
    title = title,
    x = xlab %||% x,
    y = ylab %||% y
  )
  
  # Apply theme
  p <- p + ggplot2::theme_minimal()
  
  # Convert to plotly
  plotly_obj <- plotly::ggplotly(p)
  
  return(plotly_obj)
}

#' Create a map dashboard
#'
#' @param main_map A leaflet map object for the main map
#' @param sidebar_plots A list of plots for the sidebar
#' @param title Dashboard title
#' @param sidebar_width Width of the sidebar (in pixels)
#'
#' @return A dashboard HTML widget
#' @export
#'
#' @examples
#' \dontrun{
#' # Create main map
#' main_map <- create_leaflet_map(counties, "pop_density")
#'
#' # Create sidebar plots
#' plot1 <- create_interactive_plot(counties, "area", "population")
#' plot2 <- create_interactive_plot(counties, "region", "population", type = "box")
#'
#' # Create dashboard
#' dashboard <- create_map_dashboard(main_map, list(plot1, plot2), "County Dashboard")
#' }
create_map_dashboard <- function(main_map, sidebar_plots, title = "Spatial Dashboard", sidebar_width = 400) {
  # Create HTML for dashboard
  html <- htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$style(
        htmltools::HTML(
          paste0("
          body, html {
            margin: 0;
            padding: 0;
            height: 100%;
            font-family: Arial, sans-serif;
          }
          .dashboard-container {
            display: flex;
            height: 100vh;
          }
          .main-map {
            flex: 1;
            height: 100%;
          }
          .sidebar {
            width: ", sidebar_width, "px;
            height: 100%;
            overflow-y: auto;
            background-color: #f8f9fa;
            padding: 10px;
            box-shadow: -2px 0 5px rgba(0,0,0,0.1);
          }
          .dashboard-title {
            text-align: center;
            margin: 10px 0;
            padding: 10px;
            background-color: #4285f4;
            color: white;
            font-size: 18px;
            font-weight: bold;
          }
          .plot-container {
            margin-bottom: 20px;
            background-color: white;
            border-radius: 5px;
            padding: 10px;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
          }
          ")
        )
      )
    ),
    htmltools::tags$body(
      htmltools::tags$div(
        class = "dashboard-container",
        htmltools::tags$div(class = "main-map", main_map),
        htmltools::tags$div(
          class = "sidebar",
          htmltools::tags$div(class = "dashboard-title", title),
          lapply(sidebar_plots, function(plot) {
            htmltools::tags$div(class = "plot-container", plot)
          })
        )
      )
    )
  )
  
  return(html)
}

