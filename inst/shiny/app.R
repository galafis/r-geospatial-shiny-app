#' Geospatial Shiny App
#'
#' @description A Shiny application for interactive geospatial data analysis and visualization.
#'
#' @import shiny
#' @import shinydashboard
#' @import leaflet
#' @import sf
#' @import dplyr
#' @import ggplot2
#' @import plotly
#' @import DT
#' @import RColorBrewer
#'
#' @author Gabriel Demetrios Lafis
#' @date 2023-06-01

# Load required packages
library(shiny)
library(shinydashboard)
library(leaflet)
library(sf)
library(dplyr)
library(ggplot2)
library(plotly)
library(DT)
library(RColorBrewer)
library(rgeospatialshiny)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "Geospatial Analysis"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Data Explorer", tabName = "data", icon = icon("table")),
      menuItem("Spatial Analysis", tabName = "analysis", icon = icon("chart-bar")),
      menuItem("Settings", tabName = "settings", icon = icon("cog"))
    ),
    
    # Data input section
    conditionalPanel(
      condition = "input.sidebar == 'dashboard' || input.sidebar == 'data' || input.sidebar == 'analysis'",
      h4("Data Input"),
      fileInput("file_upload", "Upload Spatial Data", accept = c(".shp", ".geojson", ".gpkg")),
      selectInput("demo_data", "Or Use Demo Data", choices = c("None", "US Counties", "World Countries", "Cities")),
      hr(),
      h4("Map Options"),
      selectInput("basemap", "Basemap", choices = c("OpenStreetMap" = "osm", 
                                                   "Satellite" = "satellite", 
                                                   "Terrain" = "terrain", 
                                                   "Toner" = "toner")),
      selectInput("color_var", "Color Variable", choices = c("None")),
      selectInput("color_palette", "Color Palette", 
                 choices = c("Viridis" = "viridis", 
                            "Blues" = "Blues", 
                            "Reds" = "Reds", 
                            "Greens" = "Greens", 
                            "Spectral" = "Spectral", 
                            "RdYlBu" = "RdYlBu"))
    )
  ),
  
  dashboardBody(
    tabItems(
      # Dashboard tab
      tabItem(tabName = "dashboard",
              fluidRow(
                box(width = 9, height = 600,
                    leafletOutput("main_map", height = 580)
                ),
                box(width = 3,
                    h4("Summary Statistics"),
                    tableOutput("summary_stats"),
                    hr(),
                    h4("Quick Analysis"),
                    selectInput("quick_analysis", "Analysis Type", 
                               choices = c("None", "Spatial Autocorrelation", "Hotspot Analysis", "Point Pattern")),
                    actionButton("run_analysis", "Run Analysis", class = "btn-primary"),
                    hr(),
                    plotOutput("quick_plot", height = 250)
                )
              ),
              fluidRow(
                box(width = 6,
                    h4("Attribute Distribution"),
                    selectInput("dist_var", "Variable", choices = c("None")),
                    plotlyOutput("dist_plot", height = 250)
                ),
                box(width = 6,
                    h4("Spatial Relationship"),
                    selectInput("x_var", "X Variable", choices = c("None")),
                    selectInput("y_var", "Y Variable", choices = c("None")),
                    plotlyOutput("scatter_plot", height = 250)
                )
              )
      ),
      
      # Data Explorer tab
      tabItem(tabName = "data",
              fluidRow(
                box(width = 12,
                    h4("Data Table"),
                    DTOutput("data_table")
                )
              ),
              fluidRow(
                box(width = 6,
                    h4("Data Summary"),
                    verbatimTextOutput("data_summary")
                ),
                box(width = 6,
                    h4("Geometry Information"),
                    verbatimTextOutput("geometry_info")
                )
              )
      ),
      
      # Spatial Analysis tab
      tabItem(tabName = "analysis",
              fluidRow(
                box(width = 3,
                    h4("Analysis Tools"),
                    selectInput("analysis_type", "Analysis Type", 
                               choices = c("Spatial Autocorrelation", "Hotspot Analysis", 
                                          "Buffer Analysis", "Spatial Join", 
                                          "Point Pattern Analysis", "Spatial Interpolation")),
                    conditionalPanel(
                      condition = "input.analysis_type == 'Spatial Autocorrelation' || input.analysis_type == 'Hotspot Analysis'",
                      selectInput("analysis_var", "Variable", choices = c("None")),
                      selectInput("weights_type", "Spatial Weights", 
                                 choices = c("Contiguity" = "contiguity", 
                                            "Distance" = "distance", 
                                            "K-Nearest Neighbors" = "knn"))
                    ),
                    conditionalPanel(
                      condition = "input.analysis_type == 'Buffer Analysis'",
                      numericInput("buffer_distance", "Buffer Distance", value = 1000),
                      checkboxInput("dissolve_buffer", "Dissolve Buffer", value = FALSE)
                    ),
                    conditionalPanel(
                      condition = "input.analysis_type == 'Spatial Join'",
                      fileInput("join_file", "Upload Join Data", accept = c(".shp", ".geojson", ".gpkg")),
                      selectInput("join_type", "Join Type", 
                                 choices = c("Intersects" = "intersects", 
                                            "Contains" = "contains", 
                                            "Within" = "within"))
                    ),
                    conditionalPanel(
                      condition = "input.analysis_type == 'Point Pattern Analysis'",
                      selectInput("ppa_method", "Method", 
                                 choices = c("Kernel Density" = "density", 
                                            "K-Function" = "kfunction", 
                                            "L-Function" = "lfunction")),
                      numericInput("ppa_nsim", "Simulations", value = 39)
                    ),
                    conditionalPanel(
                      condition = "input.analysis_type == 'Spatial Interpolation'",
                      selectInput("interp_var", "Variable", choices = c("None")),
                      selectInput("interp_method", "Method", 
                                 choices = c("IDW" = "idw", "Kriging" = "kriging")),
                      numericInput("interp_resolution", "Resolution", value = 100)
                    ),
                    actionButton("run_spatial_analysis", "Run Analysis", class = "btn-primary")
                ),
                box(width = 9,
                    h4("Analysis Results"),
                    tabsetPanel(id = "analysis_tabs",
                               tabPanel("Map", leafletOutput("analysis_map", height = 500)),
                               tabPanel("Plot", plotlyOutput("analysis_plot", height = 500)),
                               tabPanel("Summary", verbatimTextOutput("analysis_summary"))
                    )
                )
              )
      ),
      
      # Settings tab
      tabItem(tabName = "settings",
              fluidRow(
                box(width = 6,
                    h4("Map Settings"),
                    selectInput("map_projection", "Map Projection", 
                               choices = c("WGS84" = 4326, 
                                          "Web Mercator" = 3857, 
                                          "Equal Area" = 6933)),
                    sliderInput("map_opacity", "Layer Opacity", min = 0, max = 1, value = 0.7, step = 0.1),
                    checkboxInput("show_labels", "Show Labels", value = FALSE),
                    selectInput("label_field", "Label Field", choices = c("None"))
                ),
                box(width = 6,
                    h4("Export Options"),
                    selectInput("export_format", "Export Format", 
                               choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg")),
                    numericInput("export_width", "Width (px)", value = 800),
                    numericInput("export_height", "Height (px)", value = 600),
                    numericInput("export_dpi", "DPI", value = 300),
                    downloadButton("export_map", "Export Map"),
                    hr(),
                    downloadButton("export_data", "Export Data")
                )
              ),
              fluidRow(
                box(width = 12,
                    h4("About"),
                    p("R Geospatial Shiny App for Interactive Spatial Data Analysis"),
                    p("Version: 0.1.0"),
                    p("Author: Gabriel Demetrios Lafis"),
                    p("This application provides tools for loading, processing, analyzing, and visualizing spatial data with a focus on interactive maps, spatial statistics, and geoprocessing operations."),
                    p("Features include choropleth maps, point pattern analysis, spatial autocorrelation, and custom spatial queries through an intuitive user interface.")
                )
              )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  # Reactive values
  rv <- reactiveValues(
    data = NULL,
    join_data = NULL,
    analysis_result = NULL
  )
  
  # Load data
  observeEvent(input$file_upload, {
    req(input$file_upload)
    
    # Get file extension
    ext <- tools::file_ext(input$file_upload$name)
    
    # Load data based on file type
    tryCatch({
      if (ext == "shp") {
        # Create temporary directory for shapefile
        temp_dir <- tempdir()
        unzip(input$file_upload$datapath, exdir = temp_dir)
        shp_file <- list.files(temp_dir, pattern = "\\.shp$", full.names = TRUE)[1]
        rv$data <- sf::st_read(shp_file)
      } else if (ext == "geojson") {
        rv$data <- sf::st_read(input$file_upload$datapath)
      } else if (ext == "gpkg") {
        rv$data <- sf::st_read(input$file_upload$datapath)
      } else {
        showNotification("Unsupported file format", type = "error")
        return(NULL)
      }
      
      # Update UI with data columns
      updateSelectInput(session, "color_var", choices = c("None", names(rv$data)))
      updateSelectInput(session, "dist_var", choices = c("None", names(rv$data)))
      updateSelectInput(session, "x_var", choices = c("None", names(rv$data)))
      updateSelectInput(session, "y_var", choices = c("None", names(rv$data)))
      updateSelectInput(session, "analysis_var", choices = c("None", names(rv$data)))
      updateSelectInput(session, "interp_var", choices = c("None", names(rv$data)))
      updateSelectInput(session, "label_field", choices = c("None", names(rv$data)))
      
      showNotification("Data loaded successfully", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading data:", e$message), type = "error")
    })
  })
  
  # Load demo data
  observeEvent(input$demo_data, {
    req(input$demo_data != "None")
    
    if (input$demo_data == "US Counties") {
      rv$data <- sf::st_read(system.file("extdata", "us_counties.geojson", package = "rgeospatialshiny"))
    } else if (input$demo_data == "World Countries") {
      rv$data <- sf::st_read(system.file("extdata", "world_countries.geojson", package = "rgeospatialshiny"))
    } else if (input$demo_data == "Cities") {
      rv$data <- sf::st_read(system.file("extdata", "cities.geojson", package = "rgeospatialshiny"))
    }
    
    # Update UI with data columns
    updateSelectInput(session, "color_var", choices = c("None", names(rv$data)))
    updateSelectInput(session, "dist_var", choices = c("None", names(rv$data)))
    updateSelectInput(session, "x_var", choices = c("None", names(rv$data)))
    updateSelectInput(session, "y_var", choices = c("None", names(rv$data)))
    updateSelectInput(session, "analysis_var", choices = c("None", names(rv$data)))
    updateSelectInput(session, "interp_var", choices = c("None", names(rv$data)))
    updateSelectInput(session, "label_field", choices = c("None", names(rv$data)))
    
    showNotification("Demo data loaded successfully", type = "message")
  })
  
  # Load join data
  observeEvent(input$join_file, {
    req(input$join_file)
    
    # Get file extension
    ext <- tools::file_ext(input$join_file$name)
    
    # Load data based on file type
    tryCatch({
      if (ext == "shp") {
        # Create temporary directory for shapefile
        temp_dir <- tempdir()
        unzip(input$join_file$datapath, exdir = temp_dir)
        shp_file <- list.files(temp_dir, pattern = "\\.shp$", full.names = TRUE)[1]
        rv$join_data <- sf::st_read(shp_file)
      } else if (ext == "geojson") {
        rv$join_data <- sf::st_read(input$join_file$datapath)
      } else if (ext == "gpkg") {
        rv$join_data <- sf::st_read(input$join_file$datapath)
      } else {
        showNotification("Unsupported file format", type = "error")
        return(NULL)
      }
      
      showNotification("Join data loaded successfully", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading join data:", e$message), type = "error")
    })
  })
  
  # Main map
  output$main_map <- renderLeaflet({
    req(rv$data)
    
    # Get color variable
    color_var <- if (input$color_var != "None") input$color_var else NULL
    
    # Get labels
    labels <- if (input$show_labels && input$label_field != "None") input$label_field else NULL
    
    # Create map
    create_leaflet_map(
      data = rv$data,
      variable = color_var,
      palette = input$color_palette,
      basemap = input$basemap,
      title = "Main Map",
      labels = labels
    )
  })
  
  # Summary statistics
  output$summary_stats <- renderTable({
    req(rv$data)
    
    # Get numeric columns
    numeric_cols <- sapply(rv$data, is.numeric)
    numeric_data <- sf::st_drop_geometry(rv$data[, numeric_cols, drop = FALSE])
    
    # Calculate summary statistics
    if (ncol(numeric_data) > 0) {
      stats <- data.frame(
        Variable = names(numeric_data),
        Min = sapply(numeric_data, min, na.rm = TRUE),
        Mean = sapply(numeric_data, mean, na.rm = TRUE),
        Median = sapply(numeric_data, median, na.rm = TRUE),
        Max = sapply(numeric_data, max, na.rm = TRUE),
        SD = sapply(numeric_data, sd, na.rm = TRUE)
      )
      
      # Round numeric columns
      stats[, -1] <- round(stats[, -1], 2)
      
      return(stats)
    } else {
      return(data.frame(Message = "No numeric variables found"))
    }
  })
  
  # Distribution plot
  output$dist_plot <- renderPlotly({
    req(rv$data, input$dist_var != "None")
    
    # Get variable
    var <- rv$data[[input$dist_var]]
    
    # Create plot based on variable type
    if (is.numeric(var)) {
      # Histogram for numeric variables
      p <- ggplot(rv$data, aes(x = .data[[input$dist_var]])) +
        geom_histogram(fill = "#3388ff", color = "#333333", bins = 30) +
        theme_minimal() +
        labs(title = paste("Distribution of", input$dist_var),
             x = input$dist_var, y = "Count")
    } else {
      # Bar chart for categorical variables
      p <- ggplot(rv$data, aes(x = .data[[input$dist_var]])) +
        geom_bar(fill = "#3388ff", color = "#333333") +
        theme_minimal() +
        labs(title = paste("Distribution of", input$dist_var),
             x = input$dist_var, y = "Count") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    ggplotly(p)
  })
  
  # Scatter plot
  output$scatter_plot <- renderPlotly({
    req(rv$data, input$x_var != "None", input$y_var != "None")
    
    # Create scatter plot
    p <- ggplot(rv$data, aes(x = .data[[input$x_var]], y = .data[[input$y_var]])) +
      geom_point(color = "#3388ff", alpha = 0.7) +
      theme_minimal() +
      labs(title = paste(input$y_var, "vs", input$x_var),
           x = input$x_var, y = input$y_var)
    
    ggplotly(p)
  })
  
  # Data table
  output$data_table <- renderDT({
    req(rv$data)
    
    # Drop geometry column for display
    data_table <- sf::st_drop_geometry(rv$data)
    
    # Create data table
    DT::datatable(data_table, options = list(
      pageLength = 10,
      scrollX = TRUE,
      autoWidth = TRUE
    ))
  })
  
  # Data summary
  output$data_summary <- renderPrint({
    req(rv$data)
    
    # Print summary
    summary(sf::st_drop_geometry(rv$data))
  })
  
  # Geometry information
  output$geometry_info <- renderPrint({
    req(rv$data)
    
    # Get geometry type
    geom_type <- sf::st_geometry_type(rv$data, by_geometry = FALSE)
    
    # Get CRS
    crs <- sf::st_crs(rv$data)
    
    # Get bounding box
    bbox <- sf::st_bbox(rv$data)
    
    # Print information
    cat("Geometry Type:", geom_type, "\n\n")
    cat("CRS:", crs$input, "\n\n")
    cat("Bounding Box:\n")
    cat("  xmin:", bbox["xmin"], "\n")
    cat("  ymin:", bbox["ymin"], "\n")
    cat("  xmax:", bbox["xmax"], "\n")
    cat("  ymax:", bbox["ymax"], "\n\n")
    cat("Number of Features:", nrow(rv$data), "\n")
  })
  
  # Quick analysis
  observeEvent(input$run_analysis, {
    req(rv$data, input$quick_analysis != "None")
    
    if (input$quick_analysis == "Spatial Autocorrelation") {
      req(input$color_var != "None")
      
      # Create spatial weights
      weights <- spatial_weights(rv$data, "contiguity")
      
      # Calculate global Moran's I
      result <- global_autocorrelation(rv$data, input$color_var, weights, "moran")
      
      # Create plot
      output$quick_plot <- renderPlot({
        # Create Moran scatter plot
        x <- scale(rv$data[[input$color_var]])
        lag_x <- spdep::lag.listw(weights, x)
        
        plot(x, lag_x, pch = 19, col = "#3388ff",
             xlab = "Standardized Variable", ylab = "Spatial Lag",
             main = paste("Moran's I =", round(result$estimate[1], 3),
                         "\np-value =", round(result$p.value, 3)))
        abline(lm(lag_x ~ x), col = "red", lwd = 2)
        abline(h = 0, v = 0, lty = 2, col = "gray")
      })
      
    } else if (input$quick_analysis == "Hotspot Analysis") {
      req(input$color_var != "None")
      
      # Create spatial weights
      weights <- spatial_weights(rv$data, "contiguity")
      
      # Calculate local G*
      result <- local_autocorrelation(rv$data, input$color_var, weights, "getis_ord")
      
      # Update map
      output$main_map <- renderLeaflet({
        create_leaflet_map(
          data = result,
          variable = "hotspot_category",
          palette = "RdBu",
          reverse = TRUE,
          basemap = input$basemap,
          title = "Hotspot Analysis"
        )
      })
      
      # Create plot
      output$quick_plot <- renderPlot({
        # Create histogram of local G* values
        hist(result$local_g, breaks = 30, col = "#3388ff", border = "#333333",
             xlab = "Local G* Statistic", ylab = "Frequency",
             main = "Distribution of Local G* Values")
        abline(v = c(-2.58, -1.96, 1.96, 2.58), col = "red", lty = 2)
      })
      
    } else if (input$quick_analysis == "Point Pattern") {
      # Check if data is points
      if (!grepl("POINT", sf::st_geometry_type(rv$data, by_geometry = FALSE))) {
        showNotification("Point pattern analysis requires point data", type = "error")
        return(NULL)
      }
      
      # Perform K-function analysis
      result <- point_pattern_analysis(rv$data, method = "kfunction", nsim = 39)
      
      # Create plot
      output$quick_plot <- renderPlot({
        # Plot K-function with simulation envelope
        plot(result$envelope, main = "K-function with Simulation Envelope",
             xlab = "Distance", ylab = "K(r)")
      })
    }
  })
  
  # Spatial analysis
  observeEvent(input$run_spatial_analysis, {
    req(rv$data, input$analysis_type)
    
    if (input$analysis_type == "Spatial Autocorrelation") {
      req(input$analysis_var != "None")
      
      # Create spatial weights
      weights <- spatial_weights(rv$data, input$weights_type)
      
      # Calculate local Moran's I
      result <- local_autocorrelation(rv$data, input$analysis_var, weights, "moran")
      rv$analysis_result <- result
      
      # Update map
      output$analysis_map <- renderLeaflet({
        create_leaflet_map(
          data = result,
          variable = "lisa_category",
          palette = "Set1",
          basemap = input$basemap,
          title = "Local Moran's I"
        )
      })
      
      # Create plot
      output$analysis_plot <- renderPlotly({
        # Create Moran scatter plot
        x <- scale(rv$data[[input$analysis_var]])
        lag_x <- spdep::lag.listw(weights, x)
        
        p <- ggplot(data.frame(x = x, lag_x = lag_x, category = result$lisa_category),
                   aes(x = x, y = lag_x, color = category)) +
          geom_point(alpha = 0.7) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
          geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
          geom_smooth(method = "lm", se = FALSE, color = "red") +
          theme_minimal() +
          labs(title = "Moran Scatter Plot",
               x = "Standardized Variable",
               y = "Spatial Lag",
               color = "LISA Category")
        
        ggplotly(p)
      })
      
      # Create summary
      output$analysis_summary <- renderPrint({
        # Calculate global Moran's I
        global_result <- global_autocorrelation(rv$data, input$analysis_var, weights, "moran")
        
        # Print results
        cat("Global Moran's I Analysis\n")
        cat("-------------------------\n")
        cat("Moran's I statistic:", round(global_result$estimate[1], 3), "\n")
        cat("Expected value:", round(global_result$estimate[2], 3), "\n")
        cat("Variance:", round(global_result$estimate[3], 3), "\n")
        cat("Standard deviate:", round(global_result$statistic, 3), "\n")
        cat("p-value:", round(global_result$p.value, 5), "\n\n")
        
        cat("Local Moran's I Summary\n")
        cat("----------------------\n")
        cat("Number of features:", nrow(result), "\n")
        cat("LISA Categories:\n")
        print(table(result$lisa_category))
      })
      
    } else if (input$analysis_type == "Hotspot Analysis") {
      req(input$analysis_var != "None")
      
      # Create spatial weights
      weights <- spatial_weights(rv$data, input$weights_type)
      
      # Calculate local G*
      result <- local_autocorrelation(rv$data, input$analysis_var, weights, "getis_ord")
      rv$analysis_result <- result
      
      # Update map
      output$analysis_map <- renderLeaflet({
        create_leaflet_map(
          data = result,
          variable = "hotspot_category",
          palette = "RdBu",
          reverse = TRUE,
          basemap = input$basemap,
          title = "Hotspot Analysis"
        )
      })
      
      # Create plot
      output$analysis_plot <- renderPlotly({
        # Create histogram of local G* values
        p <- ggplot(result, aes(x = local_g)) +
          geom_histogram(bins = 30, fill = "#3388ff", color = "#333333") +
          geom_vline(xintercept = c(-2.58, -1.96, 1.96, 2.58), 
                    linetype = "dashed", color = "red") +
          theme_minimal() +
          labs(title = "Distribution of Local G* Values",
               x = "Local G* Statistic",
               y = "Frequency")
        
        ggplotly(p)
      })
      
      # Create summary
      output$analysis_summary <- renderPrint({
        # Print results
        cat("Getis-Ord G* Hotspot Analysis\n")
        cat("----------------------------\n")
        cat("Number of features:", nrow(result), "\n")
        cat("Hotspot Categories:\n")
        print(table(result$hotspot_category))
      })
      
    } else if (input$analysis_type == "Buffer Analysis") {
      # Create buffer
      result <- create_buffer(rv$data, input$buffer_distance, input$dissolve_buffer)
      rv$analysis_result <- result
      
      # Update map
      output$analysis_map <- renderLeaflet({
        # Create map with original data and buffer
        map <- leaflet() %>%
          add_basemap(input$basemap) %>%
          leaflet::addPolygons(
            data = result,
            fillColor = "#3388ff",
            fillOpacity = 0.5,
            weight = 1,
            color = "#333333",
            group = "Buffer"
          )
        
        # Add original data
        if (grepl("POLYGON", sf::st_geometry_type(rv$data, by_geometry = FALSE))) {
          map <- map %>%
            leaflet::addPolygons(
              data = rv$data,
              fillColor = "#ff3333",
              fillOpacity = 0.7,
              weight = 1,
              color = "#333333",
              group = "Original"
            )
        } else if (grepl("POINT", sf::st_geometry_type(rv$data, by_geometry = FALSE))) {
          map <- map %>%
            leaflet::addCircleMarkers(
              data = rv$data,
              radius = 5,
              fillColor = "#ff3333",
              fillOpacity = 0.7,
              weight = 1,
              color = "#333333",
              group = "Original"
            )
        } else if (grepl("LINE", sf::st_geometry_type(rv$data, by_geometry = FALSE))) {
          map <- map %>%
            leaflet::addPolylines(
              data = rv$data,
              color = "#ff3333",
              weight = 2,
              opacity = 0.7,
              group = "Original"
            )
        }
        
        # Add layers control
        map %>%
          leaflet::addLayersControl(
            overlayGroups = c("Original", "Buffer"),
            options = leaflet::layersControlOptions(collapsed = FALSE)
          )
      })
      
      # Create summary
      output$analysis_summary <- renderPrint({
        # Calculate area
        area <- sf::st_area(result)
        
        # Print results
        cat("Buffer Analysis\n")
        cat("--------------\n")
        cat("Buffer Distance:", input$buffer_distance, "\n")
        cat("Dissolve Buffer:", input$dissolve_buffer, "\n")
        cat("Number of Features:", nrow(result), "\n")
        cat("Total Area:", format(sum(area), scientific = FALSE), "square units\n")
      })
      
    } else if (input$analysis_type == "Spatial Join") {
      req(rv$join_data)
      
      # Perform spatial join
      result <- spatial_join(rv$data, rv$join_data, input$join_type)
      rv$analysis_result <- result
      
      # Update map
      output$analysis_map <- renderLeaflet({
        create_leaflet_map(
          data = result,
          variable = names(result)[1],
          palette = input$color_palette,
          basemap = input$basemap,
          title = "Spatial Join Result"
        )
      })
      
      # Create data table
      output$analysis_plot <- renderPlotly({
        # Create bar chart of join counts
        join_counts <- table(result[[names(result)[1]]])
        
        p <- ggplot(data.frame(category = names(join_counts), count = as.numeric(join_counts)),
                   aes(x = category, y = count)) +
          geom_bar(stat = "identity", fill = "#3388ff", color = "#333333") +
          theme_minimal() +
          labs(title = "Join Results",
               x = names(result)[1],
               y = "Count") +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        ggplotly(p)
      })
      
      # Create summary
      output$analysis_summary <- renderPrint({
        # Print results
        cat("Spatial Join Analysis\n")
        cat("-------------------\n")
        cat("Join Type:", input$join_type, "\n")
        cat("Number of Features (Original):", nrow(rv$data), "\n")
        cat("Number of Features (Join):", nrow(rv$join_data), "\n")
        cat("Number of Features (Result):", nrow(result), "\n\n")
        
        cat("Joined Attributes:\n")
        print(names(rv$join_data))
      })
      
    } else if (input$analysis_type == "Point Pattern Analysis") {
      # Check if data is points
      if (!grepl("POINT", sf::st_geometry_type(rv$data, by_geometry = FALSE))) {
        showNotification("Point pattern analysis requires point data", type = "error")
        return(NULL)
      }
      
      # Perform point pattern analysis
      result <- point_pattern_analysis(rv$data, method = input$ppa_method, nsim = input$ppa_nsim)
      rv$analysis_result <- result
      
      if (input$ppa_method == "density") {
        # Convert density to raster for mapping
        density_raster <- raster::raster(result$density)
        
        # Update map
        output$analysis_map <- renderLeaflet({
          # Create map with density raster and points
          leaflet() %>%
            add_basemap(input$basemap) %>%
            leaflet::addRasterImage(
              density_raster,
              colors = colorRampPalette(c("blue", "cyan", "green", "yellow", "red"))(100),
              opacity = 0.7
            ) %>%
            leaflet::addCircleMarkers(
              data = rv$data,
              radius = 3,
              color = "#333333",
              fillColor = "#ffffff",
              fillOpacity = 0.7,
              weight = 1
            )
        })
        
        # Create plot
        output$analysis_plot <- renderPlot({
          # Plot density
          plot(result$density, main = "Kernel Density Estimation",
               col = colorRampPalette(c("blue", "cyan", "green", "yellow", "red"))(100))
          points(sf::st_coordinates(rv$data), pch = 20, cex = 0.5)
        })
        
      } else if (input$ppa_method == "kfunction") {
        # Create plot
        output$analysis_plot <- renderPlot({
          # Plot K-function with simulation envelope
          plot(result$envelope, main = "K-function with Simulation Envelope",
               xlab = "Distance", ylab = "K(r)")
        })
        
      } else if (input$ppa_method == "lfunction") {
        # Create plot
        output$analysis_plot <- renderPlot({
          # Plot L-function with simulation envelope
          plot(result$envelope, main = "L-function with Simulation Envelope",
               xlab = "Distance", ylab = "L(r)")
        })
      }
      
      # Create summary
      output$analysis_summary <- renderPrint({
        # Print results
        cat("Point Pattern Analysis\n")
        cat("--------------------\n")
        cat("Method:", input$ppa_method, "\n")
        cat("Number of Points:", nrow(rv$data), "\n")
        
        if (input$ppa_method == "density") {
          cat("Bandwidth:", result$bandwidth, "\n")
        } else {
          cat("Number of Simulations:", input$ppa_nsim, "\n")
        }
      })
      
    } else if (input$analysis_type == "Spatial Interpolation") {
      req(input$interp_var != "None")
      
      # Check if data is points
      if (!grepl("POINT", sf::st_geometry_type(rv$data, by_geometry = FALSE))) {
        showNotification("Spatial interpolation requires point data", type = "error")
        return(NULL)
      }
      
      # Create grid for interpolation
      bbox <- sf::st_bbox(rv$data)
      grid_raster <- raster::raster(
        xmn = bbox["xmin"], xmx = bbox["xmax"],
        ymn = bbox["ymin"], ymx = bbox["ymax"],
        resolution = c(
          (bbox["xmax"] - bbox["xmin"]) / input$interp_resolution,
          (bbox["ymax"] - bbox["ymin"]) / input$interp_resolution
        )
      )
      
      # Perform interpolation
      result <- spatial_interpolation(
        rv$data, input$interp_var, grid_raster,
        method = input$interp_method,
        params = list(power = 2)
      )
      rv$analysis_result <- result
      
      # Update map
      output$analysis_map <- renderLeaflet({
        # Create map with interpolation raster and points
        leaflet() %>%
          add_basemap(input$basemap) %>%
          leaflet::addRasterImage(
            result,
            colors = colorRampPalette(c("blue", "cyan", "green", "yellow", "red"))(100),
            opacity = 0.7
          ) %>%
          leaflet::addCircleMarkers(
            data = rv$data,
            radius = 3,
            color = "#333333",
            fillColor = "#ffffff",
            fillOpacity = 0.7,
            weight = 1
          )
      })
      
      # Create plot
      output$analysis_plot <- renderPlot({
        # Plot interpolation
        plot(result, main = paste("Spatial Interpolation of", input$interp_var),
             col = colorRampPalette(c("blue", "cyan", "green", "yellow", "red"))(100))
        points(sf::st_coordinates(rv$data), pch = 20, cex = 0.5)
      })
      
      # Create summary
      output$analysis_summary <- renderPrint({
        # Print results
        cat("Spatial Interpolation\n")
        cat("-------------------\n")
        cat("Variable:", input$interp_var, "\n")
        cat("Method:", input$interp_method, "\n")
        cat("Resolution:", input$interp_resolution, "\n")
        cat("Number of Points:", nrow(rv$data), "\n")
        cat("Raster Resolution:", raster::res(result), "\n")
        cat("Raster Dimensions:", dim(result), "\n")
        cat("Value Range:", raster::minValue(result), "to", raster::maxValue(result), "\n")
      })
    }
  })
  
  # Export map
  output$export_map <- downloadHandler(
    filename = function() {
      paste("map-", Sys.Date(), ".", input$export_format, sep = "")
    },
    content = function(file) {
      # Save map based on format
      if (input$export_format == "png") {
        mapview::mapshot(leafletProxy("main_map"), file = file,
                        width = input$export_width, height = input$export_height,
                        res = input$export_dpi)
      } else if (input$export_format == "pdf") {
        mapview::mapshot(leafletProxy("main_map"), file = file,
                        width = input$export_width, height = input$export_height)
      } else if (input$export_format == "svg") {
        mapview::mapshot(leafletProxy("main_map"), file = file,
                        width = input$export_width, height = input$export_height)
      }
    }
  )
  
  # Export data
  output$export_data <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".geojson", sep = "")
    },
    content = function(file) {
      # Export data as GeoJSON
      sf::st_write(rv$data, file, delete_dsn = TRUE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)

