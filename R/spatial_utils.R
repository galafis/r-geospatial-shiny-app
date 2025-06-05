#' Spatial Utility Functions
#'
#' @description A collection of utility functions for spatial data processing and analysis.
#'
#' @importFrom sf st_read st_write st_transform st_buffer st_intersection st_union st_area
#' @importFrom sp SpatialPointsDataFrame CRS
#' @importFrom raster raster extent crop mask
#' @importFrom dplyr %>% filter select mutate group_by summarise
#' @importFrom magrittr %<>%
#'
#' @name spatial_utils
NULL

#' Load spatial data from various formats
#'
#' @param file_path Path to the spatial data file
#' @param layer Layer name (for multi-layer data sources)
#' @param crs Coordinate reference system to transform to (EPSG code or proj4string)
#'
#' @return An sf object containing the spatial data
#' @export
#'
#' @examples
#' \dontrun{
#' # Load shapefile
#' counties <- load_spatial_data("data/counties.shp")
#'
#' # Load GeoJSON and transform to WGS84
#' points <- load_spatial_data("data/points.geojson", crs = 4326)
#' }
load_spatial_data <- function(file_path, layer = NULL, crs = NULL) {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  
  # Read spatial data
  if (is.null(layer)) {
    data <- sf::st_read(file_path, quiet = TRUE)
  } else {
    data <- sf::st_read(file_path, layer = layer, quiet = TRUE)
  }
  
  # Transform CRS if specified
  if (!is.null(crs)) {
    data <- sf::st_transform(data, crs)
  }
  
  return(data)
}

#' Create a buffer around spatial features
#'
#' @param data An sf object
#' @param distance Buffer distance in the units of the data's CRS
#' @param dissolve Whether to dissolve the buffers into a single feature
#'
#' @return An sf object with buffered geometries
#' @export
#'
#' @examples
#' \dontrun{
#' # Create 1000m buffer around points
#' buffered_points <- create_buffer(points, 1000)
#'
#' # Create dissolved buffer
#' study_area <- create_buffer(points, 1000, dissolve = TRUE)
#' }
create_buffer <- function(data, distance, dissolve = FALSE) {
  # Create buffer
  buffered <- sf::st_buffer(data, dist = distance)
  
  # Dissolve if requested
  if (dissolve && nrow(buffered) > 1) {
    buffered <- sf::st_union(buffered) %>% sf::st_sf()
  }
  
  return(buffered)
}

#' Perform spatial intersection between two sf objects
#'
#' @param x First sf object
#' @param y Second sf object
#' @param keep_all Whether to keep all attributes from both datasets
#'
#' @return An sf object with the intersection
#' @export
#'
#' @examples
#' \dontrun{
#' # Intersect points with polygons
#' points_in_counties <- spatial_intersection(points, counties)
#' }
spatial_intersection <- function(x, y, keep_all = TRUE) {
  # Check if both inputs are sf objects
  if (!inherits(x, "sf") || !inherits(y, "sf")) {
    stop("Both x and y must be sf objects")
  }
  
  # Ensure both datasets have the same CRS
  if (sf::st_crs(x) != sf::st_crs(y)) {
    y <- sf::st_transform(y, sf::st_crs(x))
  }
  
  # Perform intersection
  result <- sf::st_intersection(x, y)
  
  # Filter attributes if requested
  if (!keep_all) {
    result <- result %>% dplyr::select(names(x))
  }
  
  return(result)
}

#' Calculate area of spatial features
#'
#' @param data An sf object with polygon geometries
#' @param units Units for area calculation ("m2", "km2", "ha", "acres")
#'
#' @return The input sf object with an additional column for area
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate area in square kilometers
#' counties_with_area <- calculate_area(counties, "km2")
#' }
calculate_area <- function(data, units = "m2") {
  # Check if input is an sf object
  if (!inherits(data, "sf")) {
    stop("Input must be an sf object")
  }
  
  # Calculate area in square meters
  area_m2 <- sf::st_area(data)
  
  # Convert to requested units
  area <- switch(units,
                "m2" = area_m2,
                "km2" = area_m2 / 1e6,
                "ha" = area_m2 / 10000,
                "acres" = area_m2 / 4046.86,
                stop("Unsupported units: ", units))
  
  # Add area column to data
  data$area <- as.numeric(area)
  
  return(data)
}

#' Convert coordinates to spatial points
#'
#' @param data A data frame with coordinate columns
#' @param x_col Name of the column containing x-coordinates (longitude)
#' @param y_col Name of the column containing y-coordinates (latitude)
#' @param crs Coordinate reference system (default: WGS84, EPSG:4326)
#'
#' @return An sf object with point geometries
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert coordinates to spatial points
#' points_sf <- coords_to_sf(data, "longitude", "latitude")
#' }
coords_to_sf <- function(data, x_col, y_col, crs = 4326) {
  # Check if input is a data frame
  if (!is.data.frame(data)) {
    stop("Input must be a data frame")
  }
  
  # Check if coordinate columns exist
  if (!x_col %in% names(data) || !y_col %in% names(data)) {
    stop("Coordinate columns not found in data")
  }
  
  # Convert to sf object
  sf_data <- sf::st_as_sf(data, coords = c(x_col, y_col), crs = crs)
  
  return(sf_data)
}

#' Clip raster by polygon
#'
#' @param raster_data A raster object
#' @param polygon_data An sf object with polygon geometries
#'
#' @return A clipped raster object
#' @export
#'
#' @examples
#' \dontrun{
#' # Clip raster by county boundary
#' clipped_raster <- clip_raster(elevation, county)
#' }
clip_raster <- function(raster_data, polygon_data) {
  # Check if inputs are of correct class
  if (!inherits(raster_data, "Raster")) {
    stop("raster_data must be a Raster object")
  }
  
  if (!inherits(polygon_data, "sf")) {
    stop("polygon_data must be an sf object")
  }
  
  # Convert sf to sp for raster operations
  polygon_sp <- as(polygon_data, "Spatial")
  
  # Ensure same CRS
  if (raster::compareCRS(raster::crs(raster_data), raster::crs(polygon_sp)) == FALSE) {
    polygon_sp <- sp::spTransform(polygon_sp, raster::crs(raster_data))
  }
  
  # Clip raster
  clipped <- raster::crop(raster_data, raster::extent(polygon_sp))
  clipped <- raster::mask(clipped, polygon_sp)
  
  return(clipped)
}

#' Calculate zonal statistics
#'
#' @param raster_data A raster object
#' @param zone_data An sf object with polygon geometries
#' @param stats Statistics to calculate (mean, sum, min, max, sd)
#'
#' @return The zone_data sf object with additional columns for zonal statistics
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate mean elevation by county
#' counties_elev <- zonal_statistics(elevation, counties, stats = c("mean", "max"))
#' }
zonal_statistics <- function(raster_data, zone_data, stats = c("mean")) {
  # Check if inputs are of correct class
  if (!inherits(raster_data, "Raster")) {
    stop("raster_data must be a Raster object")
  }
  
  if (!inherits(zone_data, "sf")) {
    stop("zone_data must be an sf object")
  }
  
  # Convert sf to sp for raster operations
  zone_sp <- as(zone_data, "Spatial")
  
  # Ensure same CRS
  if (raster::compareCRS(raster::crs(raster_data), raster::crs(zone_sp)) == FALSE) {
    zone_sp <- sp::spTransform(zone_sp, raster::crs(raster_data))
  }
  
  # Calculate zonal statistics
  result <- zone_data
  
  for (stat in stats) {
    if (stat == "mean") {
      values <- raster::extract(raster_data, zone_sp, fun = mean, na.rm = TRUE)
      result$mean <- values
    } else if (stat == "sum") {
      values <- raster::extract(raster_data, zone_sp, fun = sum, na.rm = TRUE)
      result$sum <- values
    } else if (stat == "min") {
      values <- raster::extract(raster_data, zone_sp, fun = min, na.rm = TRUE)
      result$min <- values
    } else if (stat == "max") {
      values <- raster::extract(raster_data, zone_sp, fun = max, na.rm = TRUE)
      result$max <- values
    } else if (stat == "sd") {
      values <- raster::extract(raster_data, zone_sp, fun = sd, na.rm = TRUE)
      result$sd <- values
    } else {
      warning("Unsupported statistic: ", stat)
    }
  }
  
  return(result)
}

