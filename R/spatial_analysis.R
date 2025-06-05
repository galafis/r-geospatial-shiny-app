#' Spatial Analysis Functions
#'
#' @description A collection of functions for spatial analysis and statistics.
#'
#' @importFrom sf st_distance st_centroid st_coordinates st_within st_intersects
#' @importFrom spdep poly2nb nb2listw moran.test geary.test localG
#' @importFrom spatstat ppp density.ppp Kest Lest envelope
#' @importFrom dplyr %>% filter select mutate group_by summarise
#'
#' @name spatial_analysis
NULL

#' Calculate distance matrix between spatial features
#'
#' @param x An sf object
#' @param y An sf object (optional, if NULL distances within x are calculated)
#' @param by_element Whether to return distances by element or as a matrix
#'
#' @return A distance matrix or vector of distances
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate distances between cities
#' dist_matrix <- distance_matrix(cities)
#'
#' # Calculate distances from points to nearest city
#' nearest_dist <- distance_matrix(points, cities, by_element = TRUE)
#' }
distance_matrix <- function(x, y = NULL, by_element = FALSE) {
  # Check if input is an sf object
  if (!inherits(x, "sf")) {
    stop("x must be an sf object")
  }
  
  if (!is.null(y) && !inherits(y, "sf")) {
    stop("y must be an sf object")
  }
  
  # Calculate distance matrix
  if (is.null(y)) {
    dist <- sf::st_distance(x)
  } else {
    dist <- sf::st_distance(x, y)
  }
  
  # Return by element if requested
  if (by_element) {
    if (is.null(y)) {
      # For each feature, find minimum distance to any other feature
      diag(dist) <- NA  # Ignore self-distances
      return(apply(dist, 1, min, na.rm = TRUE))
    } else {
      # For each feature in x, find minimum distance to any feature in y
      return(apply(dist, 1, min))
    }
  }
  
  return(dist)
}

#' Find nearest neighbors for spatial features
#'
#' @param x An sf object
#' @param y An sf object to find neighbors from (optional, if NULL neighbors within x are found)
#' @param k Number of nearest neighbors to find
#'
#' @return A list with nearest neighbor indices and distances
#' @export
#'
#' @examples
#' \dontrun{
#' # Find 5 nearest neighbors for each city
#' neighbors <- nearest_neighbors(cities, k = 5)
#'
#' # Find nearest city for each point
#' nearest_city <- nearest_neighbors(points, cities, k = 1)
#' }
nearest_neighbors <- function(x, y = NULL, k = 1) {
  # Check if input is an sf object
  if (!inherits(x, "sf")) {
    stop("x must be an sf object")
  }
  
  if (!is.null(y) && !inherits(y, "sf")) {
    stop("y must be an sf object")
  }
  
  # Calculate distance matrix
  if (is.null(y)) {
    dist <- sf::st_distance(x)
    diag(dist) <- Inf  # Ignore self-distances
  } else {
    dist <- sf::st_distance(x, y)
  }
  
  # Convert to matrix if not already
  dist_matrix <- as.matrix(dist)
  
  # Find k nearest neighbors
  nn_indices <- t(apply(dist_matrix, 1, function(row) {
    order(row)[1:min(k, length(row))]
  }))
  
  # Extract distances
  nn_distances <- matrix(NA, nrow = nrow(dist_matrix), ncol = k)
  for (i in 1:nrow(dist_matrix)) {
    nn_distances[i, ] <- dist_matrix[i, nn_indices[i, ]]
  }
  
  return(list(indices = nn_indices, distances = nn_distances))
}

#' Calculate spatial weights matrix
#'
#' @param data An sf object with polygon geometries
#' @param type Type of spatial weights ("contiguity", "distance", "knn")
#' @param k Number of neighbors (for knn)
#' @param dist_threshold Distance threshold (for distance-based weights)
#' @param style Style of weights ("W", "B", "C", "U", "S")
#'
#' @return A listw object from spdep package
#' @export
#'
#' @examples
#' \dontrun{
#' # Create contiguity-based weights
#' w_contig <- spatial_weights(counties, "contiguity")
#'
#' # Create k-nearest neighbor weights
#' w_knn <- spatial_weights(counties, "knn", k = 5)
#' }
spatial_weights <- function(data, type = "contiguity", k = 5, dist_threshold = NULL, style = "W") {
  # Check if input is an sf object
  if (!inherits(data, "sf")) {
    stop("data must be an sf object")
  }
  
  # Create neighbors based on type
  if (type == "contiguity") {
    # Queen contiguity
    nb <- spdep::poly2nb(as(data, "Spatial"), queen = TRUE)
  } else if (type == "distance") {
    # Distance-based
    coords <- sf::st_coordinates(sf::st_centroid(data))
    if (is.null(dist_threshold)) {
      stop("dist_threshold must be provided for distance-based weights")
    }
    nb <- spdep::dnearneigh(coords, 0, dist_threshold)
  } else if (type == "knn") {
    # K-nearest neighbors
    coords <- sf::st_coordinates(sf::st_centroid(data))
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k = k))
  } else {
    stop("Unsupported weights type: ", type)
  }
  
  # Create weights list
  listw <- spdep::nb2listw(nb, style = style, zero.policy = TRUE)
  
  return(listw)
}

#' Calculate global spatial autocorrelation
#'
#' @param data An sf object
#' @param variable Name of the variable to analyze
#' @param weights Spatial weights object (if NULL, queen contiguity weights are used)
#' @param method Method for autocorrelation ("moran", "geary")
#'
#' @return A list with test results
#' @export
#'
#' @examples
#' \dontrun{
#' # Test for spatial autocorrelation in population density
#' moran_result <- global_autocorrelation(counties, "pop_density")
#' }
global_autocorrelation <- function(data, variable, weights = NULL, method = "moran") {
  # Check if input is an sf object
  if (!inherits(data, "sf")) {
    stop("data must be an sf object")
  }
  
  # Check if variable exists
  if (!variable %in% names(data)) {
    stop("Variable not found in data: ", variable)
  }
  
  # Create weights if not provided
  if (is.null(weights)) {
    weights <- spatial_weights(data, "contiguity")
  }
  
  # Extract variable values
  x <- data[[variable]]
  
  # Calculate autocorrelation
  if (method == "moran") {
    result <- spdep::moran.test(x, weights, zero.policy = TRUE)
  } else if (method == "geary") {
    result <- spdep::geary.test(x, weights, zero.policy = TRUE)
  } else {
    stop("Unsupported method: ", method)
  }
  
  return(result)
}

#' Calculate local spatial autocorrelation
#'
#' @param data An sf object
#' @param variable Name of the variable to analyze
#' @param weights Spatial weights object (if NULL, queen contiguity weights are used)
#' @param method Method for local autocorrelation ("moran", "getis_ord")
#'
#' @return The input sf object with additional columns for local autocorrelation statistics
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate local Moran's I for population density
#' counties_lisa <- local_autocorrelation(counties, "pop_density")
#' }
local_autocorrelation <- function(data, variable, weights = NULL, method = "moran") {
  # Check if input is an sf object
  if (!inherits(data, "sf")) {
    stop("data must be an sf object")
  }
  
  # Check if variable exists
  if (!variable %in% names(data)) {
    stop("Variable not found in data: ", variable)
  }
  
  # Create weights if not provided
  if (is.null(weights)) {
    weights <- spatial_weights(data, "contiguity")
  }
  
  # Extract variable values
  x <- data[[variable]]
  
  # Calculate local autocorrelation
  result <- data
  
  if (method == "moran") {
    # Local Moran's I
    localmoran_result <- spdep::localmoran(x, weights, zero.policy = TRUE)
    result$local_moran <- localmoran_result[, "Ii"]
    result$local_moran_p <- localmoran_result[, "Pr(z != E(Ii))"]
    
    # Classify into LISA categories
    z_x <- scale(x)
    lag_z_x <- spdep::lag.listw(weights, z_x)
    
    result$lisa_category <- "Not Significant"
    result$lisa_category[result$local_moran_p <= 0.05 & z_x > 0 & lag_z_x > 0] <- "High-High"
    result$lisa_category[result$local_moran_p <= 0.05 & z_x < 0 & lag_z_x < 0] <- "Low-Low"
    result$lisa_category[result$local_moran_p <= 0.05 & z_x > 0 & lag_z_x < 0] <- "High-Low"
    result$lisa_category[result$local_moran_p <= 0.05 & z_x < 0 & lag_z_x > 0] <- "Low-High"
    
  } else if (method == "getis_ord") {
    # Getis-Ord G*
    localg_result <- spdep::localG(x, weights, zero.policy = TRUE)
    result$local_g <- as.numeric(localg_result)
    
    # Classify hotspots and coldspots
    result$hotspot_category <- "Not Significant"
    result$hotspot_category[result$local_g > 1.96] <- "Hotspot (95%)"
    result$hotspot_category[result$local_g > 2.58] <- "Hotspot (99%)"
    result$hotspot_category[result$local_g < -1.96] <- "Coldspot (95%)"
    result$hotspot_category[result$local_g < -2.58] <- "Coldspot (99%)"
    
  } else {
    stop("Unsupported method: ", method)
  }
  
  return(result)
}

#' Perform point pattern analysis
#'
#' @param points An sf object with point geometries
#' @param window An sf object defining the study area (optional)
#' @param method Analysis method ("density", "kfunction", "lfunction")
#' @param bw Bandwidth for kernel density estimation (if method = "density")
#' @param nsim Number of simulations for envelope calculation (if method = "kfunction" or "lfunction")
#'
#' @return A list with analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' # Calculate kernel density
#' density_result <- point_pattern_analysis(crime_points, city_boundary, "density")
#'
#' # Calculate K-function with simulation envelope
#' kfunc_result <- point_pattern_analysis(crime_points, city_boundary, "kfunction", nsim = 99)
#' }
point_pattern_analysis <- function(points, window = NULL, method = "density", bw = NULL, nsim = 39) {
  # Check if input is an sf object
  if (!inherits(points, "sf")) {
    stop("points must be an sf object")
  }
  
  if (!is.null(window) && !inherits(window, "sf")) {
    stop("window must be an sf object")
  }
  
  # Extract point coordinates
  coords <- sf::st_coordinates(points)
  
  # Create spatstat point pattern
  if (is.null(window)) {
    # Create default window from point extent
    xrange <- range(coords[, 1])
    yrange <- range(coords[, 2])
    win <- spatstat::owin(xrange = xrange, yrange = yrange)
  } else {
    # Convert sf polygon to spatstat window
    win <- as(as(window, "Spatial"), "owin")
  }
  
  # Create point pattern
  pp <- spatstat::ppp(coords[, 1], coords[, 2], window = win)
  
  # Perform analysis based on method
  if (method == "density") {
    # Kernel density estimation
    if (is.null(bw)) {
      bw <- spatstat::bw.diggle(pp)
    }
    density <- spatstat::density.ppp(pp, sigma = bw)
    return(list(density = density, bandwidth = bw))
    
  } else if (method == "kfunction") {
    # K-function
    k <- spatstat::Kest(pp)
    
    # Simulation envelope
    if (nsim > 0) {
      k_env <- spatstat::envelope(pp, fun = spatstat::Kest, nsim = nsim)
      return(list(k = k, envelope = k_env))
    }
    
    return(list(k = k))
    
  } else if (method == "lfunction") {
    # L-function
    l <- spatstat::Lest(pp)
    
    # Simulation envelope
    if (nsim > 0) {
      l_env <- spatstat::envelope(pp, fun = spatstat::Lest, nsim = nsim)
      return(list(l = l, envelope = l_env))
    }
    
    return(list(l = l))
    
  } else {
    stop("Unsupported method: ", method)
  }
}

#' Perform spatial interpolation
#'
#' @param points An sf object with point geometries and values to interpolate
#' @param variable Name of the variable to interpolate
#' @param grid An sf object or raster defining the prediction grid
#' @param method Interpolation method ("idw", "kriging")
#' @param params Additional parameters for the interpolation method
#'
#' @return A raster object with interpolated values
#' @export
#'
#' @examples
#' \dontrun{
#' # Interpolate temperature using IDW
#' temp_surface <- spatial_interpolation(weather_stations, "temperature", county_grid, "idw")
#' }
spatial_interpolation <- function(points, variable, grid, method = "idw", params = list()) {
  # Check if input is an sf object
  if (!inherits(points, "sf")) {
    stop("points must be an sf object")
  }
  
  # Check if variable exists
  if (!variable %in% names(points)) {
    stop("Variable not found in points: ", variable)
  }
  
  # Extract point coordinates and values
  coords <- sf::st_coordinates(points)
  values <- points[[variable]]
  
  # Create spatial points data frame
  sp_points <- sp::SpatialPointsDataFrame(
    coords = coords,
    data = data.frame(z = values),
    proj4string = sp::CRS(sf::st_crs(points)$proj4string)
  )
  
  # Create prediction grid
  if (inherits(grid, "sf")) {
    # Convert sf to raster grid
    bbox <- sf::st_bbox(grid)
    grid_raster <- raster::raster(
      xmn = bbox["xmin"], xmx = bbox["xmax"],
      ymn = bbox["ymin"], ymx = bbox["ymax"],
      crs = sf::st_crs(grid)$proj4string,
      resolution = params$resolution %||% c(100, 100)
    )
  } else if (inherits(grid, "Raster")) {
    grid_raster <- grid
  } else {
    stop("grid must be an sf object or Raster object")
  }
  
  # Perform interpolation
  if (method == "idw") {
    # Inverse distance weighting
    power <- params$power %||% 2
    idw <- gstat::idw(z ~ 1, sp_points, newdata = raster::rasterToPoints(grid_raster, spatial = TRUE), idp = power)
    result <- raster::raster(idw)
    
  } else if (method == "kriging") {
    # Ordinary kriging
    model <- params$model %||% "Sph"
    vgm_params <- gstat::fit.variogram(
      gstat::variogram(z ~ 1, sp_points),
      gstat::vgm(model = model)
    )
    kriging <- gstat::krige(z ~ 1, sp_points, newdata = raster::rasterToPoints(grid_raster, spatial = TRUE), model = vgm_params)
    result <- raster::raster(kriging)
    
  } else {
    stop("Unsupported method: ", method)
  }
  
  return(result)
}

#' Calculate spatial join
#'
#' @param x An sf object
#' @param y An sf object
#' @param join_type Type of spatial join ("intersects", "contains", "within", "touches")
#' @param suffix Suffixes for duplicate column names
#'
#' @return An sf object with joined attributes
#' @export
#'
#' @examples
#' \dontrun{
#' # Join points with polygons
#' points_with_attributes <- spatial_join(points, polygons)
#' }
spatial_join <- function(x, y, join_type = "intersects", suffix = c(".x", ".y")) {
  # Check if inputs are sf objects
  if (!inherits(x, "sf") || !inherits(y, "sf")) {
    stop("x and y must be sf objects")
  }
  
  # Ensure both datasets have the same CRS
  if (sf::st_crs(x) != sf::st_crs(y)) {
    y <- sf::st_transform(y, sf::st_crs(x))
  }
  
  # Perform spatial join based on join type
  if (join_type == "intersects") {
    join <- sf::st_join(x, y, join = sf::st_intersects, suffix = suffix)
  } else if (join_type == "contains") {
    join <- sf::st_join(x, y, join = sf::st_contains, suffix = suffix)
  } else if (join_type == "within") {
    join <- sf::st_join(x, y, join = sf::st_within, suffix = suffix)
  } else if (join_type == "touches") {
    join <- sf::st_join(x, y, join = sf::st_touches, suffix = suffix)
  } else {
    stop("Unsupported join type: ", join_type)
  }
  
  return(join)
}

