#' This function is used to simulated distances from a detection with 
#' adjustments. Covariates are not permitted. The detection function consists of 
#' a key function and adjustment series. The permitted key functions are: 
#' uniform ('unif'), half-normal ('hn') and hazard-rate ('hr'). The permitted 
#' adjustment series are cosine ('cos'), simple polynomial ('poly') and 
#' Hermite polynomial ('herm'). The maximum number of adjustment terms is 5, 
#' which means that the maximum number of parameters in a detection function is 
#' 7, which is when the key is hazard-rate. 
#' 
#' The permitted orders of adjustment terms depend on the key-adjustment 
#' combination:
#' Key          Adjustment          orders
#' unif         cos                 1, 2, 3, 4, 5
#' unif         poly                2, 4, 6, 8, 10
#' unif         herm                2, 4, 6, 8, 10
#' hn           cos                 2, 3, 4, 5, 6
#' hn           poly                4, 6, 8, 10, 12
#' hn           herm                4, 6, 8, 10, 12
#' hr           cos                 2, 3, 4, 5, 6
#' hr           poly                4, 6, 8, 10, 12
#' hr           herm                4, 6, 8, 10, 12
#' 
#' @param n the number of distances to be simulated
#' @param key character indicating the key function
#' @param adj character indicating the adjustment series
#' @param orders vector of orders. Minimum length is 1, max length is 5
#' @param width the width or right truncation distance
#' @param left the left truncation distance (default is 0)
#' @param point boolean for point or line transect data (default is TRUE)
#' @param pts the number of distances between min
#' 
#' 
#' 
#' 
simDists <- function(n, key, adj, orders, parameters, 
                     width, scaling = "width", left = 0, 
                     mono = TRUE, mono_strict = TRUE, point = TRUE, 
                     pts = 10000, plot_g = TRUE) {
  # n <- 100; key <- "hn"; adj <- "cos"; orders <- c(2, 3, 4); width <- 100
  # scaling <- "width"; left <- 0; mono <- TRUE; mono_strict <- TRUE
  # point <- TRUE; pts <- 10000; parameters <- c(1.5, 0.5, 1, 0.4)
  
  ## Unit testing
  if (n <= 0 | is.infinite(n)) {
    stop("n has to be positive and real.")
  }
  if (!key %in% c("unif", "hn", "hr")){
    stop("key has to be either 'unif', 'hn' or 'hr'.")
  }
  if (!adj %in% c("cos", "poly", "herm")){
    stop("key has to be either 'cos', 'poly' or 'herm'.")
  }
  if (length(orders) < 1 | length(orders) > 5) {
    stop("A minimum of 1 and maximum of 5 adjustment terms are permitted.")
  }

  ## Check validity of detection function 
  if (key == "unif" & adj == "cos" & !all(orders %in% c(1, 2, 3, 4, 5))) {
    stop("The unif+cos combination only allows orders 1, 2, 3, 4 and 5.")
  }
  if (key == "hn" & adj == "cos" & !all(orders %in% c(2, 3, 4, 5, 6))) {
    stop("The hn+cos combination only allows orders 2, 3, 4, 5 and 6.")
  }
  if (key == "hr" & adj == "cos" & !all(orders %in% c(2, 3, 4, 5, 6))) {
    stop("The hr+cos combination only allows orders 2, 3, 4, 5 and 6.")
  }
  
  if (key == "unif" & adj %in% c("poly", "herm") & 
      !all(orders %in% c(2, 4, 6, 8, 10))) {
    stop("The unif+poly and unif+herm combinations only allow orders 2, 4, 6, 8 and 10.")
  }
  
  if (key != "unif" & adj %in% c("poly", "herm") & 
      !all(orders %in% c(4, 6, 8, 10, 12))) {
    stop("The hn+poly, hn+herm, hr+poly and hr+herm combinations only allow orders 4, 6, 8, 10 and 12.")
  }
  
  ## Check the scaling
  if (key == "unif" | scaling == "scale") {
    stop("Not possible to scale distances by the scale parameter with uniform key function.")
  }
  
  ## Extract the parameters
  if (key == "unif") {
    pars <- list(adjustment = parameters)
  } else if (key == "hn") {
    pars <- list(scale = parameters[1], 
                adjustment = parameters[-1])
  } else {
    pars <- list(shape = parameters[1],
                 scale = parameters[2], 
                 adjustment = parameters[-c(1,2)])
  }
  
  ## Create a fake ddfobj to pass onto distpdf()
  distance <- seq(from = left, to = width, length = pts)

  ddfobj <- list()
  ddfobj$type <- key
  if (key %in% c("hn", "hr")) {
    ddfobj$scale$dm <- rep(1, pts)
    ddfobj$scale$parameters <- pars$scale
    }
  if (key == "hr") {
    ddfobj$shape$dm <- rep(1, pts)
    ddfobj$shape$parameters <- pars$shape
  }
  
  ddfobj$adjustment$series <- adj
  ddfobj$adjustment$scale <- scaling
  ddfobj$adjustment$order <- orders
  ddfobj$adjustment$parameters <- pars$adjustment
  ddfobj$adjustment$exp <- FALSE
  
  ## Derive non-normalised pdfs for the distances
  dist_pdf <- distpdf(distance = distance, ddfobj = ddfobj, point = point, 
                      left = left, width = width)
  
  if (any(dist_pdf < 0)) {
    cat("The pdf is negative for some distances between left and width! These were capped at 0.\n")
    dist_pdf[dist_pdf < 0] <- 0
  }
  
  sampled_dists <- sample(distance, size = n, replace = TRUE, prob = dist_pdf)
  
  if (plot_g) {
    adj_pars <- rep(NA, 5)
    
    if (key == "unif" & adj == "cos") {
      indices <- c(1:5) %in% orders
    }
    if (key != "unif" & adj == "cos") {
      indices <- c(2:6) %in% orders
    }
    if (key == "unif" & adj %in% c("poly", "herm")) {
      indices <- c(2, 4, 6, 8, 10) %in% orders
    }
    if (key != "unif" & adj %in% c("poly", "herm")) {
      indices <- c(4, 6, 8, 10, 12) %in% orders
    }
    
    adj_pars[indices] <- pars$adjustment
    
    df <- data.frame(
        key = key,
        adj = adj,
        scale = ifelse(is.null(pars$scale), NA, pars$scale),
        shape = ifelse(is.null(pars$shape), NA, pars$shape),
        adj1 = adj_pars[1],
        adj2 = adj_pars[2],
        adj3 = adj_pars[3],
        adj4 = adj_pars[4],
        adj5 = adj_pars[5],
        left = left,
        width = width,
        scaling = scaling
    )
    
    plotDetfct(df, print_request = FALSE)
  }
  
  return(sampled_dists)
}
