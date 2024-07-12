#' This function is used plot a bunch of detection functions that are key+adj.
#' Covariates are not permitted. The detection function consists of 
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
#' @param specs a data frame with the columns:
#'   'key' is the key function, 'adj' is the adjustment series, 'scale' is the 
#'   scale paramter (can be NA), 'shape' is the shape parameter (can be NA),
#'   'adj1--adj5' are the adjustment coefficients (can be NA but not all) and 
#'   the orders will be determined based on which coefficients are provided and
#'   the key+adjustment combination, 'left' is the left truncation and 'width'
#'    is the right truncation distance, and scaling. 
#' @param pts the number of distances between min
#' @param print_request whether it should request pressing 'enter' between the 
#' plotting. Defaults to TRUE and should generally not be touched. 
#' 
plotDetfct <- function(specs, pts = 10000, print_request = TRUE) {

  # specs <- data.frame(
  #   key = "unif",
  #   adj = "herm",
  #   scale = NA,
  #   shape = NA,
  #   adj1 = 0.5,
  #   adj2 = 1,
  #   adj3 = 0.4,
  #   adj4 = NA,
  #   adj5 = NA,
  #   left = 0,
  #   width = 100,
  #   scaling = 'width'
  # ); i <- 1; pts <- 10000; print_request <- TRUE

  
  for (i in 1:nrow(specs)) {
    
    if (print_request) {
      readline(prompt = paste0("Press [enter] to print plot ", i, ":"))
    }
    
    key <- specs[i, 1]
    adj <- specs[i, 2]
    left <- specs[i, "left"]
    width <- specs[i, "width"]
    scaling <- specs[i, "scaling"]
    ## Check validity of detection function 
    if (key == "unif" & adj == "cos") {
      orders <- c(1, 2, 3, 4, 5)[!is.na(specs[i, c(5:9)])]
    } else if (key != "unif" & adj == "cos") {
      orders <- c(2, 3, 4, 5, 6)[!is.na(specs[i, c(5:9)])]
    }else if (key == "unif" & adj %in% c("poly", "herm")) {
      orders <- c(2, 4, 6, 8, 10)[!is.na(specs[i, c(5:9)])]
    } else if (key != "unif" & adj %in% c("poly", "herm")) {
      orders <- c(4, 6, 8, 10, 12)[!is.na(specs[i, c(5:9)])]
    } 
    
    parameters <- unlist(specs[i, c(3:9)][!is.na(specs[i, c(3:9)])])

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
    g <- detfct(distance = distance, ddfobj = ddfobj,
                        left = left, width = width, standardize = TRUE)
    
    # col <- rep("black", pts)
    # col[g < 0] <- "red"
    plot(y = g, x = distance, col = "black", type = "l")
    abline(h = 0, col = "red", lty = 2)
    
    if (any(g < 0)) {
      cat("The pdf is negative for some distances between left and width.\n")
    }
  }
}
