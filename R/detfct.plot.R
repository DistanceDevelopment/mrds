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
detfct.plot <- function(specs, pts = 10000, print_request = TRUE, 
                        print_title = TRUE) {

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
    
    key <- specs[i, "key"]
    adj <- specs[i, "adj"]
    left <- specs[i, "left"]
    width <- specs[i, "width"]
    scaling <- specs[i, "scaling"]
    
    ## Extract the correct orders
    if (key == "unif" & adj == "cos") {
      orders <- c(1, 2, 3, 4, 5)[!is.na(specs[i, paste0("adj", 1:5)])]
    } else if (key != "unif" & adj == "cos") {
      orders <- c(2, 3, 4, 5, 6)[!is.na(specs[i, paste0("adj", 1:5)])]
    }else if (key == "unif" & adj %in% c("poly", "herm")) {
      orders <- c(2, 4, 6, 8, 10)[!is.na(specs[i, paste0("adj", 1:5)])]
    } else if (key != "unif" & adj %in% c("poly", "herm")) {
      orders <- c(4, 6, 8, 10, 12)[!is.na(specs[i, paste0("adj", 1:5)])]
    } 
    
    parameters <- unlist(specs[i, c("shape", "scale", paste0("adj", 1:5))])
    parameters <- parameters[!is.na(specs[i, c("shape", "scale", paste0("adj", 1:5))])]
    
    ## Extract the parameters
    if (key == "unif") {
      pars <- list(adjustment = parameters)
    } else if (key == "hn") {
      pars <- list(scale = parameters["scale"], 
                   adjustment = parameters[names(parameters) != "scale"])
    } else {
      pars <- list(shape = parameters["shape"],
                   scale = parameters["scale"], 
                   adjustment = parameters[!names(parameters) %in% c("shape", 
                                                                     "scale")])
    }
    
    ## Create a fake ddfobj to pass onto distpdf()
    distance <- seq(from = left, to = width, length = pts)
    
    ddfobj <- list()
    ddfobj$type <- key
    if (key %in% c("hn", "hr")) {
      ddfobj$scale$dm <- matrix(rep(1, pts), ncol = 1)
      ddfobj$scale$parameters <- pars$scale
    }
    if (key == "hr") {
      ddfobj$shape$dm <- matrix(rep(1, pts), ncol = 1)
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

    # average_p <- integrate(detfct, left, width, ddfobj = ddfobj, left =  left,
    #                        index = 1, width = width,
    #                        standardize = T)$value / (width - left)
    
    line_p <- integratepdf(ddfobj, select = c(T, rep(F, pts-1)), width = width,
                           left = left, 
                           standardize = T, point = F, 
                           int.range = c(left, width))
    point_p <- integratepdf(ddfobj, select = c(T, rep(F, pts-1)), width = width,
                            left = left,
                            standardize = T, point = T, 
                            int.range = c(left, width))
    
    # col <- rep("black", pts)
    # col[g < 0] <- "red"
    if (print_title) {
      plot(y = g, x = distance, col = "black", type = "l", 
           main = paste0("Df with key: ", key, 
                         "; parameters: ", paste(parameters, collapse = ", "), 
                         ";\nadjustment orders: ", paste(orders, collapse = ", "), 
                         ". Distances are scaled by ", scaling, "."))
      text(y = 1 * (max(g) - min(g)) + min(g), x = 0.8 * width, 
                    paste0("Line p: ", round(line_p, 3)))
      text(y = 0.9 * (max(g) - min(g)) + min(g), x = 0.8 * width, 
           paste0("Point p: ", round(point_p, 3)))
    } else {
      plot(y = g, x = distance, col = "black", type = "l")
      text(y = 1 * (max(g) - min(g)) + min(g), x = 0.8 * width, 
                    paste0("Line p: ", round(line_p, 3)))
      text(y = 0.9 * (max(g) - min(g)) + min(g), x = 0.8 * width, 
                      paste0("Point p: ", round(point_p, 3)))
    }

    abline(h = 0, col = "red", lty = 2)
    
    if (any(g < 0)) {
      cat("The pdf is negative for some distances between left and width.\n")
    }
  }
}
