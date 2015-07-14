#' Plot covariate levels for a detection function
#'
#' Utility routine to plot the levels of a covariate as lines on a plot. This only works for scale function covariates at the moment.
#'
#' @author David L Miller
ds_cov_levels <- function(x, distances){

  plot_data <- x$data
  plot_data <- get_all_vars(x$ds$aux$ddfobj$scale$formula, plot_data)

  # which are factors?
  plot_factors <- plot_data[,unlist(lapply(plot_data,is.factor)), drop=FALSE]
  # which aren't
  plot_nonfactors <- plot_data[,!unlist(lapply(plot_data,is.factor)),
                               drop=FALSE]

  # for the factors get the levels & for non-factors get quantiles
  plot_factors_q <- lapply(plot_factors, function(x)
                              factor(levels(x), levels=levels(x)))
  plot_nonfactors_q <- lapply(plot_nonfactors,
                              function(x) unique(quantile(x,
                                                          probs=c(0.25,0.5,0.75),
                                                          na.rm=TRUE)))
  # put that all in one list
  plot_all <- c(plot_factors_q, plot_nonfactors_q)

  # find the "fixed" values (when we vary other covariates)
  fixed_factors <- lapply(plot_factors, function(x)
                              factor(names(which.max(table(x))),
                                           levels=levels(x)))
  fixed_nonfactors <- lapply(plot_nonfactors, median, na.rm=TRUE)
  # put that all in one list
  fixed_all <- as.data.frame(c(fixed_factors, fixed_nonfactors))

  ## now build the data to plot with
  # initialise storage
  pd <- c()
  # loop over the covariates
  for(cov_i in seq_along(plot_all)){
    cov_name <- names(plot_all)[cov_i]
    # loop over levels of the covariates
    for(val_i in seq_along(plot_all[[cov_i]])){
      # take the covariate, make a new row with fixed values
      # then replace the covariate we want with it's plot value
      new_plot_data <- fixed_all
      new_plot_data[[cov_name]] <- plot_all[[cov_i]][val_i]
      # fold that into the data
      pd <- rbind.data.frame(pd,new_plot_data)
    }
  }
  pd <- unique(pd)

  # how many rows in the data?
  npd <- nrow(pd)

  # replicate the data
  pd <- pd[rep(1:nrow(pd),rep(length(distances),nrow(pd))),,drop=FALSE]

  # fill in the intercept
  pd[["(Intercept)"]] <- 1

  # replicate the distances
  distances <- rep(distances,npd)

  # calculate the detection function evaluations
  ddfobj2 <- x$ds$aux$ddfobj
#  ddfobj2$scale$dm <- pd
  ddfobj2$scale$dm <- model.matrix(as.formula(ddfobj2$scale$formula), pd)
  df_line <- detfct(distances, ddfobj2, width=x$meta.data$width)

  # split up into 1 list element per line
  df_line <- split(df_line, rep(1:npd,rep(length(distances)/npd,npd)))

  return(list(df_line=df_line, plot_all=plot_all))
}
