#' Add covariate levels detection function plots
#'
#' Add a line or lines to a plot of the detection function which correspond to a a given covariate combination. These can be particularly useful when there is a small number of factor levels or if quantiles of a continuous covariate are specified.
#'
#' All covariates must be specified in \code{data}. Plots can become quite busy when this approach is used. It may be useful to fix some covariates at their median level and plot set values of a covariate of interest. For example setting weather (e.g., Beaufort) to its median and plotting levels of observer, then creating a second plot for a fixed observer with levels of weather.
#'
#' Arguments to \code{\link{lines}} are supplied in \dots and aesthetics like line type (\code{lty}), line width (\code{lwd}) and colour (\code{col}) are recycled. By default \code{lty} is used to distinguish between the lines. It may be useful to add a \code{\link{legend}} to the plot (lines are plotted in the order of \code{data}).
#'
# Update Distance::add_df_covar_line when updating parameters here
#' @param ddf a fitted detection function object
#' @param data a \code{data.frame} with the covariate combination you want to plot
#' @param \dots extra arguments to give to \code{\link{line}} (\code{lty}, \code{lwd}, \code{col})
#' @param ndist number of points to evaluate the detection function at
#' @return invisibly, the values of detectability over the truncation range
#'
#' @export
#' @author David L Miller
#'
#' @examples
#' \dontrun{
#' # fit an example model
#' data(book.tee.data)
#' egdata <- book.tee.data$book.tee.dataframe
#' result <- ddf(dsmodel = ~mcds(key = "hn", formula = ~sex),
#'               data = egdata[egdata$observer==1, ], method = "ds",
#'               meta.data = list(width = 4))
#'
#' # make a base plot, showpoints=FALSE makes the plot less busy
#' plot(result, showpoints=FALSE)
#'
#' # add lines for sex one at a time
#' add_df_covar_line(result, data.frame(sex=0), lty=2)
#' add_df_covar_line(result, data.frame(sex=1), lty=3)
#'
#' # add a legend
#' legend(3, 1, c("Average", "sex==0", "sex==1"), lty=1:3)
#'
#' # alternatively we can add both at once
#' # fixing line type and varying colour
#' plot(result, showpoints=FALSE)
#' add_df_covar_line(result, data.frame(sex=c(0,1)), lty=1,
#'                   col=c("red", "green"))
#' # add a legend
#' legend(3, 1, c("Average", "sex==0", "sex==1"), lty=1,
#'        col=c("black", "red", "green"))
#' }
add_df_covar_line <- function(ddf, data, ndist=250, ...){

  # if we have a Distance object rather than mrds, use that
  if(all(class(ddf)=="dsmodel")){
    df <- ddf$ddf
  }else{
    df <- ddf
  }

  left <- df$meta.data$left
  width <- df$meta.data$width

  # distance range to plot over
  xx <- seq(left, width, length.out=ndist)

  #data <- model.frame(df$ds$aux$ddfobj$scale$formula, data)
  xm <- df$ds$aux$ddfobj$xmat

  data$object <- 1:nrow(data)
  data$distance <- rep(0, nrow(data))
  data$observer <- rep(0, nrow(data))
  data$detected <- rep(0, nrow(data))
  data$binned <- rep(df$ds$aux$ddfobj$xmat$binned[1], nrow(data))

  eval_with_covars <- function(distance, newdata, model){
    ddfobj <- model$ds$aux$ddfobj

    fpar <- model$par
    ddfobj <- assign.par(ddfobj, fpar)
    # Get integration ranges either from specified argument or from
    # values stored in the model.
    if(is.data.frame(newdata)){
      nr <- nrow(newdata)
    }else{
      nr <- 1
    }

    if(is.null(model$ds$aux$int.range)){
      int.range <- cbind(rep(0, nr), rep(width, nr))
    }else{
      int.range <- model$ds$aux$int.range
      if(is.vector(int.range)){
        int.range <- cbind(rep(int.range[1], nr),
                           rep(int.range[2], nr))
      #}else if(nrow(int.range) == (nrow(x)+1)){
      #int.range <- int.range[2:nrow(int.range), , drop=FALSE]
      }
    }

    # set the distance column to be the left truncation distance
    # this gets around an issue that Nat Kelly found where later process.data
    # will remove entires with distance < left truncation
    # BUT save the NAs!
    nas <- is.na(newdata$distance)
    newdata$distance <- left
    newdata$distance[nas] <- NA

    newdata_save <- newdata

    # get the data in the model
    model_dat <- model$data

    # counter for NAs...
    naind <- rep(FALSE, nrow(newdata))

    # do this for both scale and shape parameters
    for(df_par in c("scale", "shape")){
      # if that parameter exists...
      if(!is.null(ddfobj[[df_par]])){
        # save the column names from the design matrix
        znames <- colnames(ddfobj[[df_par]]$dm)

        # pull out the columns in the formula and the distances column
        fvars <- all.vars(as.formula(model$ds$aux$ddfobj[[df_par]]$formula))

        if(!all(fvars %in% colnames(newdata))){
          stop("columns in `newdata` do not match those in fitted model\n")
        }

        model_dat <- model_dat[, c("distance", fvars), drop=FALSE]

        if(df_par=="scale"){
          # which rows have NAs?
          naind <- naind | apply(newdata_save[, c("distance", fvars), drop=FALSE],
                                 1, function(x) any(is.na(x)))
        }

        # setup the covariate matrix, using the model data to ensure that
        # the levels are right
        newdata <- rbind(model_dat,
                         newdata_save[, c("distance", fvars), drop=FALSE])
        dm <- setcov(newdata, as.formula(ddfobj[[df_par]]$formula))

        # now check that the column names are the same for the model
        # and prediction data matrices
        if(!identical(colnames(dm), znames)){
          stop("fields or factor levels in `newdata` do not match data used in fitted model\n")
        }

        # get only the new rows for prediction
        dm <- dm[(nrow(model_dat)+1):nrow(dm), , drop=FALSE]
        # assign that!
        ddfobj[[df_par]]$dm <- dm

      }
    }

    # handle data setup for uniform key case
    if(ddfobj$type == "unif"){
      model_dat <- model_dat[, "distance", drop=FALSE]
      # which rows have NAs?
      naind <- is.na(newdata_save$distance)

      newdata <- rbind(model_dat,
                       newdata_save[, "distance", drop=FALSE])
      dm <- setcov(newdata, ~1)
      dm <- dm[(nrow(model_dat)+1):nrow(dm), , drop=FALSE]
    }

    # get the bins when you have binned data
    # use the breaks specified in the model!
    if(model$meta.data$binned){
      nanana <- apply(newdata[, c("distance", fvars), drop=FALSE],
                      1, function(x) any(is.na(x)))
      newdata_b <- create.bins(newdata[!nanana, , drop=FALSE], model$meta.data$breaks)
      newdata$distbegin <- NA
      newdata$distend <- NA
      newdata[!nanana, ] <- newdata_b
    }

    # update xmat too
    datalist <- process.data(newdata, model$meta.data, check=FALSE)
    ddfobj$xmat <- datalist$xmat[(nrow(model_dat)+1):nrow(datalist$xmat),,drop=FALSE]
    ddfobj$xmat <- ddfobj$xmat[!naind, , drop=FALSE]
    int.range <- int.range[!naind, , drop=FALSE]
    # reset newdata to be the right thing
    newdata <- newdata[(nrow(model_dat)+1):nrow(newdata), , drop=FALSE]


    detfct(distance, ddfobj, select=NULL, index=NULL, width=width,
           standardize=TRUE, stdint=FALSE, left=left)
  }

  # fiddle to get nice lty behaviour by default
  lines_args <- list(...)
  if(is.null(lines_args$lty)){
    lty <- rep_len(2:6, nrow(data))
  }else{
    lty <- rep_len(lines_args$lty, nrow(data))
  }

  # recycle width if necessary
  if(is.null(lines_args$lwd)){
    lwd <- rep_len(1, nrow(data))
  }else{
    lwd <- rep_len(lines_args$lwd, nrow(data))
  }
  # and colour
  if(is.null(lines_args$col)){
    col <- rep_len("black", nrow(data))
  }else{
    col <- rep_len(lines_args$col, nrow(data))
  }


  # storage
  linedat <- matrix(NA, nrow(data), length(xx))

  # now loop over the data rows
  for(i in 1:nrow(data)){
    # evaluate and save data
    linedat[i,] <- eval_with_covars(xx, data[i, ], df)
    # plot
    lines_args$lty <- lty[i]
    lines_args$lwd <- lwd[i]
    lines_args$col <- col[i]
    lines_args$x <- xx
    lines_args$y <- linedat[i, ]
    do.call(lines, lines_args)
  }

  # return saved data
  invisible(linedat)
}
