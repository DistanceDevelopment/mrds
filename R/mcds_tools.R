#' Run MCDS.exe as a backend for mrds
#'
#' Rather than use the R-based detection function fitting algorithms provided in 
#' `mrds`, one can also use the algorithm used by Distance for Windows, implemented
#' in the binary file `MCDS.exe`.  Note that with changes in R-based optimizer introduced 
#' in `mrds` version 3.0.0 this is unlikely to result in better estimates.  
#' The option remains available, although it may be deprecated in a future release.

#' To make use of this facility, one must first download the `MCDS.exe` binary, as
#' laid out below under `Obtaining MCDS.exe`.  Once the binary is installed, calls
#' to `ddf` will, by default, result in using the model being fit using both `MCDS.exe`
#' and the R-based algorithm, and the one with lower negative log-likelihood being selected.
#' In almost all cases, both algorithms produce the same results, but we have found edge
#' where one or other fails to find the likelihood maximum and hence trying both is useful.
#' 
#' There may also be cases where the `MCDS.exe` algorithm is faster than the R-based one.
#' Under this circumstance, you can choose to run only the `MCDS.exe` algorithm via by setting
#' the `ddf` argument \code{control=list(optimizer='MCDS')}.  For completeness, one can also
#' choose to use only the R-based algorithm by setting \code{control=list(optimizer='R')}.
#' 
#' For more information and examples comparing the R-based and `MCDS.exe` algorithms,
#' see our examples pages at https://examples.distancesampling.org/
#'
#' If you are running a non-Windows operating system, you can follow the
#' instructions below to have `MCDS.exe` run using `wine`.
#'
#' @aliases MCDS mcds_dot_exe
#' @name MCDS.exe
#' @rdname mcds_dot_exe
#'
#' @section Obtaining MCDS.exe:
#' The following code can be used to download `MCDS.exe` from the distance 
#' sampling website:
#' \code{download.file("http://distancesampling.org/R/MCDS.exe", paste0(system.file(package="mrds"),"/MCDS.exe"), mode = "wb")}
#' The MCDS binary will be installed to the main directory of your your local R 
#' mrds library. Alternatively, you can copy the `MCDS.exe` from your local 
#' Distance for Windows installation if you prefer. The location of your local 
#' mrds library main directory can be found by running the following in R:
#' \code{system.file("MCDS.exe", package="mrds")}.
#'
#' @section Running MCDS.exe on non-Windows platforms:
#' This has been tentatively tested on a mac but should currently be 
#' considered largely experimental.
#' 
#' One can still use MCDS.exe even if you are running a mac computer. To
#' do this one will need to install `wine` a Windows emulator. It is important
#' to use a version of `wine` which can run 32-bit programs.
#'
#' The package will attempt to work out which `wine` binary to use (and detect
#' if it is installed), but this doesn't always work. In this case, the
#' location of the `wine` binary can be specified in the `control` `list`
#' provided to `ddf` using the `winebin` element or supply the `winebin` 
#' argument to the `ds` function. For example, if `wine` is installed at 
#' `/usr/bin/local/wine` you can set `control$winebin` to that
#' location to use that binary.
#'
#' On macOS, this can be achieved using the `homebrew` package management
#' system and installing the `wine-crossover` package. You may need to change
#' the \code{control$winebin} to be `wine`, `wine64` or `wine32on64`, 
#' depending on your system's setup. This package tries to work out what to 
#' do, but likely doesn't handle all corner cases. Currently this is untested 
#' on Mac M1 systems.
#'
#' @section Stopping using MCDS.exe:
#' Once this feature is enabled, using `ddf` will by default run both
#' its built-in R optimizer and `MCDS.exe`. To disable this behaviour, specify which 
#' you wish to use with via the \code{optimizer=} option described above.  Alternatively,
#' if you wish to permanently stop using MCDS.exe, remove
#' the `MCDS.exe` binary file. You can find which folder it is in by running the following in R:
#' \code{system.file("MCDS.exe", package="mrds")}.
#'
#' @author David L Miller and Jonah McArthur
NULL


#' @title create.command.file
# a function to take the input as required for the ddf function in mrds and
# create a command file for the mcds engine, which will perform an equivalent
# analysis.
#' @inheritParams ddf
#' @rdname create.command.file
#' @author Jonah McArthur
#' @importFrom utils write.table
create.command.file <- function(dsmodel=call(), data,
                                method, meta.data, control) {
  # create a temporary directory
  directory <- tempdir()
  
  # create command file
  command.file.name <- tempfile(pattern="cmdtmp", 
                                tmpdir=directory,
                                fileext=".txt")
  file.create(command.file.name)
  
  
  
  out.file <- tempfile(pattern="out", 
                       tmpdir=directory,
                       fileext=".txt")
  
  log.file <- tempfile(pattern="log", 
                       tmpdir=directory,
                       fileext=".txt")
  
  stat.file <- tempfile(pattern="stat", 
                        tmpdir=directory,
                        fileext=".txt")
  
  plot.file <- tempfile(pattern="plot", 
                        tmpdir=directory,
                        fileext=".txt")
  
  data.file <- tempfile(pattern="data",
                        tmpdir=directory,
                        fileext=".txt")
  
  # Create Data File
  data.info <- create.data.file(data, dsmodel, data.file)
  # extract reformatted data
  data <- data.info$data
  # extract information about covariates 
  covar.fields <- data.info$covar.fields
  cluster <- data.info$cluster
  
  # HEADER section
  # ~~~~~~~~~~~~~~
  # specify the location of output files
  cat(out.file, file=command.file.name, "\n", append=TRUE)
  cat(log.file, file=command.file.name, "\n", append=TRUE)
  cat(stat.file, file=command.file.name, "\n", append=TRUE)
  cat(plot.file, file=command.file.name, "\n", append=TRUE)
  cat("None", file=command.file.name, "\n", append=TRUE)
  cat("None", file=command.file.name, "\n", append=TRUE)
  
  # OPTION section
  # ~~~~~~~~~~~~~~
  cat("OPTIONS;", file=command.file.name, "\n", append=TRUE)
  # specify the type of transect and its width and units
  if(meta.data$point == TRUE){
    cat(paste("DISTANCE=RADIAL /UNITS='Meters' /WIDTH=",
              meta.data$width, ";", sep=""), file=command.file.name, "\n",
        append=TRUE)
    cat("TYPE=POINT;", file=command.file.name, "\n", append=TRUE)
  }else{
    cat(paste("DISTANCE=PERP /UNITS='Meters' /WIDTH=",
              meta.data$width, ";", sep=""), file=command.file.name, "\n",
        append=TRUE)
    cat("TYPE=LINE;", file=command.file.name, "\n", append=TRUE)
  }
  
  # define whether there are clusters
  if(cluster == TRUE){
    cat("OBJECT=CLUSTER;", file=command.file.name, "\n", append=TRUE)
  }else{
    cat("OBJECT=SINGLE;", file=command.file.name, "\n", append=TRUE)
  }
  
  # managing the output options
  if(is.null(control$debug) == FALSE){
    if(control$debug == TRUE){
      cat("DEBUG=ON;", file=command.file.name, "\n", append=TRUE)
    }
  }
  if(is.null(control$showit) == FALSE){
    output_info_levels <- c("SUMMARY","RESULTS","SELECTION","ALL")
    specified_output_level <- output_info_levels[control$showit+1]
    cat(paste("PRINT=", specified_output_level, ";", sep=""),
        file=command.file.name, "\n", append=TRUE)
  }
  
  # the user will specify the adjustment term selection
  cat("SELECTION=SPECIFY;", file=command.file.name, "\n",
      append=TRUE)
  cat("END;", file=command.file.name, "\n", append=TRUE)
  
  # DATA section
  # ~~~~~~~~~~~~~~
  cat("DATA /STRUCTURE=FLAT;", file=command.file.name, "\n",
      append=TRUE)
  
  # combine data field names into single character value
  fields.comb <- paste(names(data), collapse = ",")
  cat(paste("FIELDS=", fields.comb, ";", sep=""),
      file=command.file.name, "\n", append=TRUE)
  
  # Don't tell MCDS about clusters
  # If clusters and cluster is a cov in the model
  # if(cluster && "SIZE" %in% covar.fields){
  #   cat(paste("SizeC=SIZE;", sep=""),
  #       file=command.file.name, "\n", append=TRUE)
  # }
  
  # Access the model parameters
  mod_paste <- paste(dsmodel)
  mod.vals <- try(eval(parse(text=mod_paste[2:length(mod_paste)])))
  
  # Get a list of the factor variables in the model formula
  factor.vars <- rownames(attr(terms(mod.vals$formula),"factors"))
  # Process to retain only the variable names
  for(i in seq(along = factor.vars)){
    # Find first bracket
    start <- gregexpr("\\(", factor.vars[i])[[1]][1] + 1
    # Find second bracket
    end <- gregexpr("\\)", factor.vars[i])[[1]][1] -1
    # Retain part within brackets
    factor.vars[i] <- substr(factor.vars[i], start, end)
  }
  
  # Add factor information
  factor.fields <- NULL
  for(i in seq(along = covar.fields)){
    # If either the data format is factor OR it is specified in formula as factor OR it is a character
    if(is.factor(data[[covar.fields[i]]]) || covar.fields[i] %in% factor.vars || is.character(data[[covar.fields[i]]])){
      factor.fields <- c(factor.fields, covar.fields[i])
      # Column in data may not be a factor field may only have been specified in formula
      labels <- levels(as.factor(data[[covar.fields[i]]]))
      # Reorder the labels - in R the baseline comes first but in MCDS the baseline comes last
      if(length(labels) > 1){
        labels <- c(labels[2:length(labels)], labels[1])
      }
      # Write to command file
      cat(paste("FACTOR /NAME=", toupper(covar.fields[i]),
                " /LEVELS=", length(labels), " /LABELS=",
                paste(labels, collapse = ", "), ";", sep=""), file=command.file.name, "\n",
          append=TRUE)
      
    }
  }
  
  # Data file may need adapted depending on platform
  if(Sys.info()[['sysname']]=="Darwin"){
    plat.specific.data.file <- gsub("/", "\\\\", data.file)
  }else{
    plat.specific.data.file <- data.file
  }
  
  # input the absolute path to the data file
  cat(paste("INFILE=", plat.specific.data.file, " /NOECHO;", sep=""),
      file=command.file.name, "\n", append=TRUE)
  cat("END;", file=command.file.name, "\n", append=TRUE)
  
  # ESTIMATE section
  # ~~~~~~~~~~~~~~~~
  cat("ESTIMATE;", file=command.file.name, "\n", append=TRUE)
  # we are only interested in the estimates for detection probability
  cat("DETECTION=ALL;", file=command.file.name, "\n", append=TRUE)
  
  # specify the key function
  cat("ESTIMATOR /KEY=", file=command.file.name, append=TRUE)
  if(mod.vals$key == "hn"){
    cat("HNORMAL", file=command.file.name, append=TRUE)
  }else if(mod.vals$key == "hr"){
    cat("HAZARD", file=command.file.name, append=TRUE)
  }else if(mod.vals$key == "unif"){
    cat("UNIFORM", file=command.file.name, append=TRUE)
  }else if(mod.vals$key == "gamma"){
    cat("NEXPON", file=command.file.name, append=TRUE)
  }
  
  adj.pres <- FALSE
  # check if adjustment terms are used
  if(is.null(mod.vals$adj.series) == FALSE){
    adj.pres <- TRUE
    # specify the type of adjustment term
    if(mod.vals$adj.series == "cos"){
      cat(" /ADJUST=COSINE", file=command.file.name, append=TRUE)
    }else if(mod.vals$adj.series == "herm"){
      cat(" /ADJUST=HERMITE", file=command.file.name, append=TRUE)
    }else if(mod.vals$adj.series == "poly"){
      cat(" /ADJUST=POLY", file=command.file.name, append=TRUE)
    }
    
    # specify the order of adjustment terms
    if(!is.null(mod.vals$adj.order)){
      cat(paste(" /ORDER=", paste(mod.vals$adj.order,collapse=","),sep=""),
          file=command.file.name, append=TRUE)
    }
    
    # specify the scaling of adjustment parameters
    if(mod.vals$adj.scale == "width"){
      cat(" /ADJSTD=W", file=command.file.name, append=TRUE)
    }else{
      cat(" /ADJSTD=SIGMA", file=command.file.name, append=TRUE)
    }
    
    # specify the number of adjustment parameters
    cat(paste(" /NAP=", length(mod.vals$adj.order), sep=""),
        file=command.file.name, append=TRUE)
  }
  
  # specify which fields are covariates
  if(length(covar.fields) > 0){
    cat(paste(" /COVARIATES=", paste(covar.fields,collapse=","), sep=""),
        file=command.file.name, append=TRUE)
  }
  
  # allowing for initial values for the parameters
  if(!(length(control$initial) == 1 && is.na(control$initial))){
    inits <- numeric()
    if(!mod.vals$key == "unif"){
      # first start value is the exp intercept scale param
      inits[1] <- exp(control$initial$scale[1])
    }
    if(mod.vals$key == "hr"){
      # second start value is the exp shape param
      inits[2] <- exp(control$initial$shape[1])
    }
    # Next are any covariates
    if(length(covar.fields) > 0){
      cov.inits <- control$initial$scale[-1]
      # Check if any are factors
      if(!is.null(factor.fields)){
        #Sort out ordering and may need to change scale intercept!
      }
      inits <- c(inits, cov.inits)
    }
    # If there are any adjustments
    if(!is.null(control$initial$adjustment)){
      inits <- c(inits, control$initial$adjustment)
    }
    # paste all the initial values together
    cat(paste(" /START=", paste(inits,collapse=","), sep=""),
        file=command.file.name, append=TRUE)
  }
  
  # defining upper and lower bounds for parameters
  #if(!is.null(control$lowerbounds)){
  #  cat(paste(" /LOWER=", paste(control$lowerbounds,collapse=","), sep=""),
  #      file=command.file.name, append=TRUE)
  #}
  #if(!is.null(control$upperbounds)){
  #  cat(paste(" /UPPER=", paste(control$upperbounds,collapse=","), sep=""),
  #      file=command.file.name, append=TRUE)
  #}
  
  # ending the ESTIMATOR line
  cat(";", file=command.file.name, "\n", append=TRUE)
  
  # specifying monotonicity constraint
  if(is.null(meta.data$mono.strict)){
    meta.data$mono.strict <- FALSE
  }
  if(is.null(meta.data$mono)){
    meta.data$mono <- FALSE
  }
  
  if(meta.data$mono.strict == TRUE){
    cat("MONOTONE=STRICT;", file=command.file.name, "\n", append=TRUE)
  }else if(meta.data$mono == TRUE){
    cat("MONOTONE=WEAK;", file=command.file.name, "\n", append=TRUE)
  }else{
    cat("MONOTONE=NONE;", file=command.file.name, "\n", append=TRUE)
  }
  
  # specifying the truncation distance
  cat(paste("DISTANCE /WIDTH=",meta.data$width,sep=""),
      file=command.file.name, append=TRUE)
  # dealing with grouped data
  if(is.null(meta.data$binned) == FALSE){
    if(meta.data$binned == TRUE){
      cat(paste(" /INTERVALS=", paste(meta.data$breaks, collapse=","),
                sep=""), file=command.file.name, append=TRUE)
    }
  }
  
  # specifying left trunction distance, if included in the input
  if(is.null(meta.data$left) == FALSE){
    cat(paste(" /LEFT=", meta.data$left, sep=""),
        file=command.file.name, append=TRUE)
  }
  cat(";", file=command.file.name, "\n", append=TRUE)
  cat("END;", file=command.file.name, append=TRUE)
  cat("\n", file=command.file.name, append=TRUE)
  
  ret <- list(command.file.name = command.file.name,
              stats.file.name = stat.file,
              log.file.name = log.file,
              out.file.name = out.file,
              data.file.name = data.file,
              plot.file.name = plot.file)
  
  return(ret)
}


#' @importFrom utils read.table write.table
mcds.results.and.refit <- function(statsfile, model.list, debug=FALSE){
  
  stats <- read.table(statsfile, row.names=NULL)
  # Name columns
  colnames(stats) <- c("Stratum", "Sample", "Estimator", "Module", "Statistic",
                       "Value", "CV", "Lcl", "Ucl", "Df")
  
  # Extract only parameters relating to the detection function (Module 2)
  stats <- stats[stats$Module == 2,]
  # Extract the log likelihood value (statistic 9)
  ll <- stats[stats$Statistic == 9,]$Value
  # Extract the estimated values of each parameter (statistics 100+)
  starting.values <- stats$Value[stats$Statistic > 100]
  
  # Put estimates into model as initial values
  model.list$control$initial <- list()
  mod.paste <- paste(model.list$dsmodel)
  mod.vals <- try(eval(parse(text=mod.paste[2:length(mod.paste)])))
  
  # As long as it's not a uniform there will be a scale parameter
  if(!mod.vals$key == "unif"){
    # Log scale value as mrds works on log scale for and MCDS works on natural scale for the scale parameter
    model.list$control$initial$scale <- log(starting.values[1])
    starting.values <- starting.values[-1]  
  }
  
  # If it's a hazard rate extract the shape parameter (now the 1st element)
  if(mod.vals$key == "hr"){
    model.list$control$initial$shape <- log(starting.values[1])
    starting.values <- starting.values[-1]
  }
  
  # Adjustment parameters come at the end
  if(!is.null(mod.vals$adj.order)){
    ind <- (length(starting.values)-length(mod.vals$adj.order)+1):
      length(starting.values)
    model.list$control$initial$adjustment <- starting.values[ind]
    starting.values <- starting.values[-ind]
  }
  
  # Any thing else must be covariate parameters associated with the scale parameter
  if(length(starting.values) >0){
    model.list$control$initial$scale <- c(model.list$control$initial$scale, starting.values)
  }
  
  # We are not refitting just using mrds to make the model with the parameter estimates from MCDS
  model.list$control$nofit <- TRUE
  refit <- ddf(dsmodel = model.list$dsmodel, mrmodel=NULL,
               data = model.list$data,
               method = model.list$method,
               meta.data = model.list$meta.data,
               control = model.list$control)
  
  # For debugging
  if(debug){
    message(paste0("MCDS.exe log likehood: ", round(ll,7)))
    message(paste0("MCDS.exe pars: ", paste(round(stats$Value[stats$Statistic > 100],7), collapse=", ")))
    message(paste0("mrds refitted log likehood: ", round(refit$lnl,7)))
    message(paste0("mrds refitted pars: ", paste(round(refit$par, 7), collapse=", ")))
  }
  
  return(refit)
}

run.MCDS <- function(dsmodel, data, method, meta.data, control){
  
  if(control$showit>0){
    message("Running MCDS.exe...")
  }
  # create the test file
  test.file <- create.command.file(dsmodel, data, method,
                                   meta.data, control)
  # When exiting the function tidy up any files that may have been created
  on.exit(tidy.tmp.files(unlist(test.file)))
  
  if(control$showit>0){
    message(paste("Command file written to", test.file$command.file.name))
    message(paste("Stats file written to", test.file$stats.file.name))
  }
  
  # find MCDS.exe
  path.to.MCDS.dot.exe <- system.file("MCDS.exe", package="mrds")
  
  # actually run MCDS.exe (use wine if necessary)
  if(Sys.info()[['sysname']]!="Windows"){
    if(!is.null(control$winebin)){
      winebin <- control$winebin
    }else{
      winebin <- "wine"
    }
    # macOS specific
    if(Sys.info()[['sysname']]=="Darwin"){
      osx.version <- as.numeric_version(system("sw_vers -productVersion",
                                               intern = TRUE))
      if(!is.null(control$winebin)){
        winebin <- control$winebin
      }else if(osx.version>as.numeric_version(10.15)){
        # Catalina or after removed the 32-bit support, so use wine32on64
        winebin <- "wine32on64"
      }
    }
    if(Sys.which(winebin)==""){
      stop("MCDS.exe could not be run. wine (needed to run MCDS.exe on non-Windows platforms) could not automatically be detected. See documentation for details.", call. = FALSE)
    }
    wine.call <- paste(winebin, path.to.MCDS.dot.exe, "0,", test.file$command.file.name)
    # only provide the output from MCDS.exe when we have showit>0
    w <- suppressWarnings(try(system(wine.call, intern=TRUE,
                ignore.stdout=control$showit<1, ignore.stderr=control$showit<1)
                , silent = TRUE))
  }else{
    # on Windows just execute the MCDS binary
    w <- suppressWarnings(try(system(paste0(path.to.MCDS.dot.exe,
                                            " 0, ", test.file$command.file.name), intern=TRUE,
                                     ignore.stdout=control$showit<1, ignore.stderr=control$showit<1, show.output.on.console = FALSE)
                              , silent = TRUE))
  }
  
  # Check if MCDS.exe was successful
  finish.status <- attr(w, "status")
  
  # Read in the log file
  log.file <- readLines(test.file$log.file.name)
  
  # Check for warnings or errors
  warn.index <- grep("[Ww]arning", log.file)
  err.index <- grep("[Ee]rror", log.file)
  
  # If the finish.status is not 1 and there are no warnings or errors
  if(finish.status != 1 && length(warn.index) == 0 && length(err.index) ==0){
    warning(paste("MCDS finish status indicates warnings / error but none were found in the log file. Finish status = ", finish.status, sep = ""), call. = FALSE, immediate. = TRUE)
  }
  # Display warnings from log file
  for(i in seq(along = warn.index)){
    ignore.warning <- any(grepl("Cannot estimate encounter rate variance empirically as there is only one sample.", log.file[warn.index[i]]),
                          grepl("When cluster size is a covariate, variance of the cluster size, density of individuals, and abundance estimates can only be obtained via the bootstrap.", log.file[warn.index[i]]))
    if(!ignore.warning) message(log.file[warn.index[i]], appendLF = TRUE)
  }
  # Display errors from log file
  for(i in seq(along = err.index)){
    message(log.file[err.index[i]], appendLF = TRUE)
  }
  # If there are any errors return a value of try-error
  if(length(err.index) > 0){
    ret.val <- list()
    class(ret.val) <- "try-error"
    return(ret.val)
  }
  
  # little extra parameter here to avoid infinite recursion when
  # we run ddf in mcds.results.and.refit
  control$optimizer <- "R"
  model.list <- list(dsmodel=dsmodel, data=data,
                     method=method, meta.data=meta.data, control=control)
  
  # try to refit the model
  mm <- try(mcds.results.and.refit(test.file$stats.file.name, model.list,
                                   debug=control$showit>0), silent = TRUE)
  
  return(mm)
}

# Tidy up temporary files
tidy.tmp.files <- function(files = character()){
  # Try to remove each of the files
  for(i in seq(along = files)){
    # Quietly try and remove files that were created
    tmp <- suppressWarnings(try(file.remove(files[i]), silent = TRUE))
  }
}

#' @importFrom utils write.table
create.data.file <- function(data, dsmodel, data.file){
  
  # Rename columns to those expected by MCDS.exe
  data.cols <- c("SMP_LABEL", "SMP_EFFORT", "DISTANCE")
  
  # SMP_LABEL (sample label)
  if("Sample.Label" %in% names(data)){
    index <- which(names(data) == "Sample.Label")
    names(data)[index] <- "SMP_LABEL"
  }else{
    data$SMP_LABEL <- rep(1, nrow(data))
  }
  
  # SMP_EFFORT (sample effort)
  if("Effort" %in% names(data)){
    index <- which(names(data) == "Effort")
    names(data)[index] <- "SMP_EFFORT"
  }else{
    data$SMP_EFFORT <- rep(1, nrow(data))
  }
  
  # DISTANCE
  index <- which(names(data) == "distance")
  names(data)[index] <- "DISTANCE"
  
  # STR_LABEL (stratum label)
  if("Region.Label" %in% names(data)){
    # Replicate rather than rename in case it is included as a covariate
    data$STR_LABEL <- data$Region.Label
    data.cols <- c(data.cols, "STR_LABEL")
  }
  
  # STR_AREA (stratum area)
  if("Area" %in% names(data)){
    index <- which(names(data) == "Area")
    names(data)[index] <- "STR_AREA"
    data.cols <- c(data.cols, "STR_AREA")
  }
  
  # Check if there are cluster sizes
  # Let's not tell MCDS about clusters - not needed for detection function fitting
  cluster <- FALSE
  # if("size" %in% names(data)){
  #   index <- which(names(data) == "size")
  #   names(data)[index] <- "SIZE"
  #   data.cols <- c(data.cols, "SIZE")
  #   cluster <- TRUE
  # }else{
  #   cluster <- FALSE
  # }
  
  # Check if there are any covariates in the model
  if(length(all.vars(dsmodel)) > 0){
    # Get model covariates
    covar.fields <- all.vars(dsmodel)
    data.cols <- c(data.cols, covar.fields)
  }else{
    covar.fields <- NULL
  }
  
  # remove all non-essential columns from the dataset
  data <- data[,data.cols]
  
  # Cluster Size must be named "SIZE" - but we do not want it to know about clusters so name "CS" instead
  if("size" %in% names(data)){
    index <- which(names(data) == "size")
    names(data)[index] <- "CS"
  }
  if("size" %in% covar.fields){
    index <- which(covar.fields == "size")
    covar.fields[index] <- "CS"
  }
  
  # create data file to pass to mcds
  file.create(data.file)
  write.table(data, file=data.file, col.names=FALSE,
              row.names=FALSE, sep="\t", quote = FALSE)
  
  # return reformatted dataframe and the information about covariates
  output <- list(data = data, 
                 covar.fields = covar.fields,
                 cluster = cluster)
  return(output)
}


validate.MCDS.options <- function(model){
  # Validates the options are suitable for MCDS
  
  # Extract values from the formula
  mod_paste <- paste(model)
  mod.vals <- try(eval(parse(text=mod_paste[2:length(mod_paste)])))
  
  # Check if there are covariates
  if(mod.vals$formula != ~1){
    # if so only HN and HR permitted
    if(!mod.vals$key %in% c("hn", "hr")){
      return("MCDS can only fit covariates with the half-normal and hazard rate key functions.")
    }
  }
  
  # If no issue return NULL
  return(NULL)
}

