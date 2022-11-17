#' Run MCDS.exe as a backend for mrds
#'
#' Rather than use the R code provided in `mrds`, one can also use the binary
#' of `MCDS.exe`, to reproduce the results given by Distance for Windows. There
#' is no guarantee that one approach is "better" than the other, but `mrds`
#' will select the model with the better likelihood and provide answers to
#' this.
#'
#' If you are running a non-Windows operating system, you can follow the
#' instructions below to have `MCDS.exe` run using `wine`.
#'
#' @aliases MCDS
#' @rdname mcds-dot-exe
#'
#' @section Obtaining MCDS.exe:
#' The function `download_MCDS_dot_exe` can be used to download `MCDS.exe` from
#' the distance sampling website. The MCDS binary will be installed to your
#' local R library, in the `inst/` directory of the `mrds` directory. You can
#' also copy the `MCDS.exe` from your local Distance for Windows installation
#' if you prefer.
#'
#' @section Running MCDS.exe on non-Windows platforms:
#' One can still use MCDS.exe even if you are not running a Windows computer. To
#' do this one will need to install `wine` a Windows emulator. It is important
#' to use a version of `wine` which can run 32-bit programs.
#'
#' The package will attempt to work out which `wine` binary to use (and detect
#' if it is installed), but this doesn't always work. In this case, the
#' location of the `wine` binary can be specified in the `control` `list`
#' provided to `ddf` using the `winebin` element. For example, if `wine` is
#' installed at `/usr/bin/local/wine` you can set `control$winebin` to that
#' location to use that binary.
#'
#' On Linux, one needs to install a `wine` that includes 32 bit support. On
#' Debian, that can be provided by the package `wine32`.
#'
#' On macOS, this can be achieved using the `homebrew` package management
#' system and installing the `wine-crossover` package. You may need to change
#' the `control$winebin` to be `wine`, `wine64` or `wine32on64`, depending on
#' your system's setup. This package tries to work out what to do, but likely
#' doesn't handle all corner cases. Currently this is untested on Mac M1
#' systems, please get in touch if you are trying this. It seems like this does
#' not currently work on post-Catalina versions of macOS, though this will
#' hopefully be fixed in future versions of `wine`.
#'
#' @section: Stopping using MCDS.exe
#' Once this feature is enabled (see below) using `ddf` will always run both
#' its built-in R optimizer and `MCDS.exe`. To disable this behaviour remove
#' the `MCDS.exe` binary file. You can find it by running the following in R:
#' `system.file("MCDS.exe", package="mrds")`.
#'
#' @author David L Miller and Jonah McArthur



# a function to take the input as required for the ddf function in mrds and
# create a command file for the mcds engine, which will perform an equivalent
# analysis.
# the input for this function is the same as the input for ddf
# written by Jonah McArthur
create_command_file <- function(dsmodel=call(), mrmodel=call(), data,
                                method, meta.data, control) {
  # create a temporary directory
  directory <- tempdir()

  # create command file
  command.file.name <- tempfile(pattern="cmdtmp", tmpdir=directory,
                                fileext=".txt")
  #command.file.name <- gsub("/","\\\\",command.file.name)
  file.create(command.file.name)

  # HEADER section

  # specify the location of output files
  cat(paste(directory,"/out.txt",sep=""), file=command.file.name, "\n",
      append=TRUE)
  cat(paste(directory,"/log.txt",sep=""), file=command.file.name, "\n",
      append=TRUE)
  cat(paste(directory,"/stat.txt",sep=""), file=command.file.name,
      "\n", append=TRUE)
  cat(paste(directory,"/plot.txt",sep=""), file=command.file.name,
      "\n", append=TRUE)
  cat("None", file=command.file.name, "\n", append=TRUE)
  cat("None", file=command.file.name, "\n", append=TRUE)

  # removing irrelevant data
  # combine data from multiple observers
  data <- data[data$detected==1,]
  if(TRUE %in% grepl("^object$",tolower(colnames(data)))){
    # in case the object column isn't named 'object', to remove case sensitivity
    obj_col <- grep("^object$",tolower(colnames(data)))
    colnames(data)[obj_col] <- "object"
    # identifying all objects and taking the first data point for each
    obj_nums <- unique(data$object)
    for(i in 1:length(obj_nums)){
      entries <- grep(TRUE,data$object==i)
      if(length(entries)>1){
        remove <- entries[-1]
        data <- data[-remove,]
      }
    }
  }

  # create a vector of required fields
  req_fields <- c("SMP_LABEL","SMP_EFFORT","DISTANCE")

  # change the distance column name to upper case
  colnames(data)[grep("^distance$",tolower(colnames(data)))] <- "DISTANCE"

  # find which field will be used for the effort and change the name
  # to match field name in mcds
  if(TRUE %in% grepl("SMP_EFFORT",toupper(colnames(data)))){
    colnames(data)[grep("^SMP_EFFORT",toupper(colnames(data)))] <- "SMP_EFFORT"
  }else if(TRUE %in% grepl("^effort$",tolower(colnames(data)))){
    colnames(data)[grep("^effort$",tolower(colnames(data)))] <- "SMP_EFFORT"
  }else if(TRUE %in% grepl("^Search.time$",colnames(data))){
    # !this may be a bit too specific to the example data in Distance
    colnames(data)[grep("^Search.time$",colnames(data))] <- "SMP_EFFORT"
  }else{
    data$SMP_EFFORT <- rep(1,nrow(data))
  }

  # find if SMP_LABEL is a field; if not, add it
  if(TRUE %in% grepl("^SMP_LABEL",toupper(colnames(data)))){
    colnames(data)[grep("^SMP_LABEL",toupper(colnames(data)))] <- "SMP_LABEL"
  }else if(TRUE %in% grepl("^sample.label$",tolower(colnames(data)))){
    colnames(data)[grep("^sample.label$",tolower(colnames(data)))] <- "SMP_LABEL"
  }else{
    data$SMP_LABEL <- rep(1,nrow(data))
  }

  # check if other defined fields are columns in the dataset; if they are add
  # them to the list of fields to keep in the dataset
  if(TRUE %in% grepl("^STR_LABEL$",toupper(colnames(data)))){
    colnames(data)[grep("^STR_LABEL$",toupper(colnames(data)))] <- "STR_LABEL"
    req_fields <- append(req_fields,"STR_LABEL")
  }else if(TRUE %in% grepl("^region.label$",tolower(colnames(data)))){
    colnames(data)[grep("^region.label$",tolower(colnames(data)))] <- "STR_LABEL"
    req_fields <- append(req_fields,"STR_LABEL")
  }
  if(TRUE %in% grepl("^STR_AREA$",toupper(colnames(data)))){
    colnames(data)[grep("^STR_AREA$",toupper(colnames(data)))] <- "STR_AREA"
    req_fields <- append(req_fields,"STR_AREA")
  }else if(TRUE %in% grepl("^area$",tolower(colnames(data)))){
    colnames(data)[grep("^area$",tolower(colnames(data)))] <- "STR_AREA"
    req_fields <- append(req_fields,"STR_AREA")
  }
  if(TRUE %in% grepl("^size$",tolower(colnames(data)))){
    cluster <- TRUE
    colnames(data)[grep("^size$",tolower(colnames(data)))] <- "SIZE"
    req_fields <- append(req_fields,"SIZE")
  }

  # specifying covariates in the model
  if(identical(all.vars(dsmodel),character(0)) == FALSE){
    covar_pres <- TRUE
    # extracting the list of covariates
    covars <- all.vars(dsmodel)
    # creating a list of the fields for each covariate
    covar_fields <- rep("",length(covars))
    for(i in 1:length(covars)){
      # identifying the field name for each covariate
      index <- grep(covars[i],tolower(colnames(data)))
      covar_fields[i] <- colnames(data)[index]
    }
    # the required fields cannot be covariates in the model, with the exception of size
    if(length(intersect(tolower(req_fields),tolower(covar_fields))) > 0){
      if(TRUE %in% grepl("size",tolower(covar_fields))){
        # specify whether SIZE is a covariate
        size_cov <- TRUE
      }
      # remove any required fields from the list of covariates
      covar_fields <- covar_fields[! covar_fields %in% intersect(req_fields,covar_fields)]
    }else{
      size_cov <- FALSE
    }
    # add covariates to the fields that are kept for analysis
    req_fields <- c(req_fields,covar_fields)
    # if SIZE is a covariate, add it back to the list of covariates
    if(size_cov == TRUE){
      covar_fields <- append(covar_fields,"SIZE")
    }
  }else{
    covar_pres <- FALSE
  }

  # remove all non-essential columns from the dataset
  data <- data[req_fields]

  # create data file to pass to mcds
#  data.file.name <- tempfile(pattern="data", tmpdir=directory,
#                             fileext=".txt")
  data.file.name <- paste(directory,"/data.txt",sep="")
  #data.file.name <- gsub("/","\\\\",data.file.name)
  file.create(data.file.name)
  write.table(data, file=data.file.name, col.names=FALSE,
              row.names=FALSE, sep="\t")

  # OPTION section

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

  cat("DATA /STRUCTURE=FLAT;", file=command.file.name, "\n",
      append=TRUE)

  # change all fields to upper case and combine to one string
  fields_comb <- paste(toupper(colnames(data)), collapse=",")
  cat(paste("FIELDS=", fields_comb, ";", sep=""),
      file=command.file.name, "\n", append=TRUE)

  # specifying which fields are factor covariates
  if(covar_pres == TRUE){
    factor_fields <- c()
    for(i in 1:length(colnames(data))){
      # for each covariate, check if it is a factor
      if(is.factor(data[,i]) && (TRUE %in% grepl(colnames(data)[i],covar_fields))){
        # if the covariate is a factor, specify its name, levels, and level labels
        factor_fields <- append(factor_fields,colnames(data)[i])
        labels <- paste(levels(data[,i]), collapse=",")
        cat(paste("FACTOR /NAME=", toupper(colnames(data)[i]),
                  " /LEVELS=", length(levels(data[,i])), " /LABELS=",
                  labels, sep=""), file=command.file.name, "\n",
            append=TRUE)
      }
    }
  }

  # input the absolute path to the data file
  data.file.name2 <- gsub("/", "\\\\", data.file.name)
  cat(paste("INFILE=", data.file.name2, " /NOECHO;", sep=""),
      file=command.file.name, "\n", append=TRUE)
  cat("END;", file=command.file.name, "\n", append=TRUE)

  # ESTIMATE section

  cat("ESTIMATE;", file=command.file.name, "\n", append=TRUE)

  # we are only interested in the estimates for detection probability
  cat("DETECTION ALL;", file=command.file.name, "\n", append=TRUE)

  # a messy way of accessing the model parameters
  mod_paste <- paste(dsmodel)
  mod_vals <- try(eval(parse(text=mod_paste[2:length(mod_paste)])))

  # specify the key function
  cat("ESTIMATOR /KEY=", file=command.file.name, append=TRUE)
  if(mod_vals$key == "hn"){
    cat("HNORMAL", file=command.file.name, append=TRUE)
  }else if(mod_vals$key == "hr"){
    cat("HAZARD", file=command.file.name, append=TRUE)
  }else if(mod_vals$key == "unif"){
    cat("UNIFORM", file=command.file.name, append=TRUE)
  }else{
    cat("NEXPON", file=command.file.name, append=TRUE)
  }

    adj_pres <- FALSE
  # check if adjustment terms are used
  if(is.null(mod_vals$adj.series) == FALSE){
    adj_pres <- TRUE
    # specify the type of adjustment term
    if(mod_vals$adj.series == "cos"){
      cat(" /ADJUST=COSINE", file=command.file.name, append=TRUE)
    }else if(mod_vals$adj.series == "herm"){
      cat(" /ADJUST=HERMITE", file=command.file.name, append=TRUE)
    }else if(mod_vals$adj.series == "poly"){
      cat(" /ADJUST=POLY", file=command.file.name, append=TRUE)
    }

    # specify the order of adjustment terms
    if(!is.null(mod_vals$adj.order)){
      cat(paste(" /ORDER=", paste(mod_vals$adj.order,collapse=","),sep=""),
          file=command.file.name, append=TRUE)
    }

    # specify the scaling of adjustment parameters
    if(mod_vals$adj.scale == "width"){
      cat(" /ADJSTD=W", file=command.file.name, append=TRUE)
    }else{
      cat(" /ADJSTD=SIGMA", file=command.file.name, append=TRUE)
    }

    # specify the number of adjustment parameters
    cat(paste(" /NAP=", length(mod_vals$adj.order), sep=""),
        file=command.file.name, append=TRUE)
  }

  # specify which fields are covariates
  if(covar_pres == TRUE){
    cat(paste(" /COVARIATES=", paste(covar_fields,collapse=","), sep=""),
        file=command.file.name, append=TRUE)
  }

  # allowing for initial values for the parameters
  inits <- c()
  if(is.null(control$initial) == FALSE){
    # go through covariates in order, if they are present
    if(covar_pres == TRUE){
      for(i in 1:length(covars)){
        # find the index of the covariate field
        index <- grep(toupper(covar_fields[i]),toupper(colnames(data)))
        # if the covariate is a factor, initial values must be given for each level
        if(TRUE %in% grepl(covar_fields[i],factor_fields)){
          for(j in 2:length(levels(data[,index]))){
            # create the text for the parameter that must be accessed
            access_covar <- paste("control$initial$scale$",
                                  colnames(data)[index],"[",j,"]",sep="")
            # evaluate the text in order to access the initial value
#            inits <- append(inits,eval(parse(text=access_covar)))
          }
          # the first level has to be last in MCDS
          access_covar <- paste("control$initial$scale$",
                                colnames(data)[index],"[1]",sep="")
#          inits <- append(inits,eval(parse(text=access_covar)))
        }else{
          access_covar <- paste("control$initial$scale$",
                                colnames(data)[index],sep="")
#          inits <- append(inits,eval(parse(text=access_covar)))
        }
      }
    }
    # add in shape parameter if hazard-rate used
    if(mod_vals$key == "hr"){
#      inits <- append(inits,control$initial$shape)
    }
    # add in adjustment initial values
    if(adj_pres){
      for(i in 1:length(mod_vals$adj.order)){
#        inits <- append(inits,control$initial$adjustment[i])
      }
    }
    # paste all the initial values together
#    cat(paste(" /START=", paste(inits,collapse=","), sep=""),
#        file=command.file.name, append=TRUE)
  }

  # defining upper and lower bounds for parameters
  if(is.null(control$lowerbounds) == FALSE){
    cat(paste(" /LOWER=", paste(control$lowerbounds,collapse=","), sep=""),
        file=command.file.name, append=TRUE)
  }
  if(is.null(control$upperbounds) == FALSE){
    cat(paste(" /UPPER=", paste(control$upperbounds,collapse=","), sep=""),
        file=command.file.name, append=TRUE)
  }

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

  ret <- list(command.file.name=command.file.name,
              stats.file.name = paste0(directory, "/stat.txt"))

  return(ret)
}

mcds_results_and_refit <- function(statsfile, model_list, debug=FALSE){

  stats <- read.table(statsfile, row.names=NULL)
  colnames(stats) <- c("Stratum", "Sample", "Estimator", "Module", "Statistic",
                       "Value", "CV", "Lcl", "Ucl", "Df")

  # based on the cuniform tablets, we need
  # module 2
  # statistics:
  #   9 -- log likelihood
  #   101+ -- parameter estimates

  stats <- subset(stats, Module==2)
  ll <- subset(stats, Statistic==9)
  stats <- subset(stats, Statistic>100)
  stats$Name <- ""
  stats$Name[stats$Statistic>100] <- paste0("A(",
                                       stats$Statistic[stats$Statistic>100]-100,
                                       ")")
  stats <- stats[, c("Name", "Value")]

  starting_values <- stats$Value


  model_list$control$initial <- list()
  mod_paste <- paste(model_list$dsmodel)
  mod_vals <- try(eval(parse(text=mod_paste[2:length(mod_paste)])))

  if(mod_vals$key == "hr"){
    model_list$control$initial$shape <- log(starting_values[2])
    starting_values <- starting_values[-2]
  }
  # MCDS.exe exponentiates the scale par, undo that
  model_list$control$initial$scale <- log(starting_values[1])
  starting_values <- starting_values[-1]

  # get the adjustments off the end
  if(!is.null(mod_vals$adj.order)){
    ind <- (length(starting_values)-length(mod_vals$adj.order)+1):
            length(starting_values)
    model_list$control$initial$adjustment <- starting_values[ind]
    starting_values <- starting_values[-ind]
  }

  if(length(starting_values) >0){
    model_list$control$initial$scale <- c(model_list$control$initial$scale,
                                          starting_values)
  }

  model_list$control$nofit <- TRUE

  refit <- ddf(dsmodel = model_list$dsmodel, mrmodel=NULL,
               data = model_list$data,
               method = model_list$method,
               meta.data = model_list$meta.data,
               control = model_list$control)

  if(debug){
    message(paste0("MCDS.exe log likehood: ", refit$lnl))
    message(paste0("MCDS.exe pars: ", paste(refit$par, collapse=", ")))
  }

  return(refit)
}

#' Install the MCDS.exe binary from the web
#'
#' This function will download the `MCDS.exe` binary from the distance sampling
#' website and install it.
#'
#' If you already have a copy of `MCDS.exe` on your computer you can move it to
#' the correct directory manually. You can find where this is by running the
#' following: `system.file(package="mrds")`.
#'
#' @return Nothing, just downloads and installs the binary.
#' @export
download_MCDS_dot_exe <- function(){

  url <- "http://distancesampling.org/R/"

  message("Downloading MCDS.exe...")
  success <- download.file(url, paste0(system.file(package="mrds"),"/MCDS.exe"))

  if(success==0){
    message(paste0("MCDS.exe installed at ",
            system.file(package="mrds"), "/MCDS.exe"))
  }else{
    stop("Download and installation not successful!")
  }
}

run_MCDS <- function(dsmodel, mrmodel, data, method, meta.data, control){

  # create the test file
  test_file <- create_command_file(dsmodel, mrmodel, data, method,
                                   meta.data, control)

  # find MCDS.exe
  path_to_MCDS_dot_exe <- system.file("MCDS.exe", package="mrds")

  # actually run MCDS.exe (use wine if necessary)
  if(Sys.info()[['sysname']]!="Windows"){
    if(!is.null(control$winebin)){
      winebin <- control$winebin
    }else{
      winebin <- "wine"
    }
    # macOS specific
    if(Sys.info()[['sysname']]=="Darwin"){
      osx_version <- as.numeric_version(system("sw_vers -productVersion",
                                               intern = TRUE))
      if(!is.null(control$winebin)){
        winebin <- control$winebin
      }else if(osx_version>as.numeric_version(10.15)){
        # Catalina or after removed the 32-bit support, so use wine32on64
        winebin <- "wine32on64"
      }
    }
    if(Sys.which(winebin)==""){
      stop("wine is needed to run MCDS.exe on non-Windows platforms. See documentation for details.")
    }
    wine_call <- paste0(winebin, " MCDS.exe 0, ", test_file$command.file.name)
    # only provide the output from MCDS.exe when we have showit>0
    w <- system(wine_call, intern=TRUE,
                ignore.stdout=control$showit>0, ignore.stderr=control$showit>0)

  }else{
    w <- system(paste0(path_to_MCDS_dot_exe,
                " 0, ", test_file$command.file.name), intern=TRUE,
                stdout=control$showit>0, stderr=control$showit>0)
  }

  #if(w>0){
  #  stop("Running MCDS.exe failed!")
  #}

  model_list <- list(dsmodel=dsmodel, mrmodel=mrmodel, data=golftees,
                     method=method, meta.data=meta.data, control=control)

  if(is.null(control$showit)){
    control$showit <- FALSE
  }
  mm <- mcds_results_and_refit(test_file$stats.file.name, model_list,
                               debug=control$showit)

}
