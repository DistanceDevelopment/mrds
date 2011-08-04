#' Error function
#' 
#' Writes error messages for various errors that can occur
#' 
#' 
#' @param errmsg the message to be stored/printed (optional)
#' @param errmode report or print errors (default report)
#' @return None
#' @author Dave Miller
errors <- function(errmsg=NULL,errmode="report")
#
# errors
#
# A nice, neat way of storing and printing errors/warnings.
#
# Arguments:
# 
#  errmode	- report or print errors (default report)
#  errmsg	- the message to be stored/printed (optional)
#
# Values:
#
#  We return either a formatted error message, TRUE or a list()
#  of all the errors.
#
# dlm 25-Aug-05  Initial work started. At the moment we are just
#		 able to format the data, nothing more complicated
#		 yet.
#		 By default in future it should store the errors in
#		 a global variable.
#
{

  if(errmode == "report"){

    cat(paste("\n** Warning:",errmsg,"**\n"))
  
  }else{
    cat("\nOnly report is implemented at the moment.\n")
  }

}
