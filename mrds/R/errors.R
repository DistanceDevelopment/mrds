#' Error function
#'
#' Writes error messages for various errors that can occur
#'
#'
#' @param errmsg the message to be stored/printed (optional)
#' @param errmode report or print errors (default report)
#' @param preamble character string to paste before the message 
#' @return None
#' @author Dave Miller
errors <- function(errmsg=NULL,errmode="report",preamble="Warning"){

  if(errmode == "report"){
    cat(paste("\n** ",preamble,": ",errmsg,"**\n",sep=""))
  }else{
    cat("\nOnly report is implemented at the moment.\n")
  }
}
