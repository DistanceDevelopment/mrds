#' Error function
#'
#' Writes error messages for various errors that can occur
#'
#'
#' @param errmsg the message to be stored/printed (optional)
#' @param preamble character string to paste before the message 
#' @return None
#' @author Dave Miller
errors <- function(errmsg=NULL,preamble="Warning"){
  cat(paste("\n** ",preamble,": ",errmsg,"**\n",sep=""))
}
