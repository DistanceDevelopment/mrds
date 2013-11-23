.onAttach<-function(...){
  if (!interactive()) return()
  # this now conforms with new R conventions
  # uses packageStartupMessage which can then be
  # surpressed

  version <- packageVersion("mrds")
  built <- packageDescription("mrds",fields="Built")

  hello <- paste("This is mrds ",version,"\nBuilt: ",built,sep="")
  packageStartupMessage(hello)
}
