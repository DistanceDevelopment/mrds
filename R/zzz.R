.onAttach<-function(...){
  # uses packageStartupMessage which can then be
  # surpressed

  version <- utils::packageVersion("mrds")
  built <- utils::packageDescription("mrds",fields="Built")

  hello <- paste0("This is mrds ",version,"\nBuilt: ",built)
  
  packageStartupMessage(hello)
}
