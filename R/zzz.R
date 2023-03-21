.onAttach<-function(...){
  # uses packageStartupMessage which can then be
  # surpressed

  version <- utils::packageVersion("mrds")
  built <- utils::packageDescription("mrds",fields="Built")

  hello <- paste0("This is mrds ",version,"\nBuilt: ",built)
  
  mcds.info <- ifelse(system.file("MCDS.exe", package="mrds") == "",
                      "MCDS.exe not detected, single observer analyses will only be run using optimiser in mrds R library. See ?MCDS for details.",
                      "MCDS.exe detected, by default single observer analyses will utilise both the mrds R optimiser and the MCDS.exe fortran optimiser to achieve the best fit. See ?ddf and ?MCDS for details.")
  
  hello <- paste0(hello, "\n", mcds.info)
  
  packageStartupMessage(hello)
}
