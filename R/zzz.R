.onAttach<-function(...){
  # uses packageStartupMessage which can then be
  # surpressed

  version <- utils::packageVersion("mrds")
  built <- utils::packageDescription("mrds",fields="Built")

  hello <- paste0("This is mrds ",version,"\nBuilt: ",built)
  
  varest.info <- "**Change to default variance estimator for point transects. The default encounter rate variance estimator for point transects is now 'P2' (changed from 'P3'). See 'Uncertainty' section of ?dht for more information.**"
  
  mcds.info <- ifelse(system.file("MCDS.exe", package="mrds") == "",
                      "MCDS.exe not detected, single observer analyses will only be run using optimiser in mrds R library. See ?MCDS for details.",
                      "MCDS.exe detected, by default single observer analyses will utilise both the mrds R optimiser and the MCDS.exe fortran optimiser to achieve the best fit. See ?ddf and ?MCDS for details.")
  
  hello <- paste0(hello, "\n\n", varest.info, "\n\n", mcds.info)
  
  packageStartupMessage(hello)
}
