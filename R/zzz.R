.onAttach <- function(lib, pkg)
{
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  if(interactive())
  {
    #Rfiglet::figlet(message = "raedda", font = "roman")
    packageStartupMessage(
      "
                                   .o8        .o8
                                   888        888
oooo d8b  .oooo.    .ooooo.   .oooo888   .oooo888   .oooo.         Robust and Adaptive
 888  8P 'P  )88b  d88   88b d88   888  d88   888  'P   )88b    Eigenvalue Decomposition
 888      .oP 888  888ooo888 888   888  888   888    .oP 888      Discriminant Analysis
 888     d8(  888  888    .o 888   888  888   888   d8(  888
d888b    'Y888  8o  Y8bod8P   Y8bod88P   Y8bod88P   'Y888  8o

      ", "version ", version, "\n" )
  }
  else
  { packageStartupMessage("Package 'raedda' version ", version) }

  packageStartupMessage("Type 'citation(\"raedda\")' for citing this R package in publications.")
  invisible()
}
