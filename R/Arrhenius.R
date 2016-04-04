####################################################################################################
#
#  	Author: Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Standard Arrhenius temperature function
arrhenius <- function(Tleaf.1,Tleaf.2,VJ,Ev) {
  R <- 0.008314472
  VJ.T2 <- VJ*exp((Ev*((Tleaf.2+273.15)-(Tleaf.1+273.15)))/((Tleaf.1+273.15)*R*(Tleaf.2+273.15)))
  return(VJ.T2)
}

### Apply Arrhenius Temperature Scaling
Tleaf.1 <- 13.24283482                    # Temp of Vcmax at leaf temp 1
Tleaf.2 <- 25                             # Temp to scale Vcmax to from Tleaf.1
VJ <- 55.385                              # Value of Vcmax/Jmax at Tleaf.1
Ev <- 33.861434                           # Previously determined activation energy (E)
arrhenius(Tleaf.1,Tleaf.2,VJ,Ev)
#--------------------------------------------------------------------------------------------------#
### EOF