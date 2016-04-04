#--------------------------------------------------------------------------------------------------#
#
#   C3 and C4 model representations of photosynthesis.  Used for fitting and deriving
#   key photosynthetic parameters
#
#--------------------------------------------------------------------------------------------------#


# ============================================== C3 ================================================ 


#--------------------------------------------------------------------------------------------------#
##' 
##' Farquhar functions
##' 
##' @author Shawn P. Serbin
##' 
# Vcmax, Rd, and Jmax
ACi.VJRd <- function(x) {
  Vcmax <- x[1]  ## Vcmax param
  Rd    <- x[2]  ## Resp param
  Jmax  <- x[3]  ## Jmax
  
  Ac <- ifelse(Ci<=Vcmax.cutoff,(Vcmax*(Ci-Gstar)/(Ci + Km))-Rd,9999)
  inter <- ifelse((Ci>Vcmax.cutoff & Ci<Jmax.cutoff),Photo,9999)
  Aj <- ifelse(Ci>=Jmax.cutoff,(Jmax*(Ci-Gstar)/((4.5*Ci)+(10.5*Gstar)))-Rd,9999)
  #Aj <- ifelse(Ci>=Jmax.cutoff,(Jmax*(Ci-Gstar)/((4.5*Ci)+(10.5*Gstar))),9999)
  Assim <- pmin(Ac,Aj,inter)
  RMSE <- sqrt(mean((Photo-Assim)^2))  ## RMSE cost function
  return(RMSE)
}
#==================================================================================================#