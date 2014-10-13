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
# Rd and Vcmax only
ACi.rubisco <- function(x) {
  keep <- which(Ci<Vcmax.cutoff)
  Vcmax <- x[1]  ## Vcmax param
  Rd    <- x[2]  ## Resp param
  
  Ac    <- (Vcmax*(Ci[keep]-Gstar[keep]))/(Ci[keep]+Km[keep])
  Assim <- pmin(Ac)-Rd
  RMSE <- sqrt(mean((Photo[keep]-Assim)^2))  ## RMSE cost function
  return(RMSE)
}
#==================================================================================================#


#--------------------------------------------------------------------------------------------------#
# Rd and Jmax only
ACi.rubp  <- function(x) {
  keep <- which(Ci>Jmax.cutoff)
  Jmax <- x[1]  ## Jmax param
  Rd   <- x[2]  ## Resp param
  
  Aj    <- Jmax*(Ci[keep]-Gstar[keep])/((4.5*Ci[keep])+(10.5*Gstar[keep]))
  Assim <- pmin(Aj)-Rd
  RMSE <- sqrt(mean((Photo[keep]-Assim)^2))  ## RMSE cost function
  return(RMSE)
}
#==================================================================================================#


#--------------------------------------------------------------------------------------------------#
# Vcmax, Rd, and Jmax
ACi.full <- function(x) {
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


# ============================================== C4 ================================================ 


#--------------------------------------------------------------------------------------------------#
##' 
#==================================================================================================#


# ========================================= Light Response ========================================= 


#--------------------------------------------------------------------------------------------------#
# A-PPFD parameters
A.Q <- function(x) { 
  ## Setup variables to solve for
  Amax  <- x[1]  # Amax param.  Light saturated rate of photo
  convex <- x[2]  # Convexity of the light response curve
  
  # Rectangular hyper
  #mod.Photo <- ((IaQY+Amax)-sqrt((IaQY+Amax)^2-(4*convex*IaQY*Amax)))/(2*convex)
  mod.Photo <- (((IaQY+Amax)-sqrt((IaQY+Amax)^2-(4*convex*IaQY*Amax)))/(2*convex))-Rd[i]
  An <- mod.Photo
  #An <- mod.Photo-Rd[i]
  RMSE <- sqrt(mean((Photo-An)^2))  ## RMSE cost function
  return(RMSE)
} 
#==================================================================================================#


####################################################################################################
### EOF.  End of R script file.              
####################################################################################################