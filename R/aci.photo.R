#--------------------------------------------------------------------------------------------------#
# A set of functions for processing A-Ci curve data
#--------------------------------------------------------------------------------------------------#

##' @name ACi.photo.fit
##' @title Fit A-Ci data to C3 model
##' 
##' @param in.dir Input directory or data file. Input data file should match template format
##' @param out.dir Main output directory for processed curves
##' @param Vcmax.cutoff ppm  Options: E.g. 180, 200, 220, 250
##' @param Jmax.cutoff ppm  Options: E.g. 400, 450, 500
##' @param oxygen Value for ambient O2 concentration (ppm)
##' @param
##' 
##' @references

#ACi.photo.fit <- function(data=NULL,settings.file=NULL,ACi.cutoff=250,response=c("Assim","StomCond"),
#                       op.level=c(1,2,3),op.type=c("single","multiple"),op.method=c("DEoptim"),
#                       lower.bound=c(-3,0,0),upper.bound=c(15,600,600), DEoptim.control=NULL,...){

ACi.photo.fit <- function(in.dir=NULL, out.dir=NULL, op.method=c("DEoptim"),lower.bound=c(-3,0,0),
                          upper.bound=c(15,600,600),
                          oxygen = 210) {
  
  #data(dataSpec_p4) # load requested MM constants
  #mm.constants <- list(R=R,Oxygen=Oxygen,Kc25=Kc25,Ekc=Ekc,Ko25=Ko25,Eko=Eko,Gstar25=Gstar25,
  #                     EGstar=EGstar)
}
#==================================================================================================#


####################################################################################################
### EOF.  End of R script file.              
####################################################################################################