#--------------------------------------------------------------------------------------------------#
# A set of functions for processing A-Ci curve data
#--------------------------------------------------------------------------------------------------#

##' @name ACi.photo.fit
##' @title Fit A-Ci data to C3 model
##' 
##' @param in.dir Input directory or data file. Input data file should match template format
##' @param dataset Input dataset following template format
##' @param out.dir Main output directory for processed curves
##' @param Vcmax.cutoff ppm  Options: E.g. 180, 200, 220, 250. Default 200
##' @param Jmax.cutoff ppm  Options: E.g. 400, 450, 500. Default 450
##' @param Cond.cutoff Drop measurements that have a leaf conductance to water vapour value
##' less than this value.  Default 0.09
##' @param Tleaf.cutoff Acceptable variation in leaf temperature for each group of measurements.
##' Drop any individual observation that has a Tleaf that varies by more than this amount. 
##' Default 0.5
##' @param Ci.checks Realistic range of Ci.  Useful for throwing out erroneous measurements.
##' Default c(0,2000)
##' @param max.iters <- 1000                 ## Max iterations
##' @param RMSE.min <- 0.05                  ## Min value of RMSE to be reached (VTR) during 
##' @param NP Number of population members within the DEoptim minimization. For many problems it is best 
##' to set NP to be at least 10 times the length of the parameter vector.  Default 100
##' @param DWF Differential weighting factor from interval [0,2]. Default to 0.8.
##' @param CR Crossover probability from interval [0,1]. Default to 0.6.
##' @param oxygen Value for ambient O2 concentration (ppm).  Default 210
##' 
##' @note # DEoptim optimization. DEoptim will stop when either 
##' 1) RMSE=RMSE.min or 2) max iters are reached.
##' Should be sufficiently small, e.g. 0.5. Set to 
##' a low number to ensure convergance (e.g. 0.0001)
##' or if unsure. This option can speed up convergence
##' More information see ?DEoptim
##' 
##' #@examples - add example A-Ci dataset.  Use for example
##' 
##' @import DEoptim
##' @export
##' 
##' @references Bernacchi CJ, Bagley JE, Serbin SP, Ruiz-Vera UM, Rosenthal DM, Vanloocke A (2013) Modelling C3 photosynthesis from the chloroplast to the ecosystem. Plant, Cell & Environment, 36, 1641-1657.
##' 
##' 
##' @author Shawn P. Serbin

#ACi.photo.fit <- function(data=NULL,settings.file=NULL,ACi.cutoff=250,response=c("Assim","StomCond"),
#                       op.level=c(1,2,3),op.type=c("single","multiple"),op.method=c("DEoptim"),
#                       lower.bound=c(-3,0,0),upper.bound=c(15,600,600), DEoptim.control=NULL,...){

ACi.photo.fit <- function(in.dir=NULL, dataset=NULL, out.dir=NULL, op.method=c("DEoptim"),
                          Vcmax.cutoff=200,Jmax.cutoff=450,Cond.cutoff=0.09, Tleaf.cutoff=0.5,
                          Ci.checks=c(0,2000),lower.bound=c(-3,0,0),upper.bound=c(15,600,600),
                          max.iters=1000,RMSE.min=0.05,NP=100,DWF=0.8,CR=0.6,
                          oxygen=210) {
  
  dlm = .Platform$file.sep # <--- What is the platform specific delimiter?
  
  # Check for input directory and datafile
  
  
  # Check output directory
  unlink(out.dir,recursive=T) # delete old output if rerunning
  if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
  
  #data(dataSpec_p4) # load requested MM constants
  #mm.constants <- list(R=R,Oxygen=Oxygen,Kc25=Kc25,Ekc=Ekc,Ko25=Ko25,Eko=Eko,Gstar25=Gstar25,
  #                     EGstar=EGstar)
  
  
  ##--- SECTION TEMPORARY.  NEEDS TO BE REMOVED.  PLACE MM OPTIONS IN OTHER SECTION
  # -------------------------------- Temperature adjustment settings ---------------------------------
  # Taken from Bernacchi et al. 2001. Plant, Cell and Environment. 24, 253-259
  #
  # & Long, S. P. and C. J. Bernacchi. 2003. Gas exchange measurements, what can 
  # they tell us about the underlying limitations to photosynthesis? Procedures 
  # and sources of error. Journal of Experimental Botany 54:2393-2401.
  #
  # & Medlyn et al., 2002. Plant, Cell and Environment. 25, 1167-1179
  #
  # & Bernacchi et al., 2013. Plant, Cell and Environment.
  #
  #
  #R       <- 0.008314472               ## Ideal gas constant
  R       <- 8.314                     ## Ideal gas constant. J mol-1 K-1
  Oxygen  <- 210                       ## Oxygen value (ubar)
  #Oxygen  <- 21
  Kc25    <- 404.9                     ## umol m-1
  Ekc     <- 79.430                    ## kJ mol-1
  #Ekc     <- 79430
  Ko25    <- 278.4                     ## mmol m-1
  Eko     <- 36.380                    ## kJ mol-1
  #Eko     <- 36380
  Gstar25 <- 42.75                     ## umol m-1
  EGstar  <- 37.830                    ## kJ mol-1
  #EGstar  <- 37830
  mm.constants <- list(R=R,Oxygen=Oxygen,Kc25=Kc25,Ekc=Ekc,Ko25=Ko25,Eko=Eko,Gstar25=Gstar25,
                       EGstar=EGstar)
  rm(R,Oxygen,Kc25,Ekc,Ko25,Eko,Gstar25,EGstar)
  #--------------------------------------------------------------------------------------------------#
  
  
}
#==================================================================================================#


####################################################################################################
### EOF.  End of R script file.              
####################################################################################################