####################################################################################################
#          
#   Farquhar Model Fitting Code: A-Ci curve analysis
#   Author: Shawn P. Serbin
#
#
#   Info:
#   Estimate Vcmax and Jmax using Farquhar equation optimization
#   --- Uses DEoptim global parameter optimization
#   --- Outputs Rd, Vcmax, and Jmax (when a full curve is availible) estimated parameters 
#       in a single output .csv file with assocaited information
#   --- Output A-Ci fit diagnostic figures and DEoptim parameter trace plots
#
#
#   Requirements:
#   DEoptim package
#   install.packages('DEoptim')
#
#
#  	--- Last updated:  03.30.2016 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#---------------- *User defined settings.* --------------------------------------------------------#
### Location of R scripts.  Needed for Farquhar model optimization. Contains functions.
r.functions <- '/Volumes/TEST/Projects/Leaf_Photosynthesis/NGEE-Arctic/'

### Input LI6400 dataset.  First define location of file (i.e. directory). 
#in.dir <- '/Volumes/TEST/Projects/Global_Vcmax_Synthesis_Project/'
in.dir <- '/Volumes/TEST/Projects/NGEE-Arctic/Data/Barrow/Gas_Exchange/ACi/'

### Barrow data
#dataset <- 'Rogers_GVSP_QA_for_400_errors_round2.csv'
dataset <- 'NGEE-Arctic_2012-2015_ACi_AR3.csv'

### Define input file to be processed
ge.data <- read.table(paste(in.dir,"/",dataset,sep=""), header=T,sep=",")
levels(ge.data$Plot)[1] <- "none"
ge.data$Replicate[ge.data$Replicate==-9999] <- 1
#unique(ge.data$Replicate)
ge.data[ge.data==-9999] <- NA
ge.data[ge.data==""] <- NA
summary(ge.data)          ## Summary of dataset
names(ge.data)

# Calculate deltaT
ge.data$deltaT <- ge.data$Tleaf - ge.data$Tair

### Main output directory 
out.dir <- '' # update with output directory

unlink(out.dir,recursive=T) # delete old output if rerunning
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)

# *********************************** QA/QC Options ***********************************
### Vcmax to Jmax Ci cutoff.
Vcmax.cutoff <- 400 # ppm  Options: E.g. 180, 200, 220, 250, 350
Jmax.cutoff  <- 650 # ppm  Options: E.g. 400, 450, 500
Vcmax.min.n <- 2    # minimum number of obs for calculating Vcmax  Best is 3.
Jmax.min.n <- 1     # minimum number of obs for calculating Jmax.  Best is 3.

### Sample QC checks
Cond.cutoff <- 0.05        ## Throw out observations with Cond < cutoff. E.g. 0.08
Ci.cutoff <- c(0,2000)    ## Throw out observations with Ci out of bounds
Tleaf.cutoff <- 1.6       ## How much Tleaf variation to accept in curve data. E.g. 1 
# would allow variation of 1 degree around the mean Tleaf
# *********************************** QA/QC Options ***********************************

# ********************************** DEoptim Options **********************************
### DEoptim options. Controls the parameter search space for each DEoptim iteration.
lower.bound <- c(5,-5,5)     ## Lower bound of Vcmax, Rd, & Jmax
upper.bound <- c(800,15,800)      ## Upper bound of Vcmax, Rd, & Jmax
max.iters <- 1000                 ## Max iterations
RMSE.min <- 0.05                  ## Min value of RMSE to be reached (VTR) during 
NP <- 100                         ## Number of population members. For many problems it is best to set 
## NP to be at least 10 times the length of the parameter vector.
DWF <- 0.8                        ## Differential weighting factor from interval [0,2]. Default to 0.8.
CR <- 0.6                         ## Crossover probability from interval [0,1]. Default to 0.5.
# DEoptim optimization. DEoptim will stop when either 
# 1) RMSE=RMSE.min or 2) max iters are reached.
# Should be sufficiently small, e.g. 0.5. Set to 
# a low number to ensure convergance (e.g. 0.0001)
# or if unsure. This option can speed up convergence
#see ?DEoptim
# ********************************** DEoptim Options **********************************


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
#Ekc     <- 79.430                    ## kJ mol-1
Ekc     <- 79430
Ko25    <- 278.4                     ## mmol m-1
#Eko     <- 36.380                    ## kJ mol-1
Eko     <- 36380
Gstar25 <- 42.75                     ## umol m-1
#EGstar  <- 37.830                    ## kJ mol-1
EGstar  <- 37830
mm.constants <- list(R=R,Oxygen=Oxygen,Kc25=Kc25,Ekc=Ekc,Ko25=Ko25,Eko=Eko,Gstar25=Gstar25,
                     EGstar=EGstar)
rm(R,Oxygen,Kc25,Ekc,Ko25,Eko,Gstar25,EGstar)
#--------------------------------------------------------------------------------------------------#


#==================================================================================================#
# End of user defined options. Start of program.
#==================================================================================================#


###################################### START SCRIPT ################################################
print('**************************************************')
print('**************** STARTING SCRIPT *****************')
print('**************************************************')
print(' ')
####################################################################################################


#---------------- Load required libraries ---------------------------------------------------------#
# Info: Loads required R libraries and warns if package is not availible.
ok = require(DEoptim) ; if (! ok) 
  stop("*** Package DEoptim is not available.  This is needed for model optimization ***")
rm(ok)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Load utils and Farquhar functions for data QA/QC and optimization
source(paste0(r.functions,'/','photo.processing.functions.R'))
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Setup output dirs
dlm = .Platform$file.sep # <--- What is the platform specific delimiter?
aci.fig.dir = paste(out.dir,dlm,"ACi_Diagnostics",sep="")
if (! file.exists(aci.fig.dir)) dir.create(aci.fig.dir)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
# Apply data QA/QC functions
#data <- ge.data  #For debugging
ge.data.qc <- data.qc(data=ge.data,out.dir=out.dir,Cond.cutoff=Cond.cutoff,Ci.cutoff=Ci.cutoff,
                      Tleaf.cutoff=Tleaf.cutoff)
ge.data.qc[is.na(ge.data.qc)] <- NA
out.qc.ge.data <- data.frame(ge.data.qc$Sample.Info,ge.data.qc$GE.data)
write.csv(out.qc.ge.data,file=paste(out.dir,"/","QC_GE_Data.csv",sep=""),row.names=FALSE)
rm(out.qc.ge.data,ge.data)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Apply temperature adjustments to MM constants
loc <- match("TLEAF",toupper(names(ge.data.qc$GE.data)))
TleafK  <- ge.data.qc$GE.data[,loc[1]]+273.15                                ## Leaf temperature
ge.data.qc$GE.data$Kc <- mm.constants$Kc25*exp(mm.constants$Ekc*(TleafK-298.15)/
                                                 (298.15*mm.constants$R*TleafK))       ## Kc - based on Medlyn et al
ge.data.qc$GE.data$Ko <- mm.constants$Ko25*exp(mm.constants$Eko*(TleafK-298.15)/
                                                 (298.15*mm.constants$R*TleafK))       ## Ko - based on Medlyn et al
ge.data.qc$GE.data$Gstar <- mm.constants$Gstar25*exp(mm.constants$EGstar*(TleafK-298.15)/
                                                       (298.15*mm.constants$R*TleafK)) ## Gamma star. Medlyn et al
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Get sample info and summary stats
samples <- unique(ge.data.qc$Sample.Info)

### Get data names
data.names <- names(ge.data.qc$GE.data)
remove.nms <- match(c("AREA","PRESS","CO2R","CO2S","PHOTO","CI","KC","KO","GSTAR"),toupper(data.names))
data.names <- data.names[-remove.nms]
data.names # What to summarize

index <- within(ge.data.qc$GE.data, indx <- as.numeric(interaction(ge.data.qc$Sample.Info, 
                                                                   drop=TRUE,lex.order=TRUE)))

samples <- data.frame(samples,Index=unique(index$indx))
samples <- samples[order(samples$Index),]
row.names(samples) <- seq(len=nrow(samples))
samples <- samples[,-match("Index",names(samples))]

### Obs stats
means <- aggregate(.~index$indx,data=index,mean)
sdevs <- aggregate(.~index$indx,data=index,sd)
means <- means[data.names]
sdevs <- sdevs[data.names]
CVs <- (sdevs/means)*100

# Create output dataframe of diagnostic info
names(means) <- paste0("Mean.",names(means))
names(CVs) <- paste0("PercCV.",names(CVs))
out.diagnostics <- data.frame(means[,1:2],PercCV.Tleaf=CVs$PercCV.Tleaf,means[,3:7],
                              PercCV.VpdL=CVs$PercCV.VpdL,means[,8:10],PercCV.Ci_Ca=CVs$PercCV.Ci_Ca,
                              Mean.Cond=means[,11],PercCV.Cond=CVs$PercCV.Cond)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Fit Farquhar model to data
data <- data.frame(ge.data.qc$Sample.Info,ge.data.qc$GE.data)

## Setup outputs

# Fitted params
Rd <- array(data=NA,dim=dim(samples)[1])
Vcmax <- array(data=NA,dim=dim(samples)[1])
Jmax <- array(data=NA,dim=dim(samples)[1])

# Other outputs
Meas.Photo.atCO2R.400 <- array(data=NA,dim=dim(samples)[1])
Meas.CiCa.atCO2R.400 <- array(data=NA,dim=dim(samples)[1])
Mod.Photo.atCi.388 <- array(data=NA,dim=dim(samples)[1])
Mod.Photo.atCi.388.timesCiCa <- array(data=NA,dim=dim(samples)[1]) # default 0.7
Est.Photo.stomata.limit <- array(data=NA,dim=dim(samples)[1])
Est.Ci.colimit <- array(data=NA,dim=dim(samples)[1])

# Fit
RMSE.DEoptim <- array(data=NA,dim=dim(samples)[1])
RMSE.photo <- array(data=NA,dim=dim(samples)[1])

# Main outer loop
system.time(for (i in 1:dim(samples)[1]) {
#system.time(for (i in 1:5) {
  sub.data <- merge(data,samples[i,],by=names(samples[i,]))
  #sub.data = sub.data[order(sub.data$Ci_Pcorr),]
  sub.data <- sub.data[order(sub.data$Ci),]
  print("--- Processing Sample: ")
  print(paste(as.vector(unlist(samples[i,])),collapse=" "))
  
  # Determine what level of GE processing
  chk2 <- length(which(sub.data$Ci<Vcmax.cutoff))
  chk3 <- length(which(sub.data$Ci>Jmax.cutoff))
  
  # Find merasured photo & CiCa @ CO2R 400
  keep <- which(sub.data$CO2R >= 380 & sub.data$CO2R <= 420)
  if (length(keep)>0) {
    Meas.Photo.atCO2R.400[i] <- mean(sub.data$Photo[keep],na.rm=TRUE)
    Meas.CiCa.atCO2R.400[i] <- mean(sub.data$Ci_Ca[keep],na.rm=TRUE)
  } else {
    Meas.Photo.atCO2R.400[i] <- -9999
    Meas.CiCa.atCO2R.400[i] <- 0.7 #default
  }
  
  # Check if there are enough obs. Remove curves with to few Vcmax/Jmax samples
  if (dim(sub.data)[1]<Vcmax.min.n | chk2<Vcmax.min.n & chk3<Jmax.min.n){
    Rd[i] <- -9999
    Vcmax[i] <- -9999
    Jmax[i] <- -9999
    RMSE.DEoptim[i] <- -9999
    RMSE.photo[i] <- -9999
  } else {
    Ci <- sub.data$Ci
    Oxy <- mm.constants$Oxygen
    Kc <- sub.data$Kc
    Ko <- sub.data$Ko
    Gstar <- sub.data$Gstar
    Km <- Kc*(1+(Oxy/Ko))
    Photo <- sub.data[,match("PHOTO",toupper(names(sub.data)))] # Find photo data
    f.model <- NA # Type variable for defining plots
    
    # Vcmax and Jmax
    if (chk2>=Vcmax.min.n & chk3>=Jmax.min.n){
      f.model <- 1
      fit <- DEoptim(ACi.full,lower=lower.bound, upper=upper.bound, DEoptim.control(NP=NP,
                     F=DWF,CR=CR, itermax=max.iters, VTR=RMSE.min, strategy=2, trace=FALSE))
      tempout <- data.frame(fit$optim)
      Vcmax[i] <- tempout[1,1]
      Rd[i] <- tempout[2,1]
      Jmax[i] <- tempout[3,1]
      RMSE.DEoptim[i] <- fit$optim$bestval
      fitted_photo = pmin(((Vcmax[i]*(Ci-Gstar))/(Ci+Km))-Rd[i],
                          ((Jmax[i]*(Ci-Gstar))/((4.5*Ci)+(10.5*Gstar)))-Rd[i])
      residuals = fitted_photo-Photo
      RMSE.photo[i] = sqrt(mean((residuals)^2))
      
      ### Display info to console 
      print(paste("Vcmax: ",round(Vcmax[i],2)))
      print(paste("Jmax: ",round(Jmax[i],2)))
      print(paste("Rd: ",round(Rd[i],4)))  
      print(paste("DEoptim RMSE: ",round(RMSE.DEoptim[i],2)))
      print(paste("Photo RMSE: ",round(RMSE.photo[i],2)))
      print("-----------------------------------------------------------")
      flush.console()
      
      # Vcmax only
    } else if (chk2>=Vcmax.min.n & chk3<Jmax.min.n) {
      f.model <- 2
      fit <- DEoptim(ACi.rubisco,lower=lower.bound[1:2], upper=upper.bound[1:2], DEoptim.control(NP=NP,
                    F=DWF,CR=CR, itermax=max.iters, VTR=RMSE.min, strategy=2, trace=FALSE))
      tempout <- data.frame(fit$optim)
      Vcmax[i] <- tempout[1,1]
      Rd[i] <- tempout[2,1]
      Jmax[i] <- -9999
      RMSE.DEoptim[i] <- fit$optim$bestval
      fitted_photo = ((Vcmax[i]*(Ci-Gstar))/(Ci+Km))-Rd[i]
      residuals = fitted_photo-Photo
      RMSE.photo[i] = sqrt(mean((residuals)^2))
      
      ### Display info to console 
      print(paste("Vcmax: ",round(Vcmax[i],2)))
      print(paste("Jmax: ",round(Jmax[i],2)))
      print(paste("Rd: ",round(Rd[i],4)))  
      print(paste("DEoptim RMSE: ",round(RMSE.DEoptim[i],2)))
      print(paste("Photo RMSE: ",round(RMSE.photo[i],2)))
      print("-----------------------------------------------------------")
      flush.console()
      
      # Jmax only
    } else if (chk2<Vcmax.min.n & chk3>=Jmax.min.n){
      f.model <- 3
      fit <- DEoptim(ACi.rubp,lower=lower.bound[3:2], upper=upper.bound[3:2], DEoptim.control(NP=NP,
                     F=DWF,CR=CR, itermax=max.iters, VTR=RMSE.min, strategy=2, trace=FALSE))
      tempout <- data.frame(fit$optim)
      Vcmax[i] <- -9999
      Jmax[i] <- tempout[1,1]
      Rd[i] <- tempout[2,1]
      RMSE.DEoptim[i] <- fit$optim$bestval
      fitted_photo = ((Jmax[i]*(Ci-Gstar))/((4.5*Ci)+(10.5*Gstar)))-Rd[i]
      residuals = fitted_photo-Photo
      RMSE.photo[i] = sqrt(mean((residuals)^2))
      
      ### Display info to console 
      print(paste("Vcmax: ",round(Vcmax[i],2)))
      print(paste("Jmax: ",round(Jmax[i],2)))
      print(paste("Rd: ",round(Rd[i],4)))  
      print(paste("DEoptim RMSE: ",round(RMSE.DEoptim[i],2)))
      print(paste("Photo RMSE: ",round(RMSE.photo[i],2)))
      print("-----------------------------------------------------------")
      flush.console()
      
    } # End processing if/else
    
    # Calculate photo
    if (f.model==1 | f.model==2) {
      ci.loc <- 388*Meas.CiCa.atCO2R.400[i]
      Mod.Photo.atCi.388[i] <- ((Vcmax[i]*(388-mean(Gstar)))/(388+mean(Km)))-Rd[i]
      Mod.Photo.atCi.388.timesCiCa[i] <- ((Vcmax[i]*(ci.loc-mean(Gstar)))/(ci.loc+mean(Km)))-Rd[i]
      Est.Photo.stomata.limit[i] <- Mod.Photo.atCi.388[i] - Mod.Photo.atCi.388.timesCiCa[i]
    } else {
      Mod.Photo.atCi.388[i] <- -9999
      Mod.Photo.atCi.388.timesCiCa[i] <- -9999
      Est.Photo.stomata.limit[i] <- -9999
    }
    rm(ci.loc)
    
    # Calculate co-limitation from fitted data
    if (f.model==1) {
      modCi <- seq(-2,2000,1)
      fit.vcmax <- ((Vcmax[i]*(modCi-mean(Gstar)))/(modCi+mean(Km)))-Rd[i]
      fit.jmax <- ((Jmax[i]*(modCi-mean(Gstar)))/((4.5*modCi)+(10.5*mean(Gstar))))-Rd[i]
      intersect <- which.min(abs(fit.jmax-fit.vcmax))
      Est.Ci.colimit[i] <- modCi[intersect]
      # debugging
      #plot(modCi,fit.vcmax,type="l",col="blue")
      #lines(modCi,fit.jmax,col="red")
    } else {
      Est.Ci.colimit[i] <- -9999
    }

    ### Generate fit diagnostics
    #if (length(which(sub.data$Ci_Pcorr<Vcmax.cutoff))>=3){
    if (chk2 >=Vcmax.min.n | chk3 >=Jmax.min.n){
      sample.name <- paste(unlist(samples[i,]),collapse="_")
      sample.name <- gsub(pattern="/","-",sample.name)
      params <- data.frame(Vcmax=Vcmax[i],Jmax=Jmax[i],Rd=Rd[i],RMSE.photo=RMSE.photo[i])
      plot.ge.fit(type="A-Ci",sub.data,fit,params,aci.fig.dir,sample.name,f.model)
      rm(f.model)
    }
    rm(fit)
  } # End total observations if/else check
  
  rm(chk2,chk3)
}) #End of outer main loop
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Generate output
output.data <- data.frame(samples,out.diagnostics,Meas.Photo.atCO2R.400,Meas.CiCa.atCO2R.400,
                          Mod.Photo.atCi.388,Mod.Photo.atCi.388.timesCiCa,Est.Photo.stomata.limit,
                          Est.Ci.colimit,Vcmax=Vcmax,Jmax=Jmax,Rd=Rd,RMSE.DEoptim=RMSE.DEoptim,
                          RMSE.photo=RMSE.photo)
output.dataset <- paste(strsplit(dataset,".csv"),".processed.csv",sep="")
write.csv(output.data,file=paste(out.dir,"/",output.dataset,sep=""),row.names=FALSE)

rm(list=ls(all=TRUE))   # clear workspace
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF