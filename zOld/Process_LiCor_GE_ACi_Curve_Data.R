####################################################################################################
#
#            Estimate Vcmax and Jmax using Farquhar equation optimization
#
#           --- Uses DEoptim global parameter optimization
#		        --- Outputs Rd, Vcmax, and Jmax (when a full curve is availible) estimated parameters 
#               in a single output .csv file with assocaited information
#           --- Output A-Ci fit diagnostic figures and DEoptim parameter trace plots
#
#		--- Last updated:  2.05.2013 BY SPS
####################################################################################################

#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#---------------- *User defined settings.* --------------------------------------------------------#
### Input LI6400 dataset.  First define location of file (i.e. directory). 
in.dir <- ('/Users/serbin/Data/Dropbox/Soybean_Aphid_Project/Data/LI6400_Data/')

### Define input file to be processed
ge.data <- read.table(paste(in.dir,"/",'Compiled_2012.csv',sep=""), header=T,sep=",")
summary(ge.data)          ## Summary of dataset

### Main output directory 
out.dir <- ('/Users/serbin/Data/Dropbox/Soybean_Aphid_Project/Data/Processed_LI6400_Data/Compiled_2012/')
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)

### Vcmax to Jmax Ci cutoff.
ACi.cutoff <- 220 #ppm

### Sample QC checks
Cond.cutoff <- 0.08       ## Throw out sample sets with Cond < cutoff
Ci.cutoff <- c(0,2000)    ## Throw out sample sets with Ci out of bounds

### DEoptim options. Controls the parameter search space for each DEoptim iteration.
lower.bound(-10,5,5)      ## Lower bound of Rd, Vcmax, & Jmax
upper.bound(0,500,600)    ## Upper bound of Rd, Vcmax, & Jmax

# -------------------------------- Temperature adjustment settings ---------------------------------
# Taken from Bernacchi et al. 2001. Plant, Cell and Environment. 24, 253-259
#
# & Long, S. P. and C. J. Bernacchi. 2003. Gas exchange measurements, what can 
# they tell us about the underlying limitations to photosynthesis? Procedures 
# and sources of error. Journal of Experimental Botany 54:2393-2401.
#
# & Medlyn et al., 2002. Plant, Cell and Environment. 25, 1167-1179

R       <- 0.008314472               ## Ideal gas constant
Oxygen  <- 210                       ## Oxygen value (ubar)
Kc25    <- 404.9                     ## umol m-1+
Ekc     <- 79.430                    ## kJ mol-1
Ko25    <- 278.4                     ## mmol m-1
Eko     <- 36.380                    ## kJ mol-1
Gstar25 <- 42.75                     ## umol m-1
EGstar  <- 37.830                    ## kJ mol-1
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
### Setup output dirs
dlm = .Platform$file.sep # <--- What is the platform specific delimiter?
aci.fig.dir = paste(out.dir,dlm,"ACi_Diagnostics",sep="")
if (! file.exists(aci.fig.dir)) dir.create(aci.fig.dir)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Define Farquhar models for optimization

# Rd and Vcmax only
A.1 = function(x) {
  ## Setup variables to solve for
  Rd    = x[1]  ## Resp param
  Vcmax = x[2]  ## Vcmax param
  Rubisco = ifelse(Ci<ACi.cutoff,((Vcmax*(Ci-Gstar))/(Ci + Km))-Rd,9999)
  assim = Rubisco    
  RMSE = sqrt(mean((Photo-assim)^2))  ## RMSE cost function
}

# Rd, Vcmax, and Jmax
A.2 = function(x) {
  ## Setup variables to solve for
  Rd    = x[1]  ## Resp param
  Vcmax = x[2]  ## Vcmax param
  Jmax  = x[3]  ## Jmax param
  Rubisco = ifelse(Ci<ACi.cutoff,((Vcmax*(Ci-Gstar))/(Ci + Km))-Rd,9999)
  RuBP = ifelse(Ci>ACi.cutoff,(Jmax*(Ci-Gstar)/((4.5*Ci)+(10.5*Gstar)))-Rd,9999)
  assim = pmin(Rubisco,RuBP,inter)-Rd  
  RMSE = sqrt(mean((Photo-assim)^2))  ## RMSE cost function
}

#--------------------------------------------------------------------------------------------------#


#---------------- Setup data for optimization -----------------------------------------------------#
### Remove samples not passing initial QC
loc <- match("QC",toupper(names(ge.data)))
remove <- which(ge.data[loc]==1)
if(length(remove)>0){
  ge.data <- ge.data[,-remove]
}
rm(loc,remove)

### Find sample info and data columns
pattern <- c("QC","COMMENTS")
x <- toupper(names(ge.data))
remove <- match(pattern,x)
if (length(remove)>0){
  ge.data <- ge.data[,-remove]
}
rm(pattern,x,remove)

### Find data columns
pattern <- c("Tair","Tleaf","deltaT","RH_R","RH_S","PRESS","PARi","CO2Ref","PHOTO","COND","Ci")
pattern <- toupper(pattern)
x <- toupper(names(ge.data))
keep <- match(pattern,x)

sample.info <- ge.data[,-keep]
ge.data <- ge.data[,keep]
rm(pattern,x,keep)

### Apply some quick QC filters
dims <- dim(ge.data)
loc <- match(c("COND","CI"),toupper(names(ge.data)))
check1 <- which(ge.data[,loc[1]]<0.08)
check2 <- which(ge.data[,loc[2]]<Ci.cutoff[1] | ge.data[,loc[2]]>Ci.cutoff[2])
check3 <- which(ge.data[,loc[2]]<Ci.cutoff[1] | ge.data[,loc[2]]>Ci.cutoff[2] | ge.data[,loc[1]]<0.08)
#nms <- unique(sample.info[check1,])
nms1 <- sample.info[check1,]
vals1 <- ge.data[check1,loc[1]]
nms2 <- sample.info[check2,]
vals2 <- ge.data[check2,loc[2]]

temp1 <- data.frame(nms1,Cond=vals1)
temp2 <- data.frame(nms2,Ci=vals2)
if (dim(temp1)[1]>0){
  write.csv(temp1,file=paste(out.dir,"/","Failed_QC_Conductance_Check.csv",sep=""),row.names=FALSE)
}
if (dim(temp2)[1]>0){
write.csv(temp2,file=paste(out.dir,"/","Failed_QC_Ci_Cutoff_Check.csv",sep=""),row.names=FALSE)
}

### Remove bad data
if (length(check3>0)){
  ge.data <- ge.data[-check3,]
  sample.info <- sample.info[-check3,]
}
row.names(ge.data) <- seq(len=nrow(ge.data))
row.names(sample.info) <- seq(len=nrow(sample.info))
rm(dims,loc,check1,check2,check3,nms1,vals1,nms2,vals2,temp1,temp2,Ci.cutoff,Cond.cutoff)

### Apply temperature adjustments to MM constants
loc <- match(c("CI","PRESS"),toupper(names(ge.data)))
ge.data$Ci_Pcorr = ge.data[,loc[1]]*ge.data[,loc[2]]*0.01                   ## Pressure Adj. of Ci
loc <- match("TLEAF",toupper(names(ge.data)))
TleafK  <- ge.data[,loc[1]]+273.15                                          ## Leaf temperature
ge.data$Kc = mm.constants$Kc25*exp(mm.constants$Ekc*(TleafK-298.15)/
                                     (298.15*mm.constants$R*TleafK))        ## Kc - based on Medlyn et al
ge.data$Ko = mm.constants$Ko25*exp(mm.constants$Eko*(TleafK-298.15)/
                                       (298.15*mm.constants$R*TleafK))      ## Ko - based on Medlyn et al
ge.data$Gstar = mm.constants$Gstar25*exp(mm.constants$EGstar*(TleafK-298.15)/
                                           (298.15*mm.constants$R*TleafK)) ## Gamma star. Medlyn et al

#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Fit farquhar model to data.
samples <- unique(sample.info)

which(sample.info==samples[1,])

sub = subset(data_6400_2, SAMPLE==samples[i])

#--------------------------------------------------------------------------------------------------#









