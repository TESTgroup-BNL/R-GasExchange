####################################################################################################
#          
#   Fit temperature response to Vcmax/Jmax data
#   Author: Shawn P. Serbin
#
#   Requirements:
#   DEoptim package
#   install.packages('DEoptim')
#
#  	--- Last updated:  03.30.2016 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#


#---------------- *User defined settings.* --------------------------------------------------------#
## Sample options
species <- ''     # define the species to fit or NULL for all species/samples
t.response <- FALSE  # TRUE just use data collected to estimate Ea, or FALSE to use all data
date <- NULL
year <- NULL

### Location of R scripts.  Needed for Farquhar model optimization. Contains functions.
r.functions <- ''

### Input LI6400 dataset.  First define location of file (i.e. directory). 
in.dir <- ''

### Barrow data
dataset <- ''
ge.data <- read.table(paste(in.dir,"/",dataset,sep=""), header=T,sep=",")
summary(ge.data)          ## Summary of dataset

### Main output directory 
out.dir <- ''
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)

# ********************************** Temp model options **********************************
model.type <- "Arrhenius"         ## Options: Arrhenius, Peaked, ?June?
stand.temp <- 10 #25              ## Options: 15,25, etc
Photo.Parameter <- "Jmax"        ## Options: Vcmax, Jmax
# ********************************** Temp model options **********************************

# ********************************** DEoptim Options **********************************
### DEoptim options. Controls the parameter search space for each DEoptim iteration.

## Arrehnius
lower.bound <- c(0,0)             ## Lower bounds
upper.bound <- c(500,500)         ## Upper bounds
max.iters <- 200                 ## Max iterations. E.g. 300, 500, 1000
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

## Peaked model

# ********************************** DEoptim Options **********************************
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
fig.dir = paste(out.dir,dlm,"Temp_Response_Diagnostics",sep="")
if (! file.exists(fig.dir)) dir.create(fig.dir)
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Extract sample info
ind1 <- which(toupper(names(ge.data))=="DATE")
sample.info <- ge.data[,seq(1,ind1,1)]  # keep everything to the left of Date as sample info
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Fit temp functions

##### Arrehnius
# Grab the relevant data
if (!is.null(species)){
  sub.data <- data.frame(ge.data[which(sample.info$USDA_Species_Code==species),])
  file <- unique(sub.data$USDA_Species_Code)
} else {
  sub.data <- data.frame(sample.info,ge.data)
  file <- "All_data"
}

if (t.response==TRUE) {
  sub.data <- sub.data[duplicated(sub.data$Sample_Barcode),]
}

# Remove missing
remove <- which(sub.data$Jmax==-9999)
sub.data <- droplevels(sub.data[-remove,])

Tleaf <- sub.data$Mean.Tleaf
fit.data <- sub.data[Photo.Parameter]

# sanity check
plot(Tleaf,unlist(fit.data))

# Fit
R <- 0.008314472 ## Ideal gas constant
fit <- DEoptim(arrhenius,lower=lower.bound, upper=upper.bound, DEoptim.control(NP=NP,
               F=DWF,CR=CR, itermax=max.iters, VTR=RMSE.min, strategy=2, trace=TRUE))
tempout <- data.frame(fit$optim)
Standardized.param <- tempout[1,1]
E <- tempout[2,1]
RMSE.DEoptim <- fit$optim$bestval

fitted <- Standardized.param*exp((E*((Tleaf+273.15)-(stand.temp+273.15)))/((stand.temp+273.15)*R*(Tleaf+273.15)))
residuals <- fitted-fit.data
RMSE <- sqrt(mean((residuals)^2))

params <- data.frame(Param=Standardized.param,E=E,RMSE=RMSE)
plot.temp.response(type="Arrhenius",var=Photo.Parameter,sub.data,fit,params,fig.dir,file=file)

#type="Arrhenius"
#data <- sub.data
#DEoptim.output <- fit
#var <- Photo.Parameter
#outdir <- fig.dir
#file <- unique(sub.data$USDA_Species_Code)

# Peaked model
# ** NOT YET IMPLEMENTED **
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### EOF
