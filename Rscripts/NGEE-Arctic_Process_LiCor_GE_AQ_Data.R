####################################################################################################
#          
#   Farquhar Model Fitting Code: A-Q curve analysis
#   Author: Shawn P. Serbin
#
#
#   Estimate aQY & Amax (light and CO2 sat) using Farquhar equation optimization
#
#       --- Uses DEoptim global parameter optimization
#       --- Outputs Rd, aQY, and Amax (when a full curve is availible) estimated parameters 
#               in a single output .csv file with assocaited information
#       --- Output A-Q fit diagnostic figures and DEoptim parameter trace plots
#
#
#   Requirements:
#   DEoptim package
#   install.packages('DEoptim')
#
#
#  	--- Last updated:  04.17.2017 By Shawn P. Serbin <sserbin@bnl.gov>
####################################################################################################


#---------------- Close all devices and delete all variables. -------------------------------------#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
#--------------------------------------------------------------------------------------------------#
