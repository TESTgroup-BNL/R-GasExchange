#--------------------------------------------------------------------------------------------------#
# A set of helper functions
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
##'
##' Read settings file for GE processing
##' 
##' @name ge.settings
##' @title parse settings file
##' @param input.file settings file containing information needed for GE processing
##'
##' @examples
##' \dontrun{
##' settings <- settings()
##' settings <- settings('/home/$USER/settings.xml')
##' }
##'
##' @export
##'
##' @author Shawn P. Serbin
##'
ge.settings <- function(input.file=NULL){
  settings.xml <- NULL
  
  ### Parse input settings file
  if (!is.null(input.file) && file.exists(input.file)) {
    settings.xml <- xmlParse(input.file)  
    # convert the xml to a list
    settings.list <- xmlToList(settings.xml)
    
  } else {
    print("***** WARNING: no settings file defined *****")
  }
  
  # make sure something was loaded
  if (is.null(settings.xml)) {
    #log.error("Did not find any settings file to load.")
    stop("Did not find any settings file to load.")
  }
  
  ### Remove comment or NULL fields
  settings.list <- settings.list[settings.list !="NULL" ]
  
  # Return settings file as a list
  #invisible(settings.list) # invisible option
  return(settings.list)
  
} ### End of function
#==================================================================================================#


#--------------------------------------------------------------------------------------------------#
##'
##' Arrhenius temperature scaling function
##' 
##' @name arrhenius.scaling
##' @title An Arrhenius temperature scaling function for temperature dependent parameters
##' @param temp1 the temperature (in degrees C) of the observation at the base temperature
##' @param temp2 the desired temperature (in degrees C) to scale the observation
##' @param E activation energy of the Arrhenius equation
##' @param Val.temp1 the value of the observation at the base temperature
##' 
##' @examples
##' arrhenius.scaling(temp1=15,temp2=25,E=54.08,Val.temp1=40)
##' 
##' @export
##'
##' @author Shawn P. Serbin
##' 
arrhenius.scaling <- function(temp1=NULL,temp2=NULL,E=NULL,Val.temp1=NULL){
  
  R <- 0.008314472        ## Ideal gas constant
  temp1 <- temp1+273.15   ## First temperature
  temp2 <- temp2+273.15   ## Second temperature
  
  # Arrhenius scaling function
  result <- Val.temp1*exp((E/R)*((1/temp1)-(1/temp2)))
  return(result)
  
} # End of function
#==================================================================================================#


####################################################################################################
### EOF.  End of R script file.              
####################################################################################################