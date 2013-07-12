#--------------------------------------------------------------------------------------------------#
##'
##' A function to read in formatted LiCor 6400 data for processing
##'
##' @name read.ge.data
##' @title A function to read in formatted LiCor 6400 data files
##' 
##' @export
##' 
read.ge.data <- function(file.dir=NULL,out.dir=NULL,ge.file.ext=".csv",
                         output.file.ext=".csv",ge.dataframe=FALSE,settings.file=NULL){
  
  ### Set platform specific file path delimiter.  Probably will always be "/"
  dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?

  ### Check for proper input
  if (is.null(settings.file) && is.null(file.dir)){
    print("*********************************************************************************")
    stop("******* ERROR: No input file directory given in settings file or function call. *******")
  } else if (!is.null(file.dir)){
    file.dir <- file.dir
  } else if (!is.null(settings.file$asd.dir)){
    file.dir <- settings.file$asd.dir
  } 
  
  
  ### Algorithm options
  deoptim.lower.bounds <- as.numeric(strsplit(settings$algorithm.options$deoptim.lower.bounds,split = ",")[[1]])
  
} # End of function call
#==================================================================================================#


####################################################################################################
### EOF.  End of R script file.              
####################################################################################################