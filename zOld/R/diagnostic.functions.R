#--------------------------------------------------------------------------------------------------#
# A set of helper functions for LiCor GE data processing
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
##'
##' A function plot raw A-Ci, A-Q or other analysis curves.  Shows flagged observations.
##'
##' @name plot.raw 
##' @title A function plot raw A-Ci, A-Q or other analysis curves.  Shows flagged observations.
##' @param file.dir 
##' @param data.file Input data file to process
##' @param out.dir Main output directory for diagnostic plots
##' @param data use data already in the environment
##' @param type The type of gas exchange data to plot.  Current options: A-Ci and A-Q
##' @param output.type .png or .pdf
##' 
##' @export
##' 
##' @author Shawn P. Serbin
##'
plot.raw <- function(file.dir=NULL, data.file=NULL, out.dir=NULL, data=NULL, type=c("A-Ci","A-Q")) {
  #!! NEEDS WORK !!
  
  type <- match.arg(type)
  sep <- .Platform$file.sep
  #data <- data.file # temporary.  replace with updated code
  
  ### Plot params - Make these user selectable?
  cexaxis <- 1.2
  cexlab <- 1.4
  
  if (!is.null(data)){
    ge.data <- data
  }
  
  # TODO: Need to setup loop to run over all group measurements in file
  
  if (type=="A-Ci"){
    loc1 <- match(c("CI","PHOTO"),toupper(names(data)))
    
    pdf(paste(outdir,sep,file,".pdf",sep=""),height=8,width=10)
    plot(data[,loc1[1]],data[,loc1[2]],pch=21,bg="grey70",cex=3,cex.axis=cexaxis,xlim=c(0,range(data[,loc1[1]])[2]),
         ylim=c(0,range(data[,loc1[2]])[2]),cex.lab=cexlab,xlab="Ci",ylab="Photo",main=paste(file))
    box(lwd=2.2) # update to show QC=1 points as a different color/symbol for initial assessment.
    # enable user input to select bad points?
    
    # Eg. # Interacting with a scatterplot
    #attach(mydata)
    #plot(x, y) # scatterplot
    #identify(x, y, labels=row.names(mydata)) # identify points
    #coords <- locator(type="l") # add lines
    #coords # display list 
    
    dev.off() # output the figure
    
  } else if (type=="A-Q"){
    
    loc1 <- match(c("PARI","PHOTO"),toupper(names(data)))
    pdf(paste(outdir,sep,file,".pdf",sep=""),height=8,width=10)
    plot(data[,loc1[1]],data[,loc1[2]],pch=21,bg="grey70",cex=3,cex.axis=cexaxis,xlim=c(0,range(data[,loc1[1]])[2]),
         ylim=c(0,range(data[,loc1[2]])[2]+1),cex.lab=cexlab,xlab="PARi",ylab="Photo",main=paste(file))
    box(lwd=2.2)
    
    
  } # End if/else/then
  
  
} # End of function call
#==================================================================================================#


####################################################################################################
### EOF.  End of R script file.              
####################################################################################################