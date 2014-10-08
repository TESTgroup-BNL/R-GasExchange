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
##' @param data 
##' @param type The type of gas exchange data to plot.  Current options: A-Ci and A-Q
##' @param output.type .png or .pdf
##' 
##' @author Shawn P. Serbin
##' 
plot.raw <- function(file.dir=NULL, data.file=NULL, out.dir=NULL, data, type=c("A-Ci","A-Q")) {
  type <- match.arg(type)
  sep <- .Platform$file.sep
  data <- data.file # temporary.  replace with updated code
  
  ### Plot params - Make these user selectable?
  cexaxis <- 1.2
  cexlab <- 1.4
  
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


#--------------------------------------------------------------------------------------------------#
##'
##' A work in progress !! NEEDS UPDATING !!
##'
##' GE Diagnostic plots 
##' @param type which type of data to invert. A-Ci or A-Q
##' @param data A-Ci data used for model fit
##' @param DEoptim.output DEoptim output structure
##' @param params 
##' @param outdir output directory for A-Ci diagnostic figures
##' @param file output filename
##' @param f.model full, vcmax, or jmax
##' 
##' @author Shawn P. Serbin
##' 
plot.ge.fit <- function(type=c("A-Ci","A-Q"),data,DEoptim.output,params,outdir,file,f.model){
  type <- match.arg(type)
  sep <- .Platform$file.sep
  param.length <- length(DEoptim.output$optim$bestmem)
  
  ### Plot params
  cexaxis <- 1.2
  cexlab <- 1.4
  
  if (type=="A-Ci"){
    loc1 <- match(c("CI","PHOTO"),toupper(names(data)))
    
    # Figures
    pdf(paste(outdir,sep,file,".pdf",sep=""),height=8,width=10)
    plot(data[,loc1[1]],data[,loc1[2]],pch=21,bg="grey70",cex=3,cex.axis=cexaxis,xlim=c(0,range(data[,loc1[1]])[2]),
         ylim=c(0,range(data[,loc1[2]])[2]),cex.lab=cexlab,xlab="Ci",ylab="Photo",main=paste(file))
    box(lwd=2.2)
    
    Oxygen <- mm.constants$Oxygen
    loc2 <- match(c("KC","KO","GSTAR","TLEAF"),toupper(names(data)))
    Kc <- mean(data[,loc2[1]])
    Ko <- mean(data[,loc2[2]])
    Gstar <- mean(data[,loc2[3]])
    
    plot.x <- (data[,loc1[1]]-Gstar)/(data[,loc1[1]]+(Kc*(1+Oxygen/Ko)))
    plot(plot.x,data[,loc1[2]],pch=21,bg="grey70",cex=3,cex.axis=cexaxis,xlim=c(0,range(plot.x)[2]),
         ylim=c(0,range(data[,loc1[2]])[2]),cex.lab=cexlab,xlab="Ci-Gstar/Ci+Km",ylab="Photo",main=paste(file))
    box(lwd=2.2)
    
    # DEoptim trace plots
    if (f.model==1){
      par(mfrow=c(3,1),mar=c(4,4.1,1,2)) #b, l, t, r
      plot(DEoptim.output$member$bestmemit[,1],pch=21,bg="dark grey",col="dark grey",
           cex=1.2,xlab="Iteration",ylab="Vcmax",cex.axis=cexaxis,cex.lab=cexlab)
      lines(DEoptim.output$member$bestmemit[,1],lty=2,lwd=1.8)
      box(lwd=2.2)
      plot(DEoptim.output$member$bestmemit[,2],pch=21,bg="dark grey",col="dark grey",
           cex=1.2,xlab="Iteration",ylab="Rd",cex.axis=cexaxis,cex.lab=cexlab)
      lines(DEoptim.output$member$bestmemit[,2],lty=2,lwd=1.8)
      box(lwd=2.2)
      plot(DEoptim.output$member$bestmemit[,3],pch=21,bg="dark grey",col="dark grey",
           cex=1.2,xlab="Iteration",ylab="Jmax",cex.axis=cexaxis,cex.lab=cexlab)
      lines(DEoptim.output$member$bestmemit[,3],lty=2,lwd=1.8)
      box(lwd=2.2)
    } else if (f.model==2){
      par(mfrow=c(3,1),mar=c(4,4.1,1,2)) #b, l, t, r
      plot(DEoptim.output$member$bestmemit[,1],pch=21,bg="dark grey",col="dark grey",
           cex=1.2,xlab="Iteration",ylab="Vcmax",cex.axis=cexaxis,cex.lab=cexlab)
      lines(DEoptim.output$member$bestmemit[,1],lty=2,lwd=1.8)
      box(lwd=2.2)
      plot(DEoptim.output$member$bestmemit[,2],pch=21,bg="dark grey",col="dark grey",
           cex=1.2,xlab="Iteration",ylab="Rd",cex.axis=cexaxis,cex.lab=cexlab)
      lines(DEoptim.output$member$bestmemit[,2],lty=2,lwd=1.8)
      box(lwd=2.2)
    } else if (f.model==3){
      par(mfrow=c(3,1),mar=c(4,4.1,1,2)) #b, l, t, r
      plot(DEoptim.output$member$bestmemit[,1],pch=21,bg="dark grey",col="dark grey",
           cex=1.2,xlab="Iteration",ylab="Jmax",cex.axis=cexaxis,cex.lab=cexlab)
      lines(DEoptim.output$member$bestmemit[,1],lty=2,lwd=1.8)
      box(lwd=2.2)
      plot(DEoptim.output$member$bestmemit[,2],pch=21,bg="dark grey",col="dark grey",
           cex=1.2,xlab="Iteration",ylab="Rd",cex.axis=cexaxis,cex.lab=cexlab)
      lines(DEoptim.output$member$bestmemit[,2],lty=2,lwd=1.8)
      box(lwd=2.2)
    }
    
    # A-Ci diagnostic fit fig
    plotCi = seq(-2,2200,2)
    Oxygen <- mm.constants$Oxygen
    loc2 <- match(c("KC","KO","GSTAR","TLEAF"),toupper(names(data)))
    Kc <- mean(data[,loc2[1]])
    Ko <- mean(data[,loc2[2]])
    Gstar <- mean(data[,loc2[3]])
    Tleaf <- mean(data[,loc2[4]])
    Vcmax.plot <- params[[1]]
    Jmax.plot <- params[[2]]
    Rd.plot <-params[[3]]
    RMSE.plot <- params[[4]]
    
    par(mfrow=c(1,1),mar=c(5,5,2,1)) # bot, left
    ylim <- range(data[,loc1[2]])
    plot(data[,loc1[1]],data[,loc1[2]], main=paste(file), xlab="Ci", ylab="Photo", cex.lab=2,cex=1.8,
         xlim=c(1.1,2000),ylim=c(0,ylim[2]+3))
    legend("bottomright",legend=c(paste("Tleaf =",round(Tleaf,2)),paste("Vcmax =",round(Vcmax.plot,2)),
                                  paste("Jmax = ",round(Jmax.plot,2)),paste("Rd = ",round(Rd.plot,4)),
                                  paste("RMSE = ",round(RMSE.plot,2))),
           bty="n",cex=2)
    legend("topleft",legend=c("Rubisco","RuBP","Photo"),lty=c(2,2,1),
           col=c("dark blue","dark red","dark grey"),bty="n",lwd=6.3,seg.len=3.5,cex=1.5)
    if (Vcmax.plot!=-9999){
      lines(plotCi,(Vcmax.plot*(plotCi-Gstar)/(plotCi+(Kc*(1+Oxygen/Ko))))-Rd.plot,lwd=5,col="dark blue",lty=2)
    }
    if (Jmax.plot!=-9999){
      lines(plotCi,((Jmax.plot*(plotCi-Gstar))/((4.5*plotCi)+(10.5*Gstar)))-Rd.plot,lwd=5,col="dark red",lty=2)
    }
    if (Vcmax.plot!=-9999 & Jmax.plot!=-9999){
      lines(plotCi,pmin(Vcmax.plot*(plotCi-Gstar)/(plotCi+(Kc*(1+Oxygen/Ko))),
                        (Jmax.plot*(plotCi-Gstar))/((4.5*plotCi)+(10.5*Gstar)))-Rd.plot,col="dark grey",lwd=2.0)
    }
    box(lwd=2.2)
    dev.off()
    
    
  } else if (type=="A-Q"){
    
    loc1 <- match(c("PARI","PHOTO"),toupper(names(data)))
    
    # Figures
    pdf(paste(outdir,sep,file,".pdf",sep=""),height=8,width=10)
    plot(data[,loc1[1]],data[,loc1[2]],pch=21,bg="grey70",cex=3,cex.axis=cexaxis,xlim=c(0,range(data[,loc1[1]])[2]),
         ylim=c(0,range(data[,loc1[2]])[2]+1),cex.lab=cexlab,xlab="PARi",ylab="Photo",main=paste(file))
    box(lwd=2.2)
    
    # DEoptim trace plots
    if (f.model==1){
      par(mfrow=c(2,1),mar=c(4,4.1,1,2)) #b, l, t, r
      plot(DEoptim.output$member$bestmemit[,1],pch=21,bg="dark grey",col="dark grey",
           cex=1.2,xlab="Iteration",ylab="Amax",cex.axis=cexaxis,cex.lab=cexlab)
      lines(DEoptim.output$member$bestmemit[,1],lty=2,lwd=1.8)
      box(lwd=2.2)
      plot(DEoptim.output$member$bestmemit[,2],pch=21,bg="dark grey",col="dark grey",
           cex=1.2,xlab="Iteration",ylab="Convex",cex.axis=cexaxis,cex.lab=cexlab)
      lines(DEoptim.output$member$bestmemit[,2],lty=2,lwd=1.8)
      box(lwd=2.2)
    }
    
    # A-Q diagnostic fit fig
    plotPARi <- seq(-5,2000,2)
    loc2 <- match(c("PARI","TLEAF"),toupper(names(data)))
    PARi <- data[,loc2[1]]
    Tleaf <- mean(data[,loc2[2]])
    aQY.plot <- params[[1]]   
    LCP.plot <- params[[2]]
    Amax.plot <- params[[3]]
    convex.plot <- params[[4]]
    LUE.plot <- params[[5]]
    Rd.plot <- params[[6]]
    RMSE.plot <- params[[7]]
    IaQY <- aQY.plot*plotPARi
    
    par(mfrow=c(1,1),mar=c(5,5,2,1)) # bot, left
    ylim <- range(data[,loc1[2]])
    plot(data[,loc1[1]],data[,loc1[2]], main=paste(file), xlab="PARi", ylab="Photo", cex.lab=2,cex=3,
         xlim=c(1.1,1500),ylim=c(ylim[1],ylim[2]+3))
    legend("bottomright",legend=c(paste("Tleaf =",round(Tleaf,2)),paste("LCP =",round(LCP.plot,2)),
                                  paste("aQY =",round(aQY.plot,3)),paste("Amax = ",round(Amax.plot,2)),
                                  paste("Convex = ",round(convex.plot,2)),paste("LUE.max = ",round(LUE.plot,4)),
                                  paste("Rd = ",round(Rd.plot,4)),
                                  paste("RMSE = ",round(RMSE.plot,2))),bty="n",cex=2)
    legend("topleft",legend=c("aQY","Photo"),lty=c(2,1),
           col=c("dark blue","dark grey"),bty="n",lwd=6.3,seg.len=3.5,cex=1.5)
    box(lwd=2.2)
    
    if (Amax.plot!=-9999){
      lines(plotPARi[plotPARi<=100],aQY.plot*plotPARi[plotPARi<=100] - Rd.plot,lwd=5,col="dark blue",lty=2)
      lines(plotPARi,(((IaQY+Amax.plot)-sqrt((IaQY+Amax.plot)^2-(4*convex.plot*IaQY*Amax.plot)))/(2*convex.plot))-Rd.plot,
            col="dark grey",lwd=2.0)
      box(lwd=2.2)
      
      # LUE plot
      plotPARi <- seq(0,2000,1)
      IaQY <- aQY.plot*plotPARi
      plot(plotPARi,((((IaQY+Amax.plot)-sqrt((IaQY+Amax.plot)^2-(4*convex.plot*IaQY*Amax.plot)))/(2*convex.plot))-Rd.plot)
           /plotPARi,type="l",ylab="LUE",xlab="PARi",ylim=c(0,0.1),lwd=3)
      points(which.max(((((IaQY+Amax.plot)-sqrt((IaQY+Amax.plot)^2-(4*convex.plot*IaQY*Amax.plot)))/(2*convex.plot))-Rd.plot)
                       /plotPARi),max(((((IaQY+Amax.plot)-sqrt((IaQY+Amax.plot)^2-(4*convex.plot*IaQY*Amax.plot)))/(2*convex.plot))-Rd.plot)
                                      /plotPARi),pch=20,bg="dark grey",cex=3.3)
      legend("topleft",legend="A/PARi",col="black",bty="n",lwd=6.3,seg.len=3.5,cex=1.5)
      legend("topright",legend=paste("LUE.max = ",round(LUE.plot,4)),,bty="n",cex=2)
      box(lwd=2.2)
    } else {
      lines(plotPARi[plotPARi<=100],aQY.plot*plotPARi[plotPARi<=100] - Rd.plot,lwd=5,col="dark blue",lty=2)
    }
    dev.off()
  } # End A-Ci / A-Q if/else
  
} # End of function
#==================================================================================================#


#--------------------------------------------------------------------------------------------------#
##'
##' Temperature response diagnostic figures
##' @author Shawn P. Serbin
##'  
##'  very basic, needs updating
##'  
plot.temp.response <- function(type=c("Arrhenius","Peaked","June"),var,data,DEoptim.output,params,
                               outdir,file,year){
  type <- match.arg(type)
  sep <- .Platform$file.sep
  param.length <- length(DEoptim.output$optim$bestmem)
  
  ### Plot params
  cexaxis <- 1.2
  cexlab <- 1.4
  
  if (type=="Arrhenius"){
    loc1 <- match(c("TLEAF",paste(toupper(var))),toupper(names(data)))
    pdf(paste(outdir,sep,file,".pdf",sep=""),height=8,width=10)
    plot(data[,loc1[1]],data[,loc1[2]],pch=21,bg="grey70",cex=3,cex.axis=cexaxis,xlim=c(range(data[,loc1[1]])[1]-2,
                                                                                        range(data[,loc1[1]])[2]),
         ylim=c(range(data[,loc1[2]])[1]-5,range(data[,loc1[2]])[2]+5),cex.lab=cexlab,xlab="Tleaf",
         ylab="Vcmax",main=paste(file))
    box(lwd=2.2)
    
    # Best member parameter trace plots
    par(mfrow=c(2,1),mar=c(4,4.1,1,2)) #b, l, t, r
    plot(DEoptim.output$member$bestmemit[,1],pch=21,bg="dark grey",col="dark grey",
         cex=1.2,xlab="Iteration",ylab="Vc25",cex.axis=cexaxis,cex.lab=cexlab)
    lines(DEoptim.output$member$bestmemit[,1],lty=2,lwd=1.8)
    box(lwd=2.2)
    plot(DEoptim.output$member$bestmemit[,2],pch=21,bg="dark grey",col="dark grey",
         cex=1.2,xlab="Iteration",ylab="E",cex.axis=cexaxis,cex.lab=cexlab)
    lines(DEoptim.output$member$bestmemit[,2],lty=2,lwd=1.8)
    box(lwd=2.2)
    
    # Plot output
    plotTemp = seq(10,80,0.2)
    par(mfrow=c(1,1),mar=c(5,5,2,1))
    ylim <- range(Photo.Parameter[,1])
    plot(Tleaf, Photo.Parameter[,1], main=paste0(year),xlab="Temp C", ylab=var, cex.lab=2, pch=21, cex=2.2,bg="grey50",
         xlim=c(15,40),ylim=c(ylim[1]-30,ylim[2]+10))
    legend("topleft",legend=c(paste0(var,"25 = ",round(params[1],2)), paste0("Ev = ",round(params[2],2)),
                              paste0("RMSE = ",round(params[3],2))),bty="n",cex=1.5)
    lines(plotTemp,params[1,1]*exp((params[1,2]*((plotTemp+273.15)-(stand.temp+273.15)))/((stand.temp+273.15)*R*(plotTemp+273.15))),
          ,lwd=4,col="black")
    box(lwd=2.2)
    
    dev.off()
  }
  
  
} # End of function
#==================================================================================================#


####################################################################################################
### EOF.  End of R script file.              
####################################################################################################