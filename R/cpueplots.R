
#' @title addlnorm - estimates a log-normal distribution from output of hist.
#'
#' @description  addlnorm - estiamtes a log-normal distribution from output of
#'    a histogram of a data set.
#' @param inhist - is the output from a call to 'hist' (see examples)
#' @param xdata -  is the data that is being plotted in the histogram.
#' @param inc - defaults to a value of 0.01; is the fine grain increment used to
#'    define the normal curve. The histogram will be coarse grained relative to
#'    this.
#' @return a 4 x N matrix of x and y values to be used to plot the fitted normal
#'    probability density function.Combined with estiamtes of mean(log(indata))
#'    and log(sd(indata))
#' @export addlnorm
#'
#' @examples
#' egdata <- rlnorm(200,meanlog=0.075,sdlog=0.5)
#' outh <- hist(egdata,main="",col=2,breaks=seq(0,8,0.2))
#' ans <- addlnorm(outh,egdata)
#' lines(ans[,"x"],ans[,"y"],lwd=2,col=4)
addlnorm <- function(inhist,xdata,inc=0.01) {
  lower <- inhist$breaks[1]
  upper <- tail(inhist$breaks,1)
  cw <- inhist$breaks[2]-inhist$breaks[1]
  x <- seq(lower,upper, inc) #+ (cw/2)
  avCE <- mean(log(xdata),na.rm=TRUE)
  sdCE <- sd(log(xdata),na.rm=TRUE)
  N <- length(xdata)
  ans <- cbind(x,(N*cw) * stats::dlnorm(x,avCE,sdCE),avCE,sdCE)
  colnames(ans) <- c("x","y","avCE","sdCE")
  return(ans)
} # end of addlnorm

#' @title addnorm - adds a normal distribution to a histogram of a data set.
#'
#' @description  addnorm - adds a normal distribution to a histogram of a data
#'    set. This is generally to be used to illustrate whether log-transformation
#'    normalizes a set of catch or cpue data.
#' @param inhist - is the output from a call to 'hist' (see examples)
#' @param xdata -  is the data that is being plotted in the histogram.
#' @param inc - defaults to a value of 0.01; is the fine grain increment used to
#'    define the normal curve. The histogram will be coarse grained relative to
#'    this.
#' @return a list with a vector of 'x' values and a vector of 'y' values (to be
#'    used to plot the fitted normal probability density function), and a vector
#'    used two called 'stats' containing the mean and sandard deviation of the
#'    input data
#' @export addnorm
#' @examples
#' x <- rnorm(1000,mean=5,sd=1)
#' dev.new(height=6,width=4,noRStudioGD = TRUE)
#' par(mfrow= c(1,1),mai=c(0.5,0.5,0.3,0.05))
#' par(cex=0.85, mgp=c(1.5,0.35,0), font.axis=7)
#' outH <- hist(x,breaks=25,col=3,main="")
#' nline <- addnorm(outH,x)
#' lines(nline$x,nline$y,lwd=3,col=2)
#' print(nline$stats)
addnorm <- function(inhist,xdata,inc=0.01) {
  lower <- inhist$breaks[1]
  upper <- tail(inhist$breaks,1)
  cw <- inhist$breaks[2]-inhist$breaks[1]
  x <- seq(lower,upper, inc) #+ (cw/2)
  avCE <- mean(xdata,na.rm=TRUE)
  sdCE <- sd(xdata,na.rm=TRUE)
  N <- length(xdata)
  ans <- list(x=x,y=(N*cw)*dnorm(x,avCE,sdCE),stats=c(avCE,sdCE,N))
  return(ans)
} # end of addnorm

#' @title diagnosticPlot produces diagnostic plots of the regression.
#'
#' @description  diagnosticPlot produces diagnostic plots of the regression.
#'   these include a plot of the residuals agsinst the predicted values, the
#'   normal Q-Q plot, a histogram of the Observed Log(CPUE) with a fitted
#'   normal distribution, and a histogram of the predicted Log(CPUE) with a
#'   fitted normal distribution and the normal distribution from the observed
#'   data.
#' @param inout is the list from standLM.
#' @param inmodel is a single number identifying the model inside 'inout' to
#'    plotted.
#' @param indat is the data.frame used in the standardization.
#' @param inlabel provides the means of giving a title to the 2 x 2 plot
#' @param height default = 6; the height of the output plot
#' @param width default = 8; the width of the output plot
#' @return a 2 x 2 plot of the residuals, the Q-Q plot, and the distributions
#'    of the logged observed and predicted values of CPUE.
#' @export diagnosticPlot
#' @examples
#' \dontrun{
#' data(sps)
#' splabel = "Species"
#' sps$Year <- factor(sps$Year)
#' sps$Month <- factor(sps$Month)
#' sps$Vessel <- factor(sps$Vessel)
#' labelM <- c("Year","Vessel","Month")
#' mods <- makemodels(labelM)
#' out <- standLM(mods,sps,splabel)
#' diagnosticPlot(out,3,sps,splabel)
#' }
diagnosticPlot <- function(inout, inmodel, indat, inlabel="",
                           height=6,width=8)  {
  model <- inout$Models[[inmodel]]
  resids <- model$residuals
  fitted.vals <- model$fitted.values
  labs <- 1.1
  dev.new(height=height,width=width,noRStudioGD = TRUE)
  par(mfrow= c(2,2))
  par(mai=c(0.5,0.5,0.3,0.05),oma=c(0,0,0,0))
  par(cex=0.85, mgp=c(1.5,0.35,0), font.axis=7)
  plot(fitted.vals,resids, main="",
       ylab="", xlab="")
  title(main=list(inlabel, cex=labs, font=7),
        xlab=list("Predicted Values", cex=labs, font=7),
        ylab=list("Residuals", cex=labs, font=7))
  abline(h=0.0,col=2)
  qqnorm(resids, ylab=list("Std Deviance Resid.", cex=labs, font=7),
         xlab=list("Theoretical Quantiles", cex=labs, font=7),
         main=list("Normal Q-Q Plot",cex=labs,font=7))
  qqline(resids, col=2,lwd=2)
  rug(resids,ticksize=0.02,col="royal blue")
  abline(v=c(-2.0,2.0),col="grey")
  par(mai=c(0.5,0.5,0.1,0.2))
  meance <- mean(indat$LnCE,na.rm=TRUE)
  sdce <- sd(indat$LnCE,na.rm=TRUE)
  lower <- floor(min(indat$LnCE,na.rm=TRUE))
  upper <- ceiling(max(indat$LnCE,na.rm=TRUE))
  lowerf <- floor(min(fitted.vals,na.rm=TRUE))
  upperf <- ceiling(max(fitted.vals,na.rm=TRUE))
  minx <- min(lower,lowerf)
  maxx <- max(upper,upperf)
  cebin <- 0.25
  x <- seq(minx,maxx,0.01)
  y <- dnorm(x,mean=meance,sd=sdce)
  bins <- seq(minx,maxx,cebin)
  histout <- hist(indat$LnCE,breaks=bins,main="",xlab="",ylab="")
  totn <- sum(histout$counts)
  ymax <- max(histout$counts)
  lines(x,(totn/(1/cebin))*y,col=4,lwd=2)
  title(xlab=list("Observed Log(CE)", cex=labs, font=7),
        ylab=list("Frequency", cex=labs, font=7))
  text(minx+5*cebin,0.975*ymax,paste("Mean  ",round(meance,3),sep=""),cex=labs,font=7)
  text(minx+5*cebin,0.9*ymax,paste("StDev ",round(sdce,3),sep=""),cex=labs,font=7)
  
  meancef <- mean(fitted.vals,na.rm=TRUE)
  sdcef <- sd(fitted.vals,na.rm=TRUE)
  yf <- dnorm(x,mean=meancef,sd=sdcef)
  histoutf <- hist(fitted.vals,breaks=bins,main="",xlab="",ylab="")
  totnf <- sum(histoutf$counts)
  ymaxf <- max(histoutf$counts)
  lines(x,(totnf/(1/cebin))*yf,col=2,lwd=2)
  lines(x,(totn/(1/cebin))*y,col=4,lwd=2)
  title(xlab=list("Predicted Log(CE)", cex=labs, font=7),
        ylab=list("Frequency", cex=labs, font=7))
  text(lower+6*cebin,0.975*ymaxf,paste("Mean  ",round(meancef,3),sep=""),cex=labs,font=7)
  text(lower+6*cebin,0.9*ymaxf,paste("StDev ",round(sdcef,3),sep=""),cex=labs,font=7)
}  # end of diagnosticPlot

#' @title impactplot - Plots relative contribution to CPUE trend of each Factor
#'
#' @description Calculates and plots out the contribution to a standardized CPUE
#'    trend of each of the factors included in the standardization. The number on
#'    each graph is the sum of squared differences between the previous model
#'    and the current model - as a measure of the relative effect of the
#'    factor concerned. Positive effects are shown as blue bars, negative
#'    effects are shown as red bars.
#' @param inout is the list output from the standLM function containing
#'    the standardization
#' @param mult default value is 3 - it is used to scale the bars on the plot.
#' @param FY - are the years fishing seasons or calender years; defaults to FALSE,
#'    which assumes calendar years for the x-axis.
#' @return generates a plot of the impact of each factor on the CPUE trend, in
#'   addition, the function returns, invisibly, two tables as a list, the first
#'   $result summarizes the impact of each factor included in the analysis,
#'   including the sum of absolute differences between each model and the one
#'   before it (the top plot is for the complete model). The second table
#'   $deviates lists the actual differences between the particular vaiables
#'   and the previous models.
#' @export impactplot
#' @examples
#' \dontrun{
#' data(sps)
#' splabel = "Species"
#' sps$Year <- factor(sps$Year)
#' sps$Month <- factor(sps$Month)
#' sps$Vessel <- factor(sps$Vessel)
#' labelM <- c("Year","Vessel","Month")
#' mods <- makemodels(labelM)
#' out <- standLM(mods,sps,splabel)
#' impactplot(out,mult=3)
#' }
impactplot <- function(inout,mult=3,FY=FALSE) {
  incoef <- inout[[1]]
  label <- colnames(incoef)                    # get the factor names
  nmods <- length(label)   # How many factors
  if (FY) {
    fyears <- rownames(incoef)
    ny <- length(fyears)
    years <- numeric(ny)
    for (i in 1:ny) years[i] <- as.numeric(unlist(strsplit(fyears[i],"/"))[1])
  } else {
    years <- as.numeric(rownames(incoef))        # get the years used
  }
  nyr <- length(years)
  result <- matrix(0,nrow=nmods,ncol=4,dimnames=list(label,c("SumAbsDev","%Diff","adj_r2","DeltaV")))
  deviates <- matrix(0,nrow=nyr,ncol=nmods,dimnames=list(years,label))
  par(mfrow=c(nmods,1))
  par(mai=c(0.25,0.35,0.05,0.05),oma=c(0.0,1.0,0.05,0.25))
  par(cex=0.85, mgp=c(1.5,0.35,0), font.axis=7)
  ymax <- max(incoef)*1.15   # leave enough room for the labels
  # plot the Geometric Mean trend = model 1
  plot(years,incoef[,1],type="l",ylim=c(0,ymax),ylab="",xlab="",lwd=2,yaxs="i",xaxs="r")
  lines(years,incoef[,nmods],lwd=2,col=2)
  devs <- (incoef[,nmods]-incoef[,1])        # calc diff between trends
  deviates[,1] <- devs
  pick <- which(devs < 0)                    # plot negative deviates as red
  if (length(pick)) { lines(years[pick],mult*abs(devs[pick]),type="h",lwd=4,col=2) }
  pick <- which(devs >= 0)                   # positive deviates as blue
  if (length(pick)) { lines(years[pick],mult*devs[pick],type="h",lwd=4,col=4) }
  sumdev <- sum(abs(devs))
  result[1,1] <- sumdev
  mtext(paste(label[1],round(sumdev,3),sep="  "),3,line=-1,font=6,cex=0.9)
  for (gr in 2:nmods) {                  # now step through the remaining factors
    plot(years,incoef[,gr],type="l",ylim=c(0,ymax),ylab="",xlab="",lwd=2,yaxs="i",xaxs="r")
    lines(years,incoef[,gr-1],lwd=2,col="grey")     # add previous factor
    devs <- (incoef[,gr]-incoef[,gr-1])        # calc diff between trends
    deviates[,gr] <- devs
    pick <- which(devs < 0)                    # plot negative deviates as red
    if (length(pick)) { lines(years[pick],mult*abs(devs[pick]),type="h",lwd=4,col=2) }
    pick <- which(devs >= 0)                   # positive deviates as blue
    if (length(pick)) { lines(years[pick],mult*devs[pick],type="h",lwd=4,col=4) }
    sumdev <- sum(abs(devs))
    result[gr,1] <- sumdev
    mtext(paste(label[gr],round(sumdev,3),sep="  "),3,line=-1,font=6,cex=0.9)
  }
  vlabel <- paste("Standardized Catch Rates ",inout$Label,sep="")
  mtext(vlabel,2,outer=TRUE,line=-0.5,font=6,cex=1.1)
  for (gr in 3:nmods) {                       # calc % difference between abs diff
    result[gr,2] <- 100*result[gr,1]/result[gr-1,1]
  }
  pick <- which(rownames(inout[[3]]) == "adj_r2")
  result[,3] <- inout[[3]][pick,]
  result[,4] <- inout[[3]][(pick+1),]
  deviates[,1] <- rowSums(deviates) # this sums the raw deviates to show the final
  # total change from the full model to the geometric mean
  ans <- list(result,deviates)
  names(ans) <- c("result","deviates")
  return(invisible(ans))
} # end of impactplot

#' @title inthist - a replacement for the hist function for use with integers
#'
#' @description inthist - a replacement for the hist function for use with
#'    integers because the ordinary function fails to count them correctly at
#'    the end. The function is designed for integers and if it is given real
#'    numbers it will issue a warning and then round all values before plotting.
#' @param x - the vector of integers to be counted and plotted OR a matrix of
#'     values in column 1 and counts in column 2
#' @param col - the colour of the fill; defaults to black = 1, set this to 0
#'    for an empty bar, but then give a value for border
#' @param border - the colour of the outline of each bar defaults to col
#' @param width - denotes the width of each bar; defaults to 1, should be >0
#'    and <= 1
#' @param xlabel - the label for the x axis; defaults to ""
#' @param ylabel - the label for the y axis; defaults to ""
#' @param main - the title for the individual plot; defaults to ""
#' @param lwd - the line width of the border; defaults to 1
#' @param xmin - sets the lower bound for x-axis; used to match plots
#' @param xmax - sets the upper bound for x axis; used with multiple plots
#' @param ymax - enables external control of the maximum y value; mainly of
#'    use when plotting multiple plots together.
#' @param plotout - plot the histogram or not? Defaults to TRUE
#' @param prop - plot the proportions rather than the counts
#' @param inc - sets the xaxis increment; used to customize the axis;
#'    defaults to 1.
#' @param xaxis - set to FALSE to define the xaxis outside of inthist;
#'    defaults to TRUE
#' @return a matrix of values and counts is returned invisibly
#' @export inthist
#' @examples
#' x <- trunc(runif(1000)*10) + 1
#' plotprep(width=6,height=4)
#' inthist(x,col="grey",border=3,width=0.75,xlabel="Random Uniform",
#'         ylabel="Frequency")
inthist <- function(x,col=1,border=NULL,width=1,xlabel="",ylabel="",
                    main="",lwd=1,xmin=NA,xmax=NA,ymax=NA,plotout=TRUE,
                    prop=FALSE,inc=1,xaxis=TRUE) {
  if (class(x) == "matrix") {
    counts <- x[,2]
    values <- x[,1]
  } else {
    counts <- table(x)
    if (length(counts) == 0) stop("No data provided \n\n")
    values <- as.numeric(names(counts))
  }
  
  if (sum(!(abs(values - round(values)) < .Machine$double.eps^0.5)) > 0) {
    warning("Attempting to use 'inthist' with non-integers; Values now rounded \n")
    values <- round(values,0)
  }
  if ((width <= 0) | (width > 1)) {
    warning("width values should be >0 and <= 1")
    width <- 1
  }
  counts <- as.numeric(counts)
  nct <- length(counts)
  propor <- counts/sum(counts,na.rm=TRUE)
  if (is.na(xmin)) xmin <- min(values)
  if (is.na(xmax)) xmax <- max(values)
  if (prop) {
    outplot <- propor
  } else {
    outplot <- counts
  }
  if (is.na(ymax)) {
    if (nchar(main) > 0) {
      ymax <- max(outplot) * 1.15
    } else {
      ymax <- max(outplot) * 1.05
    }
  }
  if (plotout) {
    plot(values,outplot,type="n",
         xlim=c((xmin-(width*0.75)),(xmax+(width*0.75))),
         xaxs="r",ylim=c(0,ymax),yaxs="i",xlab="",ylab="",xaxt="n")
    if (xaxis) axis(side=1,at=seq(xmin,xmax,inc),labels=seq(xmin,xmax,inc))
    if (length(counts) > 0) {
      for (i in 1:nct) {  # i <- 1
        x1 <- values[i] - (width/2)
        x2 <- values[i] + (width/2)
        x <- c(x1,x1,x2,x2,x1)
        y <- c(0,outplot[i],outplot[i],0,0)
        if (is.null(border)) border <- col
        polygon(x,y,col=col,border=border,lwd=lwd)
      }
      title(ylab=list(ylabel, cex=1.0, font=7),
            xlab=list(xlabel, cex=1.0, font=7))
      if (nchar(main) > 0) mtext(main,side=3,line=-1.0,outer=FALSE,cex=0.9)
    }
  } # end of if-plotout
  if (length(counts) > 0) {
    answer <- cbind(values,counts,propor);
    rownames(answer) <- values
    colnames(answer) <- c("values","counts","proportion")
  } else { answer <- NA  }
  class(answer) <- "inthist"
  return(invisible(answer))
}  # end of inthist

#' @title lefthist draws a histogram up the y-axis
#'
#' @description lefthist translates a histogram from along the x-axis to
#'     flow along the y-axis - it transposes a histogram.
#'
#' @param x a vector of the data to be plotted
#' @param bins the breaks from the histogram, can be a single number of a
#'     sequence of values; defaults to 25
#' @param mult the multiplier for the maximum count in the histogram. Becomes
#'     the upper limit of teh x-axis.
#' @param col the colour for the histogram polygons; default = 2
#' @param lwd the line width for each polygon; default = 1
#' @param width the width of each bar in teh histogram; default = 0.9
#' @param border the colour for the border line; default = 1 = black
#' @param xinc the step size for the x-axis (counts) labels; default= NA,
#'     which means the increment will equal the bin width.
#' @param yinc the step size for the y-axis (breaks) labels; default= 1.
#' @param title the title for the left-histogram; defaults to ""
#' @param xlabel the xlab; defaults to ""
#' @param ylabel the ylab; defaults to "Frequency"
#' @param cex the size of text in teh plot. defaults = 1.0
#' @param textout prints input data range to console; default = FALSE
#' @param hline if this has a value a horizontal line will be plotted;
#'     default = NA
#'
#' @return the output from hist but done so invisibly.
#' @export
#'
#' @examples
#' dat <- rnorm(1000,mean=5,sd=1)
#' dev.new(width=6,height=4,noRStudioGD = TRUE)
#' par(mai=c(0.45,0.45,0.05,0.05))
#' lefthist(dat)
#' lefthist(dat,textout=TRUE,width=0.8,border=3)
lefthist <- function(x,bins=25,mult=1.025,col=2,lwd=1,width=0.9,border=1,
                     xinc=1,yinc=NA,title="",xlabel="Frequency",ylabel="",
                     cex=1.0,textout=FALSE,hline=NA) {
  outh <- hist(x,breaks=bins,plot=FALSE)
  cw <- outh$breaks[2]-outh$breaks[1]
  newcount <- c(outh$counts,0)
  ymax <- max(newcount,na.rm=TRUE) * mult
  nvalues <- length(newcount)
  values <- outh$breaks
  if (is.na(yinc)) yinc <- values[2] - values[1]
  xlabs <- seq(0,(ymax+(2 * xinc)),xinc)
  xlabs <- xlabs[xlabs < ymax]
  plot(seq(0,ymax,length=nvalues),values,type="n",xlim=c(0,ymax),
       xaxs="i",ylim=c(range(values)),yaxs="r",xlab="",ylab="",xaxt="n",yaxt="n")
  grid()
  axis(side=1,at=xlabs,labels=xlabs)
  title(ylab=list(ylabel,cex=cex),xlab=list(xlabel,cex=cex))
  values1 <- seq(values[1],values[nvalues],yinc)
  axis(side=2,at=values1,labels=values1,cex.lab=cex)
  for (i in 1:nvalues) {  # i <- 1
    y1 <- values[i]
    y2 <- values[i] + (cw * width)
    yplot <- c(y1,y1,y2,y2,y1)
    xplot <- c(0,newcount[i],newcount[i],0,0)
    if (is.null(border)) border <- col
    polygon(xplot,yplot,col=col,border=border,lwd=lwd)
  }
  if (!is.na(hline)) abline(h=hline,col=(col+2))
  
  if (textout) cat("  input data range: ",range(x,na.rm=TRUE),"\n\n")
  return(invisible(outh))
}  # end of lefthist

#' @title plotdata generates graphs of CE and log-transformed CE data.
#'
#' @description  plotdata generates graphs of CE and log-transformed CE
#'   data. Included in the graph of the log-transformed cpue data is a
#'   fitted normal distribution with the associated bias-corrected average.
#' @param indat is the data.frame containing the data being standardized.
#' @param dependent variable name; defaults to LnCE.
#' @return a plot of the untransformed cpue data and of the log-transformed
#'   data. Included on the log-transformed data is a fitted normal
#'   distribution and the bias-corrected mean estiamte of average CPUE.
#' @export plotdata
#' @examples
#' \dontrun{
#' data(sps)
#' pick <- which(sps$Year == 2005)
#' sps2 <- droplevels(sps[pick,])
#' plotprep()
#' plotdata(sps2)
#' plotdata(sps2[sps2$Depth > 200,])
#' }
plotdata <- function(indat,dependent="LnCE") {
  colnames(indat) <- tolower(colnames(indat))
  dependent <- tolower(dependent)
  av1 <- mean(indat[,dependent],na.rm=TRUE)
  sd1 <- sd(indat[,dependent],na.rm=TRUE)
  n <- length(indat[,dependent])
  par(mfrow= c(2,1))
  par(mai=c(0.3,0.3,0,0.1), oma=c(0,1.0,0.0,0.0))
  par(cex=0.9, mgp=c(1.5,0.3,0), font.axis=7)
  outce <- hist(indat$ce,breaks=25,main="",ylab="",xlab="",col=2)
  outlnce <- hist(indat$lnce,breaks=25,main="",ylab="",xlab="",col=2)
  low <- outlnce$breaks[1]
  upp <- tail(outlnce$breaks,1)
  inc <- outlnce$breaks[2] - low
  ymax <- max(outlnce$counts)
  span <- seq(low,upp,0.05)
  lines(span,(n*inc)*dnorm(span,av1,sd1),lwd=3,col=1)
  abline(v=av1,col=3,lwd=2)
  mtext("Frequency",side=2,outer=TRUE,line=0.0,font=7,cex=1.0)
  label <- paste0("AvCE  ",round(exp(av1 + (sd1*sd1)/2),3))
  text(low,0.7*ymax,label,cex=1.0,font=7,pos=4)
}  # end of plotdata


#' @title plotfishery plots an array of data for the fishery.
#'
#' @description plotfishery plots an array of data for the fishery. Including
#'     the number of records by depth, the catch-by-cept by zone, the number
#'     of vessels by year, the number of records by year, the total catch by
#'     all methods, gears, and zones by year, combined with the selected
#'     catches-by year, then the selected catches-by-year to increase the
#'     y-axis scale if necessary to calrify the total catch < 30kg by year.
#'
#' @param InData2 is he selected data from selectdata, usually sps1
#' @param Datain is the datasum from makedatasum
#' @param spsname is the original splabel, usually splabel
#' @param Ldep lowest selected depth usually Ldepth, defaults to 0
#' @param Udep deepest selected depth, usually Udepth, defaults to 1000
#' @param zones defaults to 'Zone' to represent the SESSF zones, but if the
#'     fishery relates to the Orange Roughy zones then use 'ORzone', and if it
#'     relates to the Shark fishery use "SharkRegion'
#' @param legpos defaults to "L", which will put the legend on the left of the
#'     plot. Any other entry, such as "R" will place it on the right hand side
#' @param legval defaults to 0.6, which should work for most plots unless the
#'     x-axis does not start at 0, in which case a larger value may be required.
#'
#' @return currently nothing but this does generate a 3 x 2 plot
#' @export
#'
#' @examples
#' \dontrun{
#' print("still to develop a full example")
#' }
plotfishery <- function(InData2,Datain,spsname="",
                        Ldep=0,Udep=1000,zones="Zone",legpos="L",legval=0.6) {
  years <- as.numeric(rownames(Datain))
  fsize <- 0.85
  par(mfrow= c(3,2))
  par(mai=c(0.45,0.5,0.15,0.15), oma=c(0,0,0,0),xaxs="i")
  par(cex=0.85, mgp=c(1.5,0.35,0), font.axis=7, font=7)
  # plot records by Depth
  bins <- seq(Ldep,Udep+10,20)
  hist(InData2$Depth, breaks = bins, main="", xlab="", ylab="")
  title(ylab=list("Records", cex=fsize, col=1, font=7),
        xlab=list("Depth m", cex=fsize, col=1, font=7),
        main=list("Records by Depth", cex=fsize, col=1, font=7))
  # Plot the catch by depth by zone
  if (zones %in% c("Zone","ORzone","SharkRegion")) {
    table5 <- tapply(InData2$catch_kg,
                     list(InData2$DepCat,InData2[,zones]),sum,na.rm=TRUE)/1000.0
    zonein <- sort(unique(InData2[,zones]))
  } else {
    stop("Unknown name in the zones variable. Should be either 'Zone',
           'ORzone', or 'SharkRegion'")
  }
  nlist <- dim(table5)
  nline <- nlist[2]
  x <- as.numeric(rownames(table5))
  maxy <- max(table5,na.rm=T)*1.075
  plot(x,table5[,1],type="b",pch=16,cex=1.0,xlab="",ylab="",main="",
       ylim=c(0,maxy),xlim=c(Ldep,Udep),yaxs="i",lwd=2)
  if (nline > 1) {
    for (i in 2:nline) {
      points(x,table5[,i],pch=16,cex=1.0,col=i)
      lines(x,table5[,i],type="l",col=i,lwd=2)
    }
  }
  title(ylab=list("Catch t", cex=fsize, col=1, font=7),
        xlab=list("Depth m", cex=fsize, col=1, font=7),
        main=list("Catch by Depth by Zone", cex=fsize, col=1, font=7))
  legX <- Ldep
  if (legpos != "L") legX <- legval * Udep
  legend(legX,maxy*0.975,zonein,col=c(1:nline),lwd=2,bty="n")
  # Number of Vessels by Year for selected gear and depths
  maxy <- max(Datain[,4],na.rm=T)*1.025
  plot(years,Datain[,4], type="l", col=1, ylim=c(0,maxy),
       lwd=2, xlab="", ylab="", yaxs="i")
  abline(v=c(1991.5,2006.5),col="grey")
  title(ylab=list("Vessel Numbers", cex=fsize, col=1, font=7),
        xlab=list("Year", cex=fsize, col=1, font=7),
        main=list("Vessels by Year", cex=fsize, col=1, font=7))
  # Plot the number of records selected1
  maxy <- max(Datain[,2])*1.025
  plot(years,Datain[,2], type="l", col=1, ylim=c(0,maxy),
       lwd=2, xlab="", ylab="", yaxs="i")
  abline(v=c(1991.5,2006.5),col="grey")
  title(ylab=list("Records", cex=fsize, col=1, font=7),
        xlab=list("Year", cex=fsize, col=1, font=7),
        main=list("Records by Year", cex=fsize, col=1, font=7))
  # Plot the catch history of the fishery, total and selected
  maxy <- max(Datain[,1])*1.025
  plot(years,Datain[,1], type="l", col=1, ylim=c(0,maxy),
       lwd=2, xlab="", ylab="", yaxs="i")
  lines(years,Datain[,3],col=4,lwd=2)
  lines(years,Datain[,8],col=2,lwd=2)
  title(ylab=list("All Catches T", cex=fsize, col=1, font=7),
        xlab=list("Year", cex=fsize, col=1, font=7),
        main=list("Total & Selected Catches", cex=fsize, col=1, font=7))
  # Plot the particular catches from teh selected fishery for clarity
  maxy <- max(Datain[,3])*1.025
  plot(years,Datain[,3], type="l", col=4, ylim=c(0,maxy),
       lwd=2, xlab="", ylab="", yaxs="i")
  lines(years,Datain[,8],col=2,lwd=2)
  abline(v=c(1991.5,2006.5),col="grey")
  title(ylab=list("Catches T", cex=fsize, col=1, font=7),
        xlab=list("Year", cex=fsize, col=1, font=7),
        main=list("Selected Catches", cex=fsize, col=1, font=7))
}    # END OF PLOT FISHERY



#' @title plotprep: sets up a window and the par values for a single plot
#'
#' @description plotprep: sets up a window and the par values for a single plot.
#'   it checks to see if a graphics device is open and opens a new one if not.
#'   This is simply a utility function to save typing the standard syntax.
#'   Some of the defaults can be changed. Typing the name without () will
#'   provide a template for modification. If 'windows' is called repeatedly this
#'   will generate a new active graphics device each time leaving the older ones
#'   inactive but present. For quick exploratory plots this behaviour is not
#'   wanted, hence the check if an active device exists already or not.
#' @param width defaults to 6 inches = 15.24cm - width of plot
#' @param height defaults to 3 inches = 7.62cm - height of plot
#' @param usefont default is 7 (bold Times); 1 = sans serif, 2 = sans serif bold
#' @param cex default is 0.85, the size of font used for text within the plots
#' @param newdev reuse a previously defined graphics device or make a new one;
#'    defaults to TRUE
#' @param filename defaults to "" = do not save to a filename. If a filename is
#' @param resol the resolution of the png file, if there is one, default=300
#' @param verbose comments to the connsole, default = TRUE
#' 
#' @return Checks for and sets up a graphics device and sets the 
#'     default plotting par values. oldpar values returned invisibly 
#' @export plotprep
#' @examples
#' x <- rnorm(1000,mean=0,sd=1.0)
#' plotprep()
#' hist(x,breaks=30,main="",col=2)
plotprep <- function (width=6,height=3.6,usefont=7,cex=0.85,newdev=TRUE, 
                      filename = "", resol = 300, verbose = TRUE) 
{
  if ((names(dev.cur()) != "null device") & (newdev)) 
    suppressWarnings(dev.off())
  lenfile <- nchar(filename)
  if (lenfile > 3) {
    end <- substr(filename, (lenfile - 3), lenfile)
    if (end != ".png") 
      filename <- paste0(filename, ".png")
    png(filename = filename, width = width, height = height, 
        units = "in", res = resol)
  }
  else {
    if (names(dev.cur()) %in% c("null device", "RStudioGD")) 
      dev.new(width = width, height = height, noRStudioGD = TRUE)
  }
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(1, 1), mai=c(0.45,0.45,0.1,0.05),oma=c(0,0, 0, 0))
  par(cex = cex, mgp = c(1.35, 0.35, 0), font.axis = usefont, 
      font = usefont, font.lab = usefont)
  if ((lenfile > 0) & (verbose)) 
    cat("\n Remember to place 'graphics.off()' after plot \n")
  return(invisible(oldpar))
} # end of plotprep



#' @title plotstand - plot optimum model from standLM  vs Geometric mean
#'
#' @description plot optimum model from standLM  vs Geometric mean.
#'   Has options that allow for log-normls P% intervals around each time
#'   period's parameter estimate. Also can rescale the graph to have an average
#'   the same as the geometric mean average of the original time series of data.
#' @param stnd is the list output from standLM
#' @param bars is a logical T or F determining whether to put confidence bounds
#'   around each estimate; defaults to FALSE
#' @param geo is an estimate of the original geometric mean catch rate across
#'   all years. If this is > 0.0 it is used to rescale the graph to the
#'   nominal scale, otherwise the mean of each time-series will be 1.0, which
#'   simplifies visual comparisons. geo defaults to 0.0.
#' @param P is the percentile used for the log-normal confidence bounds, if
#'   they are plotted; defaults to 95.
#' @param catch if it is desired to plot the catch as well as the CPUE
#'   then a vector of catches needs to be input here
#' @param usefont enables the font used in the plot to be modified. Most
#'   publications appear to prefer usefont=1; defaults to 7 - Times bold
#' @param devpar define the plot par values, default=TRUE
#' @return a plot of the model with the smallest AIC (solid line) and the
#'   geometric mean (model 1, always = LnCE ~ Year, the dashed line). 'Year'
#'   could be some other time step.
#' @export plotstand
#' @examples
#' \dontrun{
#' data(sps)
#' splabel = "SpeciesName"
#' labelM <- c("Year","Vessel","Month")
#' sps1 <- makecategorical(labelM[1:3],sps)
#' mods <- makemodels(labelM)
#' out <- standLM(mods,sps1,splabel)
#' plotprep()
#' plotstand(out, bars=TRUE, P=90,geo=100.0,usefont=1)
#' plotstand(out)
#' }
plotstand <- function(stnd,bars=FALSE,geo=0.0,P=95,catch=NA,usefont=7,
                      devpar=TRUE) {
  result <- stnd$Results
  if (geo > 0.0) result <- stnd$Results*geo
  sterr <- stnd$StErr
  whichM <- stnd$WhichM
  optimum <- stnd$Optimum
  years <- rownames(result)
  fishyr <- FALSE
  if (nchar(years[1]) > 4) {
       yrs <- 1:length(years)
       fishyr <- TRUE
    } else {
       yrs <- as.numeric(rownames(result))
    }
  laby <- paste(stnd$Label," CPUE",sep="")
  if (bars) {
    Zmult <- -qnorm((1-(P/100))/2.0)
    lower <- result[,optimum] * exp(-Zmult*sterr[,optimum])
    upper <- result[,optimum] * exp(Zmult*sterr[,optimum])
    ymax <- max(result[,1],result[,optimum],upper,na.rm=TRUE)*1.025
  } else {
    ymax <- max(result[,1],result[,optimum],na.rm=TRUE)*1.025
  }
  if (devpar) {
    if (length(catch) > 1) par(mfrow= c(2,1)) else par(mfrow= c(1,1))
    par(mai=c(0.4,0.5,0,0), oma=c(0,0,0.25,0.25))
    par(cex=0.85, mgp=c(1.5,0.3,0), font.axis=usefont)
  }
  plot(yrs,result[,1],type="l",lty=2,lwd=2,ylim=c(0,ymax),yaxs="i",
       ylab="",xlab="",xaxs="r",panel.first=grid())
  lines(yrs,result[,optimum],lwd=3)
  if (bars) {
    arrows(x0=yrs[-1],y0=lower[-1],x1=yrs[-1],y1=upper[-1],
           length=0.035,angle=90,col=2,lwd=2,code=3)
  }
  title(ylab=list(laby, cex=1.0, col=1, font=usefont))
  if (geo > 0.0) {
    abline(h=geo,col="grey")
  } else {
    abline(h=1.0,col="grey")
  }
  if (length(catch) > 1) {
    if (length(catch) != length(yrs)) {
      stop("input catch data has incorrect number of years")
    }
    ymax <- max(catch,na.rm=TRUE) * 1.05
    plot(yrs,catch,type="b",pch=16,cex=0.8,lwd=2,ylim=c(0,ymax),yaxs="i",ylab="",
         xlab="",xaxs="r")
    grid()
    title(ylab=list("Catch", cex=1.0, col=1, font=usefont))
  }
} # End of plotstand


#' @title plotstandFY - optimum model & Geometric mean for fishing years
#'
#' @description plotstandFY - optimum model & Geometric mean for fishing years.
#'    Standardizations from standLM and standLnCE can be used.
#'    Has options that allow for log-normls P% intervals around each time
#'    period's parameter estimate. Also can rescale the graph to have an average
#'    the same as the geometric mean average of the original time series of data.
#' @param stnd is the list output from standLM
#' @param yrtick defaults to 4 denoting the gap between years on the
#'     x-axis
#' @param bars is a logical T or F determining whether to put confidence bounds
#'    around each estimate; defaults to FALSE
#' @param P is the percentile used for the log-normal confidence bounds, if
#'    they are plotted; defaults to 95.
#' @param lstyr is last years time-series of the optimal model. This can be
#'    plotted on top of the current year's optimum to chack on repeatability.
#' @param alterpar is a logical that, if TRUE, which is the default state, it
#'    will prevent this function from altering the par settings.
#'    If FALSE this will allow for additions to be made to the resulting plot.
#' @return a plot of the model with the smallest AIC (solid line) and the
#'    geometric mean (model 1, always = LnCE ~ Year, the dashed line). 'Year'
#'    could be some other time step.
#' @export plotstandFY
#' @examples
#' \dontrun{
#' data(sps)
#' splabel <- "Species2"
#' labelM <- c("Year","Vessel","Month")
#' sps1 <- makecategorical(labelM,sps)
#' mods <- makemodels(labelM)
#' out <- standLM(mods,sps1,splabel)
#' plotstandFY(out, bars=TRUE, P=90)
#' plotstandFY(out)
#' }
plotstandFY <- function(stnd,yrtick=4,bars=FALSE,P=95,lstyr=NA,
                        alterpar=TRUE) {
  result <- stnd$Results
  sterr <- stnd[[2]]
  whichM <- stnd$WhichM
  optimum<- stnd$Optimum
  years <- rownames(result)
  nyrs <- length(years)
  yrs <- seq(1,nyrs,1)
  laby <- paste(stnd[[6]]," CPUE",sep="")
  opar <- par(no.readonly=TRUE)
  if (bars) {
    Zmult <- -qnorm((1-(P/100))/2.0)
    lower <- result[,optimum] * exp(-Zmult*sterr[,optimum])
    upper <- result[,optimum] * exp(Zmult*sterr[,optimum])
    ymax <- max(result[,1],result[,optimum],upper,na.rm=TRUE)*1.025
    } else {
    ymax <- max(result[,1],result[,optimum],na.rm=TRUE)*1.025
  }
  par(mfrow= c(1,1))
  par(mai=c(0.3,0.5,0,0.1), oma=c(0,0,0.0,0.0))
  par(cex=0.9, mgp=c(1.5,0.3,0), font.axis=7)
  plot(yrs,result[,1],type="l",lty=2,lwd=2,ylim=c(0,ymax),yaxs="i",ylab="",xlab="",
       xaxs="r",xaxt="n")
  lines(yrs,result[,optimum],col=1,lwd=3)
  if (bars) { segments(yrs,lower,yrs,upper,lwd=2)  }
  abline(h=1.0,col="grey")
  if (length(lstyr) > 1) {
    optres <- result[,optimum]  # optimum result
    ratioy <- mean(optres[1:(nyrs-1)])
    lines(yrs[1:nyrs-1],lstyr*ratioy,col=4,lwd=2)
  }
  xloc <-seq(1,nyrs,yrtick)       ## generates tick marks
  axis(side=1,xloc,labels=years[xloc],las=1)
  title(ylab=list(laby, cex=1.0, col=1, font=7))
  par(opar)
} # end of plotstandFY

#' @title qqdiag generates a qqplot with a histogram of residuals
#'
#' @description qqdiag generates a qqplot with a complementary histogram of
#'     the residuals to illustrate the proportion of all residuals along the
#'     qqline. If the qqline deviates from the expected straigt line, which
#'     is red i colour to make for simpler comparisons, then the histogram
#'     enables one to estiamte what proportion of records deviate from
#'     normality. The zero point is identified with a line, as are the
#'     approximate 5% and 95% percentiles. In both cases > 5% is above or
#'     below the blue lines, with < 90% in between depending on the
#'     proportions in each class. To get a more precise estimate use the
#'     invisibly returned histogram values.
#'
#' @param inmodel the optimum model being considered
#' @param plotrug a logical term determining whether a rug is plotted on the
#'     qqplot.
#' @param bins defaults to NA, but can be set to a given series
#' @param hline Include some horizontal lines on the histogram. defaults to 0.
#' @param xinc the increment for tick marks on the xaxis of the histogram
#' @param yinc the increment for tick marks on the y-axis of the histogram
#' @param ylab the y-axis label for the histogram, defaults to 'residuals'
#'
#' @return plots a graph and invisibly returns the output from the histogram
#' @export
#'
#' @examples
#' y <- rep(1:100,2)
#' x <- rnorm(200,mean=10,sd=1)
#' model <- lm(y ~ x)
#' dev.new(width=6,height=3.5,noRStudioGD = TRUE)
#' par(mai=c(0.45,0.45,0.15,0.05),font.axis=7)
#' qqdiag(model,xinc=1,yinc=10,bins=seq(-55,50,2.5))
qqdiag <- function(inmodel,plotrug=FALSE,bins=NA,hline=0.0,
                   xinc=100,yinc=1,ylab="residuals") {
  layout(matrix(c(1,2),ncol=2),widths=c(5,2.5))
  par(mai=c(0.45,0.45,0.15,0.05),oma=c(0.0,0,0.0,0.0))
  par(cex=0.85, mgp=c(1.35,0.35,0), font.axis=7,font=7,font.lab=7)
  resids <- inmodel$residuals
  qs <- quantile(resids,probs=c(0.01,0.05,0.95,0.99))
  if (!is.numeric(bins)) {
    loy <- min(resids); hiy <- max(resids)
    scale <- trunc(100*(hiy - loy)/35) / 100
    loy <- round(loy - (scale/2),2); hiy <- round(hiy + scale,2)
    bins <- seq(loy,hiy,scale)
  } else {
    loy <- min(bins); hiy <- max(bins)
  }
  qqplotout(inmodel,plotrug=plotrug,ylow=loy,yhigh=hiy)
  abline(h=qs,lwd=c(1,2,2,1),col=4)
  outL <- lefthist(resids,bins=bins,hline=0.0,yinc=yinc,xinc=xinc,
                   ylabel=ylab,width=0.9,border=1)
  abline(h=qs,lwd=c(1,2,2,1),col=4)
  ans <- addnorm(outL,resids)
  lines(ans$y,ans$x,lwd=2,col=3)
  return(invisible(outL))
}  # end of qqdiag



#' @title qqplotout plots up a single qqplot for a lm model
#'
#' @description qqplotout generates a single qqplot in isolation from the
#'     plot of a model's diagnostics. It is used with lefthist to
#'     illustrate how well a model matches a normal distribution
#'
#' @param inmodel the optimum model from standLM or dosingle
#' @param title a title for the plot, defaults to 'Normal Q-Q Plot'
#' @param cex the size of the font used, defaults to 0.9
#' @param ylow the lower limit of the residuals
#' @param yhigh he upper limit of the residuals
#' @param plotrug a logical value determinning whether a rug is included
#'
#' @return currently nothing, but it does generate a qqplot to the current
#'     device
#' @export
#'
#' @examples
#' y <- rep(1:100,2)
#' x <- rnorm(200,mean=10,sd=1)
#' model <- lm(y ~ x)
#' dev.new(width=6,height=3.5,noRStudioGD = TRUE)
#' par(mai=c(0.45,0.45,0.15,0.05),font.axis=7)
#' qqplotout(model,ylow=-50,yhigh=50)
qqplotout <- function(inmodel, title="Normal Q-Q Plot", cex=0.9,
                      ylow=-5,yhigh=5,plotrug=FALSE)  {
   resids <- inmodel$residuals
   labs <- cex
   qqnorm(resids, ylab=list("Standardized Residuals", cex=labs, font=7),
          xlab=list("Theoretical Quantiles", cex=labs, font=7),
          main=list(title,cex=labs,font=7),ylim=c(ylow,yhigh))
   qqline(resids, col=2,lwd=2)
   grid()
   if (plotrug) rug(resids)
   abline(v=c(-2.0,2.0),col="grey")
}  # end of qqplotout


#' @title yearBubble: Generates a bubbleplot of x against Year.
#'
#' @description yearBubble: Generates a bubbleplot of x against Year.
#' @param x - a matrix of variable * Year; although it needn't be year
#' @param xlabel - defaults to nothing but allows a custom x-axis label
#' @param ylabel - defaults to nothing but allows a custom y-axis label
#' @param diam - defaults to 0.1, is a scaling factor to adjust bubble size
#' @param vline - defaults to NA but allows vertical ablines to higlight regions
#' @param txt - defaults are lines to vessel numbers, catches, catches, maximumY
#' @param Fyear - defaults to FALSE, if TRUE generates a fishing year x-axis
#' @param xaxis - defaults to TRUE, allows for a custom x-axis if desired by
#'    using something like axis(1,at=years,labels=years).
#' @param yaxis - defaults to TRUE, allows for a custom y-axis if desired by
#'    using something like axis(side=2,at=years,labels=years).
#' @param hline - defaults to FALSE
#' @param nozero - defaults to FALSE, if TRUE replaces all zeros with NA so they
#'    do not appear in the plot
#' @return - invisible, vectors of catch and vessels by year, and radii matrix
#' @export yearBubble
#' @examples
#' \dontrun{
#' data(sps)
#' cbv <- tapply(sps$catch_kg,list(sps$Vessel,sps$Year),sum,na.rm=TRUE)/1000
#' dim(cbv)
#' early <- rowSums(cbv[,1:6],na.rm=TRUE)
#' late <- rowSums(cbv[,7:14],na.rm=TRUE)
#' cbv1 <- cbv[order(late,-early),]
#' plotprep(width=7,height=6)
#' yearBubble(cbv1,ylabel="Catch by Trawl",vline=2006.5,diam=0.2)
#' }
yearBubble <- function(x,xlabel="",ylabel="",diam=0.1,vline=NA,txt=c(4,6,9,11),
                       Fyear=FALSE,xaxis=TRUE,yaxis=TRUE,hline=FALSE,nozero=FALSE) {
  nyrs <- dim(x)[2]
  if (Fyear) {
    tyrs <- colnames(x)  # assumes a yyyy/yyyy format
    if (nchar(tyrs[1]) != 9) warning("Wrong fishing year format for yearBubble \n")
    years <- as.numeric(substr(tyrs,1,4))
  } else { years <- as.numeric(colnames(x)) # assumes columns are years
  }
  nves <- length(rownames(x))
  yvar <- seq(1,nves,1)
  if (nozero) {
    pick <- which(x == 0)
    x[pick] <- NA
  }
  radii <- sqrt(x)
  biggest <- max(radii,na.rm=TRUE)
  catch <- colSums(x,na.rm=TRUE)   # total annual catches
  numves <- apply(x,2,function(x1) length(which(x1 > 0))) # num vess x year
  answer <- list(catch,numves,radii) # generate output
  names(answer) <- c("Catch","Vessels","Radii")
  xspace <- 0.3
  if (nchar(xlabel) > 0) xspace <- 0.45
  par(mfrow= c(1,1))
  par(mai=c(xspace,0.45,0.1,0.1), oma=c(0.0,0.0,0.0,0.0))
  par(cex=0.85, mgp=c(1.5,0.3,0), font.axis=7,font=7)
  xt <- "s"
  yt <- "s"
  if (!xaxis) xt <- "n"
  if (!yaxis) yt <- "n"
  plot(years,years,type="n",xlab="",ylab="",ylim=c(0,(nves+txt[4])),yaxs="r",
       yaxt=yt,xaxt=xt,xaxs="r")
  if (hline) abline(h=yvar,col="grey")
  for (x in 1:nyrs) {
    yr <- years[x]
    odd.even<-x%%2
    if (odd.even == 0) text(yr,nves+txt[3],round(catch[x],0),cex=0.65,font=7)
    else text(yr,nves+txt[2],round(catch[x],0),cex=0.65,font=7)
    text(yr,nves+txt[1],numves[x],cex=0.8,font=7)
    mult <- max(radii[,x],na.rm=TRUE)/biggest
    symbols(rep(yr,nves),yvar,circles=radii[,x],inches=diam*mult,
            bg=rgb(1, 0, 0, 0.5), fg = "black",xlab="",ylab="",add=TRUE)
  }
  
  if (length(vline) > 0) abline(v=c(vline),col="grey")
  title(ylab=list(ylabel, cex=1.0, col=1, font=7))
  return(invisible(answer))
} # end of YearBubble

