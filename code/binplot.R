#!/usr/bin/R
# author: Yen Hoang
# German Rheumatism Research Center (DRFZ) Berlin - 2019



fcs = new.env()


####################################################################
###### helpful print function
printf <- function(...) invisible(print(sprintf(...)))


######## CHOOSE COLORS #############################################
# color         name of color options
fcs$colrs <- function(color) {
  ### light to dark
  
  if (color == "rainbow") {
    cols = c("#001AC0","#001AC0","#00A7D0","#00F0EC","#00F884","#81FF42","#D1FF00","#FFFF00","#FFDB00FF","#FF6D00FF","#FF0000FF","#FF0000FF")
    cols = col2rgb(cols)
  } else if (color == "rainbow.pale") {
    cols = c("#868EC0","#868EC0","#93C5D2","#A8F1EF","#AEF8D6","#DEFFCD","#F2FFB3","#FFFFB3","#FFFACD","#FFD4B3","#FFB3B3")
    cols = col2rgb(cols)
  } else if (color == "green") {
    # green
    cols = c("#F4FBEF","#E0F3D3","#C3E1AD","#ACD88D","#93CC6B","#7DC24C","#71B441","#63AB30","#56972B","#4F8726")
    cols = col2rgb(cols)
  } 

  cols
}


######## ADD CUTOFF LINES ##########################################
##### @param
# cutoffs       cutoff vector (x,y)
fcs$addProdline <- function (cutoffs=c(0,0)) {
  this=fcs
  # add production cutoff line if provided
  if ( cutoffs[1] > 0 ) {
    abline(v=cutoffs[1],col="darkgrey")
  }
  if ( cutoffs[2] > 0 ) {
    abline(h=cutoffs[2],col="darkgrey")
  }
}


######## CREATE BINPLOT TABLE this$tab3D w/ OR w/o PLOTTING ##############
##### @param
# data   			      asinh transformed data file
# feat.X            name of feature X
# feat.Y            name of feature Y
# feat.Z1           name of feature Z1
# calc              calculation method, options="density", "freq", "MFI", "MFI+"
# cutoffs 			    cutoffs for X,Y,Z1, default=c(0,0,0)
# col             	color palette for calculation of Z1, default="rainbow"
# binsize           bin size, default = 0.2
# mincells          minimum amount of cells in a bin, default = 10
# plot.range        plot range, first x axis, second y axis, default= c(3,12,3,12)
# manual.Z1range	  set only if a global range is preferred for sample comparison
# plotting 			    plot binplot with TRUE OR get binplot table with FALSE, default=TRUE
fcs$binplot_table <- function(
  data, 
  feat.X, 
  feat.Y, 
  feat.Z1, 
  calc,
  cutoffs=c(0,0,0), 
  col="rainbow", 
  binsize=0.2, 
  mincells=10, 
  plot.range=c(3,12,3,12),
  manual.Z1range=NA,
  plotting=TRUE)
{ 
  this=fcs 
  
  ### remove all signs and write anything with capitals
  feat.X.clean = gsub("[^[:alnum:]]","",toupper(feat.X))
  feat.Y.clean = gsub("[^[:alnum:]]","",toupper(feat.Y))

  
  ### axes range
  xmin.val = plot.range[1]
  xmax.val = plot.range[2]
  ymin.val = plot.range[3]
  ymax.val = plot.range[4]
  
  ### change legend numbers in red if range is small
  col.legend = "black"
  
  features = colnames(data) 
  
  ### remove all signs and write anything with capitals
  features.clean = gsub("[^[:alnum:]]","",make.unique(unlist(lapply(features,function(x) {
    len = length(strsplit(x,"[.]")[[1]])
    y = toupper(strsplit(x,"[.]")[[1]])
    paste(y, collapse="")
  }))))
  
  ### get indices for feature X and Y
  idx.X = which(features.clean == feat.X.clean)
  idx.Y = which(features.clean == feat.Y.clean)
  
  
  if ( is.na(idx.X) | is.na(idx.Y)) stop(sprintf("Either feat.X (%s) or feat.Y (%s) not in data",feat.X, feat.Y))
  
  
  if (plotting) {
    ### initiate plot frame
    plot(1, type='n', frame.plot=FALSE,
         xlim=c(xmin.val,xmax.val+10*binsize),axes=FALSE,
         ylim=c(ymin.val-2.5*binsize,ymax.val+5*binsize),
         xlab="",ylab="")
    box(lwd=0.5,col="darkgrey")
    ### draw axis label on the bottom
    asinh.scale = c(asinh(-1000),asinh(-100),asinh(-10),asinh(0),asinh(10),asinh(100),asinh(1000),asinh(10000),asinh(100000),asinh(1000000),asinh(10000000),asinh(100000000))
    axis(side=1, at=asinh.scale,labels=c("1e-3","1e-2","1e-1","0","1e1","1e2","1e3","1e4","1e5","1e6","1e7","1e8"),
         las=1,cex.axis=0.7,col="darkgrey")
    ### draw axis label on the left
    axis(side=2, at=asinh.scale,labels=c("1e-3","1e-2","1e-1","0","1e1","1e2","1e3","1e4","1e5","1e6","1e7","1e8"),
         las=3,cex.axis=0.7,col="darkgrey")
    
    ### PERCENTAGE COLORS
    quadrants.color="black"
    prodcells.color="red"
    prodpluscells.color="chartreuse4"

    ### calculate quadrant percentages (black)
    ncells.total = nrow(data)
    
    # q1 = bottom left
    tdata.q1 = data[which( data[,idx.X]<cutoffs[1] &  data[,idx.Y]<cutoffs[2] ),]
    # q2 = bottom right
    tdata.q2 = data[which( data[,idx.X]>=cutoffs[1] &  data[,idx.Y]<cutoffs[2] ),]
    # q3 = top right
    tdata.q3 = data[which( data[,idx.X]>=cutoffs[1] &  data[,idx.Y]>=cutoffs[2] ),]
    # q4 = top left
    tdata.q4 = data[which( data[,idx.X]<cutoffs[1] &  data[,idx.Y]>=cutoffs[2] ),]
    
    this$tdata.q1=tdata.q1
    ### % cells in quadrant to total cells (black)
    q1.total = abs(100 * length( tdata.q1 ) / ncells.total)
    q2.total = abs(100 * length( tdata.q2 ) / ncells.total)
    q3.total = abs(100 * length( tdata.q3 ) / ncells.total)
    q4.total = abs(100 - q1.total - q2.total - q3.total)
    # printf("q1.total=%0.1f, idx.Z1=%s, col(tdata.q1)[idx.Z1]=%s",q1.total,idx.Z1,colnames(tdata.q1)[idx.Z1])
    
    
    if ( calc == "MFI+" ) {
      quadrants.color = "blue"
    }
    
    ####### START PLOT PERCENTAGES
    ### quadrant bottom left
    text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.02*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q1.total),col=quadrants.color,cex=0.9,pos=4,xpd=TRUE)
    ### quadrant bottom right
    text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.02*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q2.total),col=quadrants.color,cex=0.9,pos=2,xpd=TRUE)
    ### quadrant top right
    text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q3.total),col=quadrants.color,cex=0.9,pos=2,xpd=TRUE)
    ### quadrant top left
    text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q4.total),col=quadrants.color,cex=0.9,pos=4,xpd=TRUE)
    ####### STOP PLOT PERCENTAGES
  }
  
    
    
  ### get indices for feature Z1
  if ( calc != "density" ) {
    feat.Z1.clean = gsub("[^[:alnum:]]","",toupper(feat.Z1))
    idx.Z1 = which(features.clean == feat.Z1.clean)
    
    ### set negative values of Z1 (colnum=3) to zero
    data[which(data[,idx.Z1]<0),idx.Z1] = 0
    
    if ( is.na(idx.Z1) ) stop(sprintf("feat.Z1 (%s) not in data.",feat.Z1))
    
    
    if (plotting) {
      # number of cells double producers in quadrants
      q1.prodcells.num = nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[3]),])
      q2.prodcells.num = nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[3]),])
      q3.prodcells.num = nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[3]),])
      q4.prodcells.num = nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[3]),])
      
      # % cells double producers to quadrant cells
      q1.prodcells = round( (100 * q1.prodcells.num / nrow(tdata.q1) ), 2)
      q2.prodcells = round( (100 * q2.prodcells.num / nrow(tdata.q2) ), 2)
      q3.prodcells = round( (100 * q3.prodcells.num / nrow(tdata.q3) ), 2)
      q4.prodcells = round( (100 * q4.prodcells.num / nrow(tdata.q4) ), 2)
      
      # % cells double producers to total cells
      q1.prodcells.total = round( 100 * nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[3]),]) / ncells.total,2)
      q2.prodcells.total = round( 100 * nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[3]),]) / ncells.total,2)
      q3.prodcells.total = round( 100 * nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[3]),]) / ncells.total,2)
      q4.prodcells.total = round( 100 * nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[3]),]) / ncells.total,2)
      
      ####### START PLOT PERCENTAGES
      ### quadrant bottom left
      text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.06*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q1.prodcells),col=prodcells.color,cex=0.9,pos=4,xpd=TRUE)
      text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.10*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q1.prodcells.total),col=prodpluscells.color,cex=0.9,pos=4,xpd=TRUE)
      
      ### quadrant bottom right
      text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.06*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q2.prodcells),col=prodcells.color,cex=0.9,pos=2,xpd=TRUE)
      text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.10*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q2.prodcells.total),col=prodpluscells.color,cex=0.9,pos=2,xpd=TRUE)
      
      ### quadrant top right
      text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.07*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q3.prodcells),col=prodcells.color,cex=0.9,pos=2,xpd=TRUE)
      text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.11*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q3.prodcells.total),col=prodpluscells.color,cex=0.9,pos=2,xpd=TRUE)
      
      ### quadrant top left
      text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.07*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q4.prodcells),col=prodcells.color,cex=0.9,pos=4,xpd=TRUE)
      text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.11*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.1f%%",q4.prodcells.total),col=prodpluscells.color,cex=0.9,pos=4,xpd=TRUE)
      ####### STOP PLOT PERCENTAGES
    }
  }
  
  ### construct bin table with number of cells per bin
  fX = cut(data[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  fY = cut(data[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  fXY=as.factor(paste(fX,fY))
  
  tab = table(fX,fY)
  colnames(tab)=seq(ymin.val,ymax.val-binsize,by=binsize)
  rownames(tab)=seq(xmin.val,xmax.val-binsize,by=binsize)
  tab.origin = tab
  
  if (calc == "MFI+") {
    ### cut data if calc method is MFI+
    # filter data with Z1 positive
    data = data[which(data[,idx.Z1]>cutoffs[3]),]
    
    ### construct new bin table with number of cells per bin
    fX = cut(data[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
    fY = cut(data[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
    fXY=as.factor(paste(fX,fY))
    
    tab = table(fX,fY)
    colnames(tab)=seq(ymin.val,ymax.val-binsize,by=binsize)
    rownames(tab)=seq(xmin.val,xmax.val-binsize,by=binsize)
  }
  
  ########## START CALCULATE BIN PROPERTIES
  if ( calc == "density") {
    ########## CALCULATE NUMBER OF CELLS
    my.calc=aggregate(data[,idx.X],by=list(fXY),length)
    idx.len = which(my.calc$x >= mincells)
    
    min.range=floor(min(my.calc[idx.len,'x'])*10)/10
    max.range=max(tab)
    
    # get steps for Z1
    #printf("step.Z1=%s",round(diff(range(max.range,min.range))/10,2))
    step.Z1=round(diff(range(max.range,min.range))/10,2)
    steps.Z1=seq(min.range,max.range,by=step.Z1)
    
    ### legend steps
    legend.steps = round(steps.Z1,-2)
    
    # bin color factor Z1
    my.calc.fac.Z1=cut(my.calc$x,breaks=steps.Z1,labels=1:10,include.lowest=TRUE)
    names(my.calc.fac.Z1) = my.calc$x
  } else if ( grepl("MFI",calc) ) {
    ########## CALCULATE MFI
    my.lengths = aggregate(data[,idx.Z1],by=list(fXY),length)
    idx.len = which(my.lengths$x >= mincells)
    
    ### start with MFI calculation of Z1
    my.calc = aggregate(data[,idx.Z1],by=list(fXY),mean)
    
    if (is.na(manual.Z1range)) {
      # dynamic range
      min.range.Z1 = floor(min(my.calc[idx.len,'x'])*10)/10
      max.range.Z1 = ceiling(max(my.calc[idx.len,'x'])*10)/10
    } else{
      # manual range
      min.range.Z1 = manual.Z1range[1]
      max.range.Z1 = manual.Z1range[2]
    }
    
    ### get steps for Z1
    #printf("step.Z1=%s",round(diff(range(max.range.Z1,min.range.Z1))/10,2) )
    step.Z1=round(diff(range(max.range.Z1,min.range.Z1))/10,2) 
    steps.Z1=seq(min.range.Z1,max.range.Z1,by=step.Z1)
    
    ### legend steps
    legend.steps = round(steps.Z1,1)
    
    ### bin color factor Z1
    my.calc.fac.Z1=cut(my.calc$x,breaks=steps.Z1,labels=1:10,include.lowest=TRUE)
    names(my.calc.fac.Z1) = my.calc$x
    
  } else if ( calc == "freq" ) {
    ########## CALCULATE FREQUENCIES
    ### frequency of feature Z1
    my.calc = aggregate(data[,idx.Z1],by=list(fXY),function(x) {
      y= round( 100 * length(which(x >= cutoffs[1])) / length(x))
      return(y)
    })
    my.lengths = aggregate(data[,idx.Z1], by=list(fXY), length)
    
    ### legend steps
    step.Z1 = 10
    steps.Z1=seq(0,100,by=step.Z1)
    legend.steps = round(steps.Z1)
    
    # bin color factor Z1
    my.calc.fac.Z1 = cut(my.calc$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
    names(my.calc.fac.Z1) = my.calc$x
  }
  ### combine all values in one table
  my.calc = cbind(my.calc,fac.Z1=as.numeric(my.calc.fac.Z1))
  
  ########## DONE CALCULATE BIN PROPERTIES
  
  # help vector
  brackets.open = c("(","[")
  
  if (plotting) {
    cols = fcs$colrs(col)
    ########## PLOT BINPLOT 
    for (x in rownames(tab)) {
      for (y in colnames(tab)) {
        if (tab[x,y]>=mincells) {
          brackets.idx.x = brackets.idx.y = 1
          if (x==0) brackets.idx.x = 2
          if (y==0) brackets.idx.y = 2
          
          fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+binsize,'] ',
                                 brackets.open[brackets.idx.y],y,',',as.numeric(y)+binsize,']',sep=''))
          
          idx = which(as.character(fact)==as.character(my.calc$Group.1))
          if (length(idx) != 0) {
            rect(x,y,as.numeric(x)+binsize,as.numeric(y)+binsize,
                 col=eval(parse(text=paste0("rgb(",paste0(cols[,my.calc[idx,'fac.Z1']+1],collapse=","),",maxColorValue=255)"))),
                 border=NA)
            }
          ### if calc=="MFI+"
        } else if ( calc == "MFI+" & tab.origin[x,y] >= mincells) {
          rect(x,y,as.numeric(x)+binsize,as.numeric(y)+binsize,
               col="gray",border=NA)
        }
      }
    }
    
    this$addProdline(cutoffs[c(1,2)])
    
    ### x axis label
    text(x = 0.5*(par()$usr[1]+par()$usr[2]),
         y = par()$usr[3] - 1,
         label = feat.X,
         xpd = TRUE
    )
    ### y axis label
    text(x = par()$usr[1] - 0.8,
         y = 0.5*(par()$usr[3]+par()$usr[4]),
         label = feat.Y,
         xpd = TRUE,
         srt = 90
    )
    
    
    ########## PLOT LEGEND START
    
    ### start position for x and y
    legend.pos.x = par()$usr[2]-0.16*(par()$usr[2]-par()$usr[1])
    legend.pos.y = par()$usr[3]+0.14*(par()$usr[4]-par()$usr[3])
    
    ### from bottom to top
    for (i in 1:10) {
      # display legend rectangles
      rect(xleft = legend.pos.x,
           ybottom = legend.pos.y+0.03*(par()$usr[4]-par()$usr[3])*i,
           xright = legend.pos.x+(0.02*(par()$usr[2]-par()$usr[1])),
           ytop = legend.pos.y+0.03*(par()$usr[4]-par()$usr[3])*(i+1),
           col=eval(parse(text=paste0("rgb(",paste0(cols[,i+1],collapse=","),",maxColorValue=255)"))),
           border=NA, xpd=TRUE
      )
      
      ### DISPLAY LEGEND NUMBERS
      if (i==1 | i==6) {
        # display lowest and middle legend number
        if (grepl("MFI",calc)) {
          display.label = sprintf("%.1f",(as.numeric(legend.steps[i])))
        } else if (calc=="density") {
          display.label = sprintf("%s",round(as.numeric(legend.steps[i]),-2))
        } else {
          display.label = sprintf("%.0f",as.numeric(legend.steps[i]))
        }
        text(x=legend.pos.x+(0.007*(par()$usr[2]-par()$usr[1])),
             y=legend.pos.y+0.03*(par()$usr[4]-par()$usr[3])*i,
             label=display.label,col=col.legend,cex=0.8,pos=4, xpd=TRUE)
      } else if (i==3 | i==8) {
        # display second lowest and second highest legend number
        if (grepl("MFI",calc)) {
          display.label = sprintf("%.1f",(as.numeric(legend.steps[i]) + (step.Z1/2)))
        } else if (calc=="density") {
          display.label = sprintf("%.0f",round(as.numeric(legend.steps[i]) + (step.Z1/2)),-2)
        } else {
          display.label = sprintf("%.0f",(as.numeric(legend.steps[i]) + (step.Z1/2)))
        }
        text(x=legend.pos.x+(0.007*(par()$usr[2]-par()$usr[1])),
             y=legend.pos.y+0.03*(par()$usr[4]-par()$usr[3])*(i+0.5),
             label=display.label,col=col.legend,cex=0.8,pos=4, xpd=TRUE)
      } else if (i==10) {
        # display highest legend number
        display.label = sprintf("%.1f",as.numeric(legend.steps[i+1]))
        display.label = sprintf("%s",as.numeric(legend.steps[i+1]))
        text(x=legend.pos.x+(0.007*(par()$usr[2]-par()$usr[1])),
             y=legend.pos.y+0.03*(par()$usr[4]-par()$usr[3])*(i+1),
             label=display.label,col=col.legend,cex=0.8,pos=4, xpd=TRUE)
      }
    }
    
    ### display legend title Z1
    if (calc=="density") {
      legend.title = "# cells"
    } else {
      legend.title = feat.Z1
    }
    
    text(x = par()$usr[2]-0.20*(par()$usr[2]-par()$usr[1]),
         y = legend.pos.y+0.03*(par()$usr[4]-par()$usr[3])*(i+2.8),
         label = legend.title,
         cex=0.9, pos=4, xpd=TRUE
    )
    ########## PLOT LEGEND END
    
    
  } else {
    ########## OR CREATE BINPLOT MATRIX ONLY
    for (x in rownames(tab)) {
      for (y in colnames(tab)) {
        if ( tab[x,y]>=mincells ) {
          brackets.idx.x = brackets.idx.y = 1
          if (x==0) brackets.idx.x = 2
          if (y==0) brackets.idx.y = 2
          
          fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+binsize,'] ',
                                 brackets.open[brackets.idx.y],y,',',as.numeric(y)+binsize,']',sep=''))
          
          idx=which(as.character(fact)==as.character(my.calc$Group.1))
          
          tab[x,y] = my.calc[idx,'x']
          #} else if ( tab[x,y]>0 & calc!="freq") {
          #    tab[x,y] = min.range.Z1
        } else {
          tab[x,y] = NA
        }
      }
    }
    
    #this$tab3D = round(tab,2)
    round(tab,2)
  }
}

attach(fcs)