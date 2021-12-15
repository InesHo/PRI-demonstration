#!/usr/bin/R
# author: Yen Hoang
# German Rheumatism Research Center (DRFZ) Berlin - 2019

####################################################################
###### helpful print function
printf <- function(...) invisible(print(sprintf(...)))


######## CHOOSE COLORS #############################################
# color         name of color options
fcs$colors <- function(color) {
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



######## CREATE BINPLOT BIN SCAFFOLD #######
### NO FREQUENCY CALCULATION HERE
##### @param
# data              asinh transformed data file
# feat.X            name of feature X
# feat.Y            name of feature Y
# binsize           bin size, default=0.2
# mincells          minimum amount of cells in a bin, default=10
# plot.range        plot range, first x axis, second y axis, default= c(3,12,3,12)
fcs$binplot_construct <- function(
  data, 
  feat.X, 
  feat.Y, 
  binsize=0.2, 
  mincells=10, 
  plot.range=c(3,12,3,12)) 
{ 
  this=fcs 
  
  ### remove all signs and write anything with capitals
  feat.X.clean = gsub("[^[:alnum:]]","",toupper(feat.X))
  feat.Y.clean = gsub("[^[:alnum:]]","",toupper(feat.Y))
  
  xmin.val = plot.range[1]
  xmax.val = plot.range[2]
  ymin.val = plot.range[3]
  ymax.val = plot.range[4]
  
  features = colnames(data) 
  
  ### remove all signs and write anything with capitals
  features.clean = gsub("[^[:alnum:]]","",make.unique(unlist(lapply(features,function(x) {
    len = length(strsplit(x,"[.]")[[1]])
    y = toupper(strsplit(x,"[.]")[[1]][1])
    paste(y, collapse=".")
  }))))
  
  idx.X = which(features.clean == feat.X.clean)
  idx.Y = which(features.clean == feat.Y.clean)
  
  if ( is.na(idx.X) | is.na(idx.Y) ) stop(sprintf("Either feat.X (%s) or feat.Y (%s) not in data.",feat.X, feat.Y))
  
  fX = cut(data[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  fY = cut(data[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  
  ### if first for loop run
  ### construct bin table with number of cells per bin
  tab = table(fX,fY)
  colnames(tab)=seq(ymin.val,ymax.val-binsize,by=binsize)
  rownames(tab)=seq(xmin.val,xmax.val-binsize,by=binsize)
  
  fXY=as.factor(paste(fX,fY))
  
  ### start plot frame
  plot(1,type='n',frame.plot=FALSE,xlim=c(xmin.val,xmax.val+10*binsize),axes=FALSE
       ,ylim=c(ymin.val-2.5*binsize,ymax.val+5*binsize),xlab=NA,ylab=NA,cex.lab=1,cex.axis=1,mgp=c(1.7,0.4,0))
  box(lwd=0.8,col="darkgrey")
  
  ### draw an axis on the bottom
  ### draw an axis on the left
  asinh.scale = c(asinh(-1000),asinh(-100),asinh(-10),asinh(0),asinh(10),asinh(100),asinh(1000),asinh(10000),asinh(100000),asinh(1000000),asinh(10000000),asinh(100000000))
  axis(side=1, at=asinh.scale,labels=c("1e-3","1e-2","1e-1","0","1e1","1e2","1e3","1e4","1e5","1e6","1e7","1e8"),
       las=1,cex.axis=0.7,col="darkgrey")
  axis(side=2, at=asinh.scale,labels=c("1e-3","1e-2","1e-1","0","1e1","1e2","1e3","1e4","1e5","1e6","1e7","1e8"),
       las=3,cex.axis=0.7,col="darkgrey")
  
  ### add grid
  xgrid.steps=seq((xmin.val),(xmax.val),by=2)
  ygrid.steps=seq((ymin.val),(ymax.val),by=2)
  abline(h=ygrid.steps,v=xgrid.steps,col="grey",lty=3)
  
  ##### plot bin construct in grey first
  for (x in rownames(tab)) {
    for (y in colnames(tab)) {
      if ( tab[x,y]>=mincells ) {
        rect(x,y,as.numeric(x)+binsize,as.numeric(y)+binsize,col="lightgrey",border=NA)
      }
    }
  }
}


######## CREATE BINPLOT BINS WITH FREQUENCY OF DOUBLE POSITIVES ####
##### @param
# data              asinh transformed data file
# feat.X            name of feature X
# feat.Y            name of feature Y
# feat.Z1           name of feature Z1
# feat.Z2           name of feature Z2
# cutoffs 			    cutoffs for X,Y,Z1, default=c(0,0,0,0)
# col               color palette for frequency of Z1+/Z2+, default="green"
# binsize           bin size, default=0.2
# mincells          minimum amount of cells in a bin, default=10
# maxfreq           maximum frequency for scaling the color, default=100
# plot.range        plot range, first x axis, second y axis, default= c(3,12,3,12)  
# sepa="."          separater for the feature names 
fcs$binplot_freq_doublepos <- function(
  data, 
  feat.X, 
  feat.Y, 
  feat.Z1, 
  feat.Z2, 
  cutoffs=c(0,0,0,0),
  col="green", 
  binsize=0.2, 
  mincells=10, 
  maxfreq=100, 
  plot.range=c(3,12,3,12)) 
{ 
  this=fcs 
  
  ### axes range
  xmin.val = plot.range[1]
  xmax.val = plot.range[2]
  ymin.val = plot.range[3]
  ymax.val = plot.range[4]
  
  ### remove all signs and write anything with capitals
  feat.X.clean = gsub("[^[:alnum:]]","",toupper(feat.X))
  feat.Y.clean = gsub("[^[:alnum:]]","",toupper(feat.Y))
  feat.Z1.clean = gsub("[^[:alnum:]]","",toupper(feat.Z1))
  feat.Z2.clean = gsub("[^[:alnum:]]","",toupper(feat.Z2))
  
  features = colnames(data) 
  
  ### remove all signs and write anything with capitals
  features.clean = gsub("[^[:alnum:]]","",make.unique(unlist(lapply(features,function(x) {
    len = length(strsplit(x,"[.]")[[1]])
    y = toupper(strsplit(x,"[.]")[[1]][1])
    paste(y, collapse=".")
  }))))
  
  ### get indices for feature X, Y, Z1, Z2
  idx.X = which(features.clean == feat.X.clean)
  idx.Y = which(features.clean == feat.Y.clean)
  idx.Z1 = which(features.clean == feat.Z1.clean)
  idx.Z2 = which(features.clean == feat.Z2.clean)
  
  if ( is.na(idx.X) | is.na(idx.Y) ) {
    stop(sprintf("Either feat.X (%s) or feat.Y (%s) not in data.",feat.X, feat.Y))
  }
  if ( is.na(idx.Z1) | is.na(idx.Z2) ) {
    stop(sprintf("Either feat.Z1 (%s) or feat.Z2 (%s) not in data.",feat.Z1, feat.Z2))
  }
  
  ### calculate quadrant percentages
  ncells.total = nrow(data)
  
  # cut all cells which are not producing cells
  # q1 = bottom left
  tdata.q1 = data[which( data[,idx.X]<cutoffs[1] &  data[,idx.Y]<cutoffs[2] ),]
  # q2 = bottom right
  tdata.q2 = data[which( data[,idx.X]>=cutoffs[1] &  data[,idx.Y]<cutoffs[2] ),]
  # q3 = top right
  tdata.q3 = data[which( data[,idx.X]>=cutoffs[1] &  data[,idx.Y]>=cutoffs[2] ),]
  # q4 = top left
  tdata.q4 = data[which( data[,idx.X]<cutoffs[1] &  data[,idx.Y]>=cutoffs[2] ),]
  
  # number of cells double producers in quadrants
  q1.prodcells.num = nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[3] & tdata.q1[,idx.Z2]>=cutoffs[4]),])
  q2.prodcells.num = nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[3] & tdata.q2[,idx.Z2]>=cutoffs[4]),])
  q3.prodcells.num = nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[3] & tdata.q3[,idx.Z2]>=cutoffs[4]),])
  q4.prodcells.num = nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[3] & tdata.q4[,idx.Z2]>=cutoffs[4]),])
  
  # % cells double producers to quadrant cells
  q1.prodcells = round( (100 * q1.prodcells.num / nrow(tdata.q1) ), 2)
  q2.prodcells = round( (100 * q2.prodcells.num / nrow(tdata.q2) ), 2)
  q3.prodcells = round( (100 * q3.prodcells.num / nrow(tdata.q3) ), 2)
  q4.prodcells = round( (100 * q4.prodcells.num / nrow(tdata.q4) ), 2)
  
  # % cells double producers to total cells
  q1.prodcells.total = round( 100 * nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[3] & tdata.q1[,idx.Z2]>=cutoffs[4]),]) / ncells.total,2)
  q2.prodcells.total = round( 100 * nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[3] & tdata.q2[,idx.Z2]>=cutoffs[4]),]) / ncells.total,2)
  q3.prodcells.total = round( 100 * nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[3] & tdata.q3[,idx.Z2]>=cutoffs[4]),]) / ncells.total,2)
  q4.prodcells.total = round( 100 * nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[3] & tdata.q4[,idx.Z2]>=cutoffs[4]),]) / ncells.total,2)
  
  # printf("cells quadrant:                %s %s %s %s",nrow(tdata.q1),nrow(tdata.q2),nrow(tdata.q3),nrow(tdata.q4))
  # printf("double producing cells in quadrant: %s %s %s %s",q1.prodcells.num,q2.prodcells.num,q3.prodcells.num,q4.prodcells.num)
  # printf("%% double producing cells to cells in quadrant (red): %s %s %s %s",q1.prodcells,q2.prodcells,q3.prodcells,q4.prodcells)
  # printf("%% double producing cells to total cells (green): %s %s %s %s",q1.prodcells.total,q2.prodcells.total,q3.prodcells.total,q4.prodcells.total)
  
  ####### START PLOT PERCENTAGES
  prodcells.color="red"
  prodpluscells.color="chartreuse4"
  ### quadrant bottom left
  text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q1.prodcells),col=prodcells.color,cex=0.9,pos=4,xpd=TRUE)
  text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.09*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q1.prodcells.total),col=prodpluscells.color,cex=0.9,pos=4,xpd=TRUE)
  
  ### quadrant bottom right
  text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q2.prodcells),col=prodcells.color,cex=0.9,pos=2,xpd=TRUE)
  text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.09*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q2.prodcells.total),col=prodpluscells.color,cex=0.9,pos=2,xpd=TRUE)
  
  ### quadrant top right### calculate quadrant percentages
  ncells.total = nrow(data)
  
  # cut all cells which are not producing cells
  # q1 = bottom left
  tdata.q1 = data[which( data[,idx.X]<cutoffs[1] &  data[,idx.Y]<cutoffs[2] ),]
  # q2 = bottom right
  tdata.q2 = data[which( data[,idx.X]>=cutoffs[1] &  data[,idx.Y]<cutoffs[2] ),]
  # q3 = top right
  tdata.q3 = data[which( data[,idx.X]>=cutoffs[1] &  data[,idx.Y]>=cutoffs[2] ),]
  # q4 = top left
  tdata.q4 = data[which( data[,idx.X]<cutoffs[1] &  data[,idx.Y]>=cutoffs[2] ),]
  
  # number of cells double producers in quadrants
  q1.prodcells.num = nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[3] & tdata.q1[,idx.Z2]>=cutoffs[4]),])
  q2.prodcells.num = nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[3] & tdata.q2[,idx.Z2]>=cutoffs[4]),])
  q3.prodcells.num = nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[3] & tdata.q3[,idx.Z2]>=cutoffs[4]),])
  q4.prodcells.num = nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[3] & tdata.q4[,idx.Z2]>=cutoffs[4]),])
  
  # % cells double producers to quadrant cells
  q1.prodcells = round( (100 * q1.prodcells.num / nrow(tdata.q1) ), 2)
  q2.prodcells = round( (100 * q2.prodcells.num / nrow(tdata.q2) ), 2)
  q3.prodcells = round( (100 * q3.prodcells.num / nrow(tdata.q3) ), 2)
  q4.prodcells = round( (100 * q4.prodcells.num / nrow(tdata.q4) ), 2)
  
  # % cells double producers to total cells
  q1.prodcells.total = round( 100 * nrow(tdata.q1[which(tdata.q1[,idx.Z1]>=cutoffs[3] & tdata.q1[,idx.Z2]>=cutoffs[4]),]) / ncells.total,2)
  q2.prodcells.total = round( 100 * nrow(tdata.q2[which(tdata.q2[,idx.Z1]>=cutoffs[3] & tdata.q2[,idx.Z2]>=cutoffs[4]),]) / ncells.total,2)
  q3.prodcells.total = round( 100 * nrow(tdata.q3[which(tdata.q3[,idx.Z1]>=cutoffs[3] & tdata.q3[,idx.Z2]>=cutoffs[4]),]) / ncells.total,2)
  q4.prodcells.total = round( 100 * nrow(tdata.q4[which(tdata.q4[,idx.Z1]>=cutoffs[3] & tdata.q4[,idx.Z2]>=cutoffs[4]),]) / ncells.total,2)
  
  # printf("cells quadrant:                %s %s %s %s",nrow(tdata.q1),nrow(tdata.q2),nrow(tdata.q3),nrow(tdata.q4))
  # printf("double producing cells in quadrant: %s %s %s %s",q1.prodcells.num,q2.prodcells.num,q3.prodcells.num,q4.prodcells.num)
  # printf("%% double producing cells to cells in quadrant (red): %s %s %s %s",q1.prodcells,q2.prodcells,q3.prodcells,q4.prodcells)
  # printf("%% double producing cells to total cells (green): %s %s %s %s",q1.prodcells.total,q2.prodcells.total,q3.prodcells.total,q4.prodcells.total)
  
  ####### START PLOT PERCENTAGES
  prodcells.color="red"
  prodpluscells.color="chartreuse4"
  ### quadrant bottom left
  text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q1.prodcells),col=prodcells.color,cex=0.9,pos=4,xpd=TRUE)
  text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.09*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q1.prodcells.total),col=prodpluscells.color,cex=0.9,pos=4,xpd=TRUE)
  
  ### quadrant bottom right
  text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.03*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q2.prodcells),col=prodcells.color,cex=0.9,pos=2,xpd=TRUE)
  text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[3]+0.09*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q2.prodcells.total),col=prodpluscells.color,cex=0.9,pos=2,xpd=TRUE)
  
  ### quadrant top right
  text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q3.prodcells),col=prodcells.color,cex=0.9,pos=2,xpd=TRUE)
  text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.10*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q3.prodcells.total),col=prodpluscells.color,cex=0.9,pos=2,xpd=TRUE)
  
  ### quadrant top left
  text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q4.prodcells),col=prodcells.color,cex=0.9,pos=4,xpd=TRUE)
  text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.10*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q4.prodcells.total),col=prodpluscells.color,cex=0.9,pos=4,xpd=TRUE)
  ####### STOP PLOT PERCENTAGES
  
  text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q3.prodcells),col=prodcells.color,cex=0.9,pos=2,xpd=TRUE)
  text(par()$usr[2]+0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.10*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q3.prodcells.total),col=prodpluscells.color,cex=0.9,pos=2,xpd=TRUE)
  
  ### quadrant top left
  text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.04*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q4.prodcells),col=prodcells.color,cex=0.9,pos=4,xpd=TRUE)
  text(par()$usr[1]-0.026*(par()$usr[2]-par()$usr[1]),par()$usr[4]-0.10*(par()$usr[4]-par()$usr[3]),label=sprintf("%0.2f%%",q4.prodcells.total),col=prodpluscells.color,cex=0.9,pos=4,xpd=TRUE)
  ####### STOP PLOT PERCENTAGES
  
  ### construct bin table with number of cells per bin
  fX = cut(data[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  fY = cut(data[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  fXY=as.factor(paste(fX,fY))
  
  tab = table(fX,fY)
  colnames(tab)=seq(ymin.val,ymax.val-binsize,by=binsize)
  rownames(tab)=seq(xmin.val,xmax.val-binsize,by=binsize)
  tab3D = tab3D_z2 = tab
  
  ########## START CALCULATE BIN PROPERTIES
  ### frequency of feature Z1
  my.calc = aggregate(data[,idx.Z1], 
                      by=list(fXY),
                      function(x) {
                        y = round( 100 * length(which(x >= cutoffs[3])) / length(x))
                        return(y)
                      })
  # bin color factor Z1
  my.calc.fac.Z1 = cut(my.calc$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
  names(my.calc.fac.Z1) = my.calc$x
  
  ### frequency of feature Z2        
  freq.Z2 = aggregate(data[,idx.Z2],
                      by=list(fXY),
                      function(x) {
                        y = round( 100 * length(which(x >= cutoffs[4])) / length(x))
                        return(y)
                      })
  # bin color factor Z2
  my.calc.fac.Z2 = cut(freq.Z2$x,breaks=seq(0,100,by=10),labels=1:10,include.lowest=TRUE)
  names(my.calc.fac.Z2) = freq.Z2$x
  
  ### get only data table where Z1 and Z2 is produced
  data.double = data[ which( (data[,idx.Z1]>=cutoffs[3]) & (data[,idx.Z2]>=cutoffs[4]) ),]
  ### construct bin table with number of cells per bin
  fX.double = cut(data.double[,idx.X],breaks=seq(xmin.val,xmax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  fY.double = cut(data.double[,idx.Y],breaks=seq(ymin.val,ymax.val,by=binsize),include.lowest=TRUE,dig.lab=5)
  fXY.double = as.factor(paste(fX.double,fY.double))
  
  length.double = aggregate(data.double[,1],by=list(fXY.double),length)
  rownames(length.double)=length.double$Group.1
  
  length.all = aggregate(data[,idx.Z1], by=list(fXY), length)
  rownames(length.all)=length.all$Group.1
  
  freq.double = merge(length.all,length.double,by="row.names",all.x=TRUE)
  freq.double = freq.double[,-c(2,4)]
  freq.double = cbind( freq.double, round(freq.double[,3]/freq.double[,2] * 100))
  # bin color factor double producer Z1+Z2
  my.calc.fac.double = cut(freq.double[,4],breaks=seq(0,maxfreq,by=maxfreq/10),labels=1:10,include.lowest=TRUE)
  names(my.calc.fac.double) = freq.double[,4]
  
  
  ### combine all frequencies in one table
  my.calc = cbind(my.calc,fac.Z1=as.numeric(my.calc.fac.Z1))
  my.calc = cbind(my.calc,freq.Z2=freq.Z2$x)
  my.calc = cbind(my.calc,fac.Z2=as.numeric(my.calc.fac.Z2))
  my.calc = cbind(my.calc,freq.double=freq.double[,4])
  my.calc = cbind(my.calc,fac.double=as.numeric(my.calc.fac.double))
  my.calc = cbind(my.calc,ncells=length.all$x)
  my.calc = cbind(my.calc,ncells.double=freq.double[,3])
  
  
  ### check the maximum frequency 
  # maxfreq.real = max(this$my.calc[which(!is.na(this$my.calc$freq.double) & this$my.calc$ncells>9),6])
  # printf("MAXIMUM FREQUENCY = %s",maxfreq.real)
  # if ( maxfreq.real > maxfreq ) print("    !!!!! MAXFREQ SHOULD BE HIGHER !!!!!")
  
  ########## DONE CALCULATE BIN PROPERTIES
  
  
  cols = fcs$colors(col)
  cols.heat = cols[,as.numeric(my.calc.fac.double)]
  this$cols.heat = cols.heat
  
  brackets.open = c("(","[")
  ########## PLOT DOUBLE FREQUENCY BINS
  for (x in rownames(tab)) {
    for (y in colnames(tab)) {
      if ( tab[x,y]>=mincells ) {
        brackets.idx.x = brackets.idx.y = 1
        if (x==0) brackets.idx.x = 2
        if (y==0) brackets.idx.y = 2
        
        fact = as.factor(paste(brackets.open[brackets.idx.x],x,',',as.numeric(x)+binsize,'] ',
                               brackets.open[brackets.idx.y],y,',',as.numeric(y)+binsize,']',sep=''))
        
        idx = which(as.character(fact)==as.character(my.calc$Group.1))
        
        if (length(cols.heat[,idx])!=0) {
          if ( !is.na(cols.heat[1,idx]) ) {  
            
            rect(x,y,as.numeric(x)+binsize,as.numeric(y)+binsize,
                 col=eval(parse(text=paste0("rgb(",paste0(cols.heat[,idx],collapse=","),",maxColorValue=255)"))),
                 border=NA
            )
            
            tab3D[x,y] = my.calc[idx,'x']
            tab3D_z2[x,y] = my.calc[idx,'freq.Z2']
          }
        }
      } else {
        tab3D[x,y] = tab3D_z2[x,y] = NA
      }
    }
  }
  ########### DONE PLOT BINS
  
  this$addProdline(cutoffs[c(1,2)])
  
  ########## PLOT LEGEND
  start.legend = par()$usr[2] - 0.28*(par()$usr[2]-par()$usr[1])
  for (step in 1:10) {
    ###### from bottom to top
    ### legend Z1/Z2
    rect(xleft = start.legend+0.023*(par()$usr[2]-par()$usr[1])*step,
         xright = start.legend+0.023*(par()$usr[2]-par()$usr[1])*(step+1),
         ybottom = par()$usr[4]+0.10*(par()$usr[3]+par()$usr[4])-3*binsize,
         ytop = par()$usr[4]+0.10*(par()$usr[3]+par()$usr[4])-2*binsize,
         col=eval(parse(text=paste0("rgb(",paste0(cols[,step],collapse=","),",maxColorValue=255)"))),
         xpd=TRUE,border=NA
    )
    
    if (step==1) {
      text(x=start.legend+0.023*(par()$usr[2]-par()$usr[1])*step,
           y=par()$usr[4]+0.10*(par()$usr[3]+par()$usr[4])-2*binsize,
           label="0",cex=0.8,pos=1,xpd=TRUE)
    } else if (step==10) {
      text(x=start.legend+0.023*(par()$usr[2]-par()$usr[1])*(step+1),
           y=par()$usr[4]+0.10*(par()$usr[3]+par()$usr[4])-2*binsize,
           label=sprintf("%.0f",as.numeric(maxfreq)),cex=0.8,pos=1,xpd=TRUE)
    }
  }
  
  ### x axis label
  text(x = 0.5*(par()$usr[1]+par()$usr[2]),
       y = par()$usr[3] - 1,
       label =feat.X,
       xpd = TRUE
  )
  ### y axis label
  text(x = par()$usr[1] - 0.8,
       y = 0.5*(par()$usr[3]+par()$usr[4]),
       label = feat.Y,
       xpd = TRUE,
       srt=90
  )
  
  ######## from bottom to top
  ### legend title Z1
  text(x = start.legend + 0.8*binsize,
       y = par()$usr[4]+0.075*(par()$usr[3]+par()$usr[4])-3*binsize,
       label=feat.Z1,cex=0.8,pos=2,xpd=TRUE
  )
  ### legend title Z2
  text(x = start.legend + 0.8*binsize,
       y = par()$usr[4]+0.075*(par()$usr[3]+par()$usr[4])-1*binsize,
       label=feat.Z2,cex=0.8,pos=2,xpd=TRUE
  )
  ########## DONE PLOT LEGEND
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
    cols = fcs$colors(col)
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
