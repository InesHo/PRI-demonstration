# Demonstration for Data analysis with PRI - Pattern Recognition with Immune cells

# author: "Yen Hoang & Ines Hoppe"
# date: "11/30/2021"

##################################################n
#################################

# Load flowCore library and source PRI functions 
rm(list = ls())
library(flowCore)
source("../code/binplot.R")

### 1 process data and create plots for general PRI example
## 1.1 Data preparation

# READ IN UNTREATED SAMPLE 
path.folder <- "../data/"
file.untr <- "CD4_TIN_BLD1_Untreated_Day3.fcs"
data.untr <- read.FCS(file=file.path(path.folder,file.untr),which.lines=1)
keys.untr <- keyword(data.untr)
aliases <- unlist(keyword(data.untr, c("$P5S","$P10S","$P28S","$P29S","$P33S")))
channels <- unlist(keyword(data.untr,c("$P5N","$P10N","$P28N","$P29N","$P33N")))
map.untr <- data.frame(alias = aliases, channels = channels)
# load data
exprs.untr <- as.data.frame(exprs(read.FCS(file=file.path(path.folder,file.untr),
                                           channel_alias = map.untr, transformation=NULL)))
# rename column names for channels we are interested in
for ( i in 1:ncol(exprs.untr)) {
  colnames(exprs.untr)[i] = strsplit(colnames(exprs.untr)[i]," ")[[1]][1]
}
# transform data with asinh()
exprs.untr <- asinh(exprs.untr)

## 1.2 Data visualization
# INITIATING PARAMETERS
featX = "CD90"
featY = "CD44"
featZ1 = "CD27"
cutoffs.man = c(6,4,2)
binSize = 0.2

# PLOTTING PARAMETERS 
png("../images/general_example.png")
par(mfrow = c(1,1),
    cex.lab = 1.4, cex.axis = 1.4,
    pin = c(4,4),
    mgp = c(2.0, 0, 0))

# CALL PLOT FUNCTIONS
binplot_table(
  data = exprs.untr, 
  feat.X = featX, 
  feat.Y = featY, 
  feat.Z1 = featZ1, 
  calc = "MFI+",
  binsize = binSize,
  cutoffs = cutoffs.man,
  plot.range=c(0.05,12,0.05,12)
)
dev.off()

### 2 Use of PRI to Identify Differing Activation Patterns 
## 2.1  Data Preparation

# READ IN TREATED SAMPLE 
file.b6ab <- "CD4_TIN_BLD2_B6Antibodies_Day3.fcs"
data.b6ab <- read.FCS(file=file.path(path.folder,file.b6ab),which.lines=1)
keys.b6ab <- keyword(data.b6ab)
aliases.b6ab <- unlist(keyword(data.b6ab,c("$P5S","$P10S","$P28S","$P29S","$P33S")))
channels.b6ab <- unlist(keyword(data.b6ab,c("$P5N","$P10N","$P28N","$P29N","$P33N")))
map.b6ab <- data.frame(alias = aliases.b6ab, channels = channels.b6ab)

# load data
exprs.b6ab <- as.data.frame(exprs(read.FCS(file=file.path(path.folder,file.b6ab),
                                           channel_alias = map.b6ab, transformation=NULL)))
# rename column names for channels we are interested in 
for ( i in 1:ncol(exprs.b6ab)) {
  colnames(exprs.b6ab)[i] = strsplit(colnames(exprs.b6ab)[i]," ")[[1]][1]
}
# transform data with asinh()
exprs.b6ab <- asinh(exprs.b6ab)

## 2.2   Binplots Visualizing MFI+ of CD27 for T-cells From Blood Samples

# INITIATING PARAMETERS
featX = "CD90"
featY = "CD44"
featZ1 = "CD27"
cutoffs.man = c(6,4,2)
binSize = 0.2


# PLOTTING PARAMETERS FOR TWO PLOTS
png("../images/triplot_CD27.png")
par(mfrow = c(1,2),
    cex.lab = 1.4, cex.axis = 1.4,
    mgp = c(2.0, 0, 0),
    mar = c(2.4,2,2.5,0),
    oma = c(1,1,1,1))

# CALL PLOT FUNCTIONS
binplot_table(
  data = exprs.untr, 
  feat.X = featX, 
  feat.Y = featY, 
  feat.Z1 = featZ1, 
  calc = "MFI+",
  binsize = binSize,
  cutoffs = cutoffs.man,
  plot.range=c(0.05,12,0.05,12)
)

binplot_table(
  data = exprs.b6ab, 
  feat.X = featX, 
  feat.Y = featY, 
  feat.Z1 = featZ1, 
  calc = "MFI+",
  binsize = binSize,
  cutoffs = cutoffs.man,
  plot.range=c(0.05,12,0.05,12)
)
dev.off()

## 2.3   Binplots For Cell Characterisation in Treated Blood Sample using T.bet and Foxp3
# PLOTTING PARAMETERS FOR TWO PLOTS
png("../images/triplot_Tbet_Foxp3.png")
par(mfrow = c(1,2),
    cex.lab = 1.4, cex.axis = 1.4,
    mgp = c(2.0, 0, 0),
    mar = c(2.4,2,2.5,0),
    oma = c(1,1,1,1))

# CALL PLOT FUNCTIONS
binplot_table(
  data = exprs.b6ab, 
  feat.X = featX, 
  feat.Y = featY, 
  feat.Z1 = "T.bet", 
  calc = "MFI+",
  binsize = binSize,
  cutoffs = c(6,4,0.4),
  plot.range=c(0.05,12,0.05,12)
)

binplot_table(
  data = exprs.b6ab, 
  feat.X = featX, 
  feat.Y = featY, 
  feat.Z1 = "Foxp3", 
  calc = "MFI+",
  binsize = binSize,
  cutoffs = c(6,4,0),
  plot.range=c(0.05,12,0.05,12)
)
dev.off()
