if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")

# Set path to working directory.This is were the repository will be stored
# sets Users home directory to working directory 
setwd(dir = path.expand("~")) 

# download a .zip file of the repository
URL="https://github.com/InesHo/PRI-demonstration/archive/refs/heads/main.zip"
download.file(URL, destfile = "PRI-demo.zip")

# unzip file
unzip(zipfile = "~/PRI-demo.zip")
