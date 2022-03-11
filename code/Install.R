if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!requireNamespace("flowCore", quietly = TRUE))
    BiocManager::install("flowCore")

# Set path to working directory.This is were the repository will be stored
# sets Users home directory to working directory 
setwd(dir = path.expand("~")) 

# download a .zip file of the repository
URL="https://github.com/InesHo/PRI-demonstration/archive/refs/heads/main.zip"
download.file(URL, destfile = "PRI-demo.zip")

# unzip file
unzip(zipfile = "~/PRI-demo.zip")

setwd("~/PRI-demonstration-main")

cat("\nThe PRI repository folder 'PRI-demonstration-main' is stored in your home folder and now set as the current working directory.\n\n")
cat("You can now follow the tutorial under\n\n http://htmlpreview.github.io/?https://github.com/InesHo/Pri-demonstration/blob/master/code/demonstration_code.html \n\nto use PRI.")
