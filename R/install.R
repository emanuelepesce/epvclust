# ------------------------------------------------------------------------
# Install file of epvclust

# Clean workspace
rm(list = ls())

# Compiling
setwd(dir = "./epvclust")
library(devtools)
library(roxygen2)
document()


# Installation
setwd(dir = "..")
install("epvclust")


# Run examples
library(epvclust)

ptm <- proc.time()

res <- epvdriver(filename = "./../exampleData/medium.csv", nboot = 10)
plot(res);

proc.time() - ptm


