#Lets import the necessary packages


# List of packages to check
required_packages <- c(
  "ppcor", "mgcv", "gplm", "FOCI", "dplyr", "tidyr", 
  "ggplot2", "XICOR", "forecast", "minerva", "tsDyn"
)

# Identify which packages are missing
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# Install missing packages if there are any
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
  message("Installed packages: ", paste(missing_packages, collapse = ", "))
} else {
  message("All packages are already installed!")
}

library(ppcor)
library(mgcv)
library(gplm)
library(FOCI) #Codec
library(dplyr)
library(tidyr)
library(ggplot2)
library(XICOR)
library(forecast)
library(minerva)
library(tsDyn) 

### Creating a function for partial correlation-------------------------------


partcor<-function(datos, method="kendall"){
  if (method=="xi"){a<-matrix(0,ncol=ncol(datos),nrow=ncol(datos))
  for (i in 1:ncol(datos)){
    for (j in 1:ncol(datos)){
      a[i,j]<-XICOR::xicor(datos[i],datos[j])
    }
  }}else{a<-cor(datos,method=method)}
  
  b<-solve(a)
  c<-matrix(0,ncol = ncol(a), nrow = ncol(a))
  for (i in 1:ncol(a)){
    for (j in 1:ncol(a)){
      c[i,j]=-b[i,j]/(sqrt(b[i,i])*sqrt(b[j,j]))
      if (i==j){c[i,j]=-c[i,j]}
    }
  }
  return(c)}

# Function to show the progress bar of the simulations
update_progress <- function(current, total) {
  message <- paste("|--", current, "/", total, "--|")
  cat("\r", message)
  flush.console()
}

