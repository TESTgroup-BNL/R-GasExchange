#!/usr/bin/env Rscript
#
# Install required packages for R-GasExchange
#
list.of.packages <- c('XML','DEoptim') 

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  print("installing : ")
  print(new.packages)
  install.packages(new.packages, repos="http://cran.us.r-project.org")
}

# Look for suggested packages
if(!('testthat') %in% installed.packages()[,"Package"]) {
  print("Suggested install testthat not installed")
  print("You can install the driver using the following command in R.")
  print("install.packages('testthat', repos='http://cran.rstudio.com/')")
} else if (!('devtools') %in% installed.packages()[,"Package"]) {
  print("Suggested install devtools not installed")
  print("You can install the driver using the following command in R.")
  print("install.packages('devtools', repos='http://cran.rstudio.com/')")
}