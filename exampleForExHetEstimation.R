#Example use case for using est.rf.exHet() as part of the Truong and McCormick et al. (2014)
#   manuscript. This example requires the R package devtools to download and compile the hetexp
#   branch of the modified R/qtl code base (originally forked from Prof. Karl Broman's R/qtl repository).
#
#Written by Ryan McCormick and Sandra Truong at Texas A&M University, July 2014.
#
#Much of the code originates or is inspired by Karl Broman's 
#   document from http://www.rqtl.org/tutorials/geneticmaps.pdf
#
#Provided as is, without warranty, without guarantee.
# *     This program is free software; you can redistribute it and/or
# *     modify it under the terms of the GNU General Public License,
# *     version 3, as published by the Free Software Foundation.
# *
# *     This program is distributed in the hope that it will be useful,
# *     but without any warranty; without even the implied warranty of
# *     merchantability or fitness for a particular purpose.  See the GNU
# *     General Public License, version 3, for more details.
# *
# *     A copy of the GNU General Public License, version 3, is available
# *     at http://www.r-project.org/Licenses/GPL-3

library(devtools); install_github("qtl", "MulletLab", "hetexp")   #The Rqtl branch with heterzygosity functionality
library(qtl)

#Generation interval of the selfed cross:
generation_interval=7

#Change the file path to reflect system if necessary.
input_file_directory <- file.path("./")
setwd(input_file_directory)

#This function writes a BCsFt cross to an output file; it is a workaround for an in RQTL (may have been fixed since this was written).
write.cross.BCsFt <- function(cross_cross, str_filestem="BCsFTout", str_format="csvr")
{
  class(cross_cross)[1] <- "f2"
  write.cross(cross_cross, str_format, filestem=str_filestem)
}

###################################
###################################
###################################
#Beginning with Simulated (n=1000, t=7, h=0.6373) – No errors, no missing data
#
#
#

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – No errors, no missing data
##### Analyzing with est.rf.exHet(h=0.6373)

input_file_name <- file.path("./simulated_genotypes_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.6373027617 #e^(ln(big_h)/(generation_interval-1))

cross_filteredcross <- est.rf.exHet(cross_filteredcross, het=little_h)
newmap <- mapFromRF(cross_filteredcross, mapfunc="haldane")
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – No errors, no missing data
##### Analyzing with est.rf.exHet(h=0.5)

input_file_name <- file.path("./simulated_genotypes_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

cross_filteredcross <- est.rf.exHet(cross_filteredcross, het=little_h)
newmap <- mapFromRF(cross_filteredcross, mapfunc="haldane")
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – No errors, no missing data
##### Analyzing with est.rf() (Which uses a Mendelian model where h=0.5)

input_file_name <- file.path("./simulated_genotypes_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

cross_filteredcross <- est.rf(cross_filteredcross)
newmap <- mapFromRF(cross_filteredcross, mapfunc="haldane")
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – No errors, no missing data
##### Analyzing with est.map(error.prob=0.01) (Which uses a Mendelian model where h=0.5, and performs
#####   multipoint estimation using a hidden Markov model).

input_file_name <- file.path("./simulated_genotypes_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

newmap <- est.map(cross_filteredcross, error.prob=0.01, map.function="haldane") #This uses the hidden Markov model
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – No errors, no missing data
##### Analyzing as a fixed RIL with est.map(error.prob=0.01) (Which uses a Mendelian model where h=0.5, and performs
#####   multipoint estimation using a hidden Markov model).

input_file_name <- file.path("./simulated_genotypes_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))
cross_filteredcross <- convert2riself(cross_filteredcross)
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

newmap <- est.map(cross_filteredcross, error.prob=0.01, map.function="haldane") #This uses the hidden Markov model
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

###################################
###################################
###################################
#Beginning with Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%)
#
#
#

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%)
##### Analyzing with est.rf.exHet(h=0.6373)

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.6373027617 #e^(ln(big_h)/(generation_interval-1))

cross_filteredcross <- est.rf.exHet(cross_filteredcross, het=little_h)
newmap <- mapFromRF(cross_filteredcross, mapfunc="haldane")
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)
write.cross.BCsFt(cross_filteredcross, str_filestem="simulated_genotypes_errors_and_missing_rqtl_mapEstimated", str_format="csvr")

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%)
##### Analyzing with est.rf.exHet(h=0.5)

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

cross_filteredcross <- est.rf.exHet(cross_filteredcross, het=little_h)
newmap <- mapFromRF(cross_filteredcross, mapfunc="haldane")
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%)
##### Analyzing with est.rf() (Which uses a Mendelian model where h=0.5)

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

cross_filteredcross <- est.rf(cross_filteredcross)
newmap <- mapFromRF(cross_filteredcross, mapfunc="haldane")
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%)
##### Analyzing with est.map(error.prob=0.01) (Which uses a Mendelian model where h=0.5, and performs
#####   multipoint estimation using a hidden Markov model).

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

newmap <- est.map(cross_filteredcross, error.prob=0.01, map.function="haldane") #This uses the hidden Markov model
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%)
##### Analyzing as fixed RIL with est.map(error.prob=0.01) (Which uses a Mendelian model where h=0.5, and performs
#####   multipoint estimation using a hidden Markov model).

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl.csv")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("aa","ab","bb","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))
cross_filteredcross <- convert2riself(cross_filteredcross)
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

newmap <- est.map(cross_filteredcross, error.prob=0.01, map.function="haldane") #This uses the hidden Markov model
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

###################################
###################################
###################################
#Beginning Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%) - No SDCO
#
#
#

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%) - No SDCO
##### Analyzing with est.rf.exHet(h=0.6373)

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl_mapEstimated.csv-SDCOsToMissingV02")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("AA","AB","BB","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.6373027617 #e^(ln(big_h)/(generation_interval-1))

cross_filteredcross <- est.rf.exHet(cross_filteredcross, het=little_h)
newmap <- mapFromRF(cross_filteredcross, mapfunc="haldane")
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%) - No SDCO
##### Analyzing with est.rf.exHet(h=0.5)

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl_mapEstimated.csv-SDCOsToMissingV02")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("AA","AB","BB","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

cross_filteredcross <- est.rf.exHet(cross_filteredcross, het=little_h)
newmap <- mapFromRF(cross_filteredcross, mapfunc="haldane")
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%) - No SDCO
##### Analyzing with est.rf() (Which uses a Mendelian model where h=0.5)

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl_mapEstimated.csv-SDCOsToMissingV02")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("AA","AB","BB","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

cross_filteredcross <- est.rf(cross_filteredcross)
newmap <- mapFromRF(cross_filteredcross, mapfunc="haldane")
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%) - No SDCO
##### Analyzing with est.map(error.prob=0.01) (Which uses a Mendelian model where h=0.5, and performs
#####   multipoint estimation using a hidden Markov model).

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl_mapEstimated.csv-SDCOsToMissingV02")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("AA","AB","BB","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

newmap <- est.map(cross_filteredcross, error.prob=0.01, map.function="haldane") #This uses the hidden Markov model
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)

#################
#################
##### Simulated (n=1000, t=7, h=0.6373) – Errors (1%) and missing data (5%) - No SDCO
##### Analyzing as fixed RIL with est.map(error.prob=0.01) (Which uses a Mendelian model where h=0.5, and performs
#####   multipoint estimation using a hidden Markov model).

input_file_name <- file.path("./simulated_genotypes_errors_and_missing_rqtl_mapEstimated.csv-SDCOsToMissingV02")
cross_inputcross <- read.cross("csvr", input_file_directory, input_file_name, BC.gen=0, F.gen=generation_interval, genotypes=c("AA","AB","BB","D","C"))
cross_filteredcross <- cross_inputcross
print(summary(cross_filteredcross))
cross_filteredcross <- convert2riself(cross_filteredcross)
print(summary(cross_filteredcross))

#Calculated heterozygosity maintained per generation
little_h = 0.5 #e^(ln(big_h)/(generation_interval-1))

newmap <- est.map(cross_filteredcross, error.prob=0.01, map.function="haldane") #This uses the hidden Markov model
cross_filteredcross <- replace.map(cross_filteredcross, newmap)

map <- pull.map(cross_filteredcross)
map_size <- 0
for (j in 1:nchr(cross_filteredcross)) {
  last_marker_index <- length(map[[j]])
  print(c(last_marker_index, map[[j]][last_marker_index]))
  map_size <- map_size + map[[j]][last_marker_index]
}
print(map_size)