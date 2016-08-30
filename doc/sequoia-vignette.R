## ---- echo=TRUE, eval=FALSE, results='asis'------------------------------
## # install the package (only required once)
## install.packages("sequoia")
## 
## # set the working directory
## setwd("E:/Sequoia/test")
## 
## # load the package
## library(sequoia)
## 
## # copy the example pedigree and associated life history file to the working
## # directory.
## file.copy(system.file("Ped_HSg5.txt", package="sequoia"), getwd())
## file.copy(system.file("LH_HSg5.txt", package="sequoia"), getwd())
## 
## # simulate genotype data for 200 SNPs, and use otherwise default values
## SimGeno(PedFile = "Ped_HSg5.txt", nSnp = 200)
## 
## # run the preparation step, duplicate checking and parentage assignment,
## # but not yet the slower sibship clustering. Iterate as necessary,
## # weeding out duplicated and erroneous samples (using PLINK's toolkit),
## # and adding estimated birth years to the life history file.
## sequoia(GenoFile = "SimGeno.txt", LifeHistFile = "LH_HSg5.txt", Sibships = FALSE)
## 
## # compare the assigned parents to those in the true pedigree
## PedCompare(PedIN = "Ped_HSg5.txt", PedOUT = "Parents_assigned.txt")
## 
## # run sibship clustering
## sequoia(Prep = FALSE, CheckDup = FALSE, Parentage = FALSE, Sibships = TRUE)
## 
## # compare the assigned real and dummy parents to the true pedigree
## PedCompare(PedIN = "Ped_HSg5.txt", PedOUT = "PedSeq.txt")

## ---- echo=TRUE, eval=FALSE, results='asis'------------------------------
## plink --file mydata --maf 0.4 --indep 50 5 2

## ---- echo=TRUE, eval=FALSE, results='asis'------------------------------
## plink --file mydata --extract plink.prune.in --recodeA --out inputfile_for_sequoia

## ---- echo=TRUE, eval=FALSE, results='asis'------------------------------
## sequoia(RawFile = NULL, GenoFile = NULL, LifeHistFile = NULL,
## 				Prep = TRUE, CheckDup = TRUE, Parentage = TRUE, Sibships = TRUE)

