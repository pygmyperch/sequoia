# sequoia
Pedigree inference based on SNP data

Sequoia provides a method to reconstruct multi-generational pedigrees based on SNP data, as described in the manuscript ``Pedigree reconstruction using SNP data: parentage assignment, sibship clustering, and beyond''. The bulk of the algorithm is written in Fortran, to minimise computation times.

If devtools::install_github("jiscah/sequoia") does not work, the following might: download the zip, install the package from the zip file using 
install.packages("E:/Sequoia/test/sequoia_0.3.zip", repos=NULL, type="binary")
, and if necessary rename the folder in R\library from `sequoia-master' to `sequoia'.

For further information, please contact  jisca.huisman@gmail.com


# A quick useage example with simulated data  
For details, see the R vignette and ?sequoia

setwd("E:/Sequoia/test")  
library(sequoia)  

copy the example pedigree and associated life history file to the working directory:  
file.copy(system.file("Ped_HSg5.txt", package="sequoia"), getwd())  
file.copy(system.file("LH_HSg5.txt", package="sequoia"), getwd())  

simulate genotype data for 200 SNPs, and use otherwise default values:  
SimGeno(PedFile = "Ped_HSg5.txt", nSnp = 200)

Run Sequoia:  
sequoia(GenoFile = "SimGeno.txt", LifeHistFile = "LH_HSg5.txt")

Compare the assigned parents to those in the true pedigree:  
PedCompare(PedIN = "Ped_HSg5.txt", PedOUT = "PedSeq.txt")
