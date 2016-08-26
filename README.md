# sequoia
Pedigree inference based on SNP data

Sequoia provides a method to reconstruct multi-generational pedigrees based on SNP data, as described in the manuscript ``Pedigree reconstruction using SNP data: parentage assignment, sibship clustering, and beyond''. The bulk of the algorithm is written in Fortran, to minimise computation times.

Useage is explained in the R vignette, or see ?sequoia 

If devtools::install_github("jiscah/sequoia") does not work, the following might: download the zip, install the package from the zip file using 
install.packages("E:/Sequoia/test/sequoia_0.3.zip", repos=NULL, type="binary")
, and if necessary rename the folder in R\library from `sequoia-master' to `sequoia'.

For further information, please contact  jisca.huisman@gmail.com
