# sequoia
Sequoia provides a method to reconstruct multi-generational pedigrees based on SNP data, as described in the manuscript `Pedigree reconstruction using SNP data: parentage assignment, sibship clustering, and beyond`. The bulk of the algorithm is written in Fortran, to minimise computation times.

The package has been submitted to CRAN and should hopefully be available there soon.

If devtools::install_github("jiscah/sequoia") does not work, the following might: download the zip, install the package from the zip file using 
install.packages("E:/Sequoia/test/sequoia_0.3.zip", repos=NULL, type="binary")
, and if necessary rename the folder in R\library from `sequoia-master` to `sequoia`.

An working example is provided under ?sequoia, with details in the R vignette. Fortran source code and a binary for the package for 64-bit windows with R 3.3.1 are available in https://github.com/JiscaH/Sequoia-source-code

For further information, questions or comments, please contact me at jisca.huisman@gmail.com
