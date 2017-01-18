# sequoia
Sequoia provides a method to reconstruct multi-generational pedigrees based on SNP data, as described in the manuscript `Pedigree reconstruction using SNP data: parentage assignment, sibship clustering, and beyond`. The bulk of the algorithm is written in Fortran, to minimise computation times.

The package has been submitted to CRAN and should hopefully be available there soon.

Provided you have a Fortran compiler on your computer, you can install the package using  
`library(devtools)`    
`install_github("JiscaH/sequoia")`

Alternatively, a binary for the package for 64-bit windows with R 3.3.1 is available in https://github.com/JiscaH/Sequoia-source-code, 
as well as the Fortran source code. 

For further information, questions or comments, please contact me at jisca.huisman@gmail.com
