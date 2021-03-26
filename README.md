# sequoia
Sequoia provides a method to reconstruct multi-generational pedigrees based on SNP data, as described in the manuscript `Pedigree reconstruction using SNP data: parentage assignment, sibship clustering, and beyond` ( http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12665/full ). 

The latest stable and thoroughly bug and performance checked version of the package is available from CRAN, and can be installed using 
`install.packages("sequoia")`

The version here may be a beta version with new features, or contain bug fixes, or be the current CRAN version. Previous versions can be found in /sequoia_archives .

The bulk of the algorithm is written in Fortran, to minimise computation times, and needs compilation. A pre-compiled .zip binary file for Windows and the current R version can be found in /sequoia_archives , as well as .tar.gz source archives. You can install these using
`install.packages("C:/file/to/path/sequoia_1.1.1.zip",  repos = NULL)` (sometimes you need to turn R off & on before it works). 

Detailed instructions on how to use the package are available via the R command
`vignette("sequoia")` , after loading the package with `library(sequoia)`. 

## Not R 
A stand-alone Fortran version can be found in /sequoia_notR . It does not offer all the features of the R package, but does not require the genetic data to pass through R. Its version will typically be functionally identical to the latest R package beta version. 


For further information, questions or comments, feel free to contact me at jisca.huisman @ gmail.com
