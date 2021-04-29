# sequoia

## General
Sequoia provides a method to reconstruct multi-generational pedigrees based on SNP data, as described in the manuscript `Pedigree reconstruction using SNP data: parentage assignment, sibship clustering, and beyond` ( http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12665/full ). It performs parentage assignment (parent genotyped), sibship clustering (parent not genotyped), and grandparent assignment to link sibship clusters to the rest of the pedigree. 

The R package also includes various functions to compare pedigrees and to check if a pedigree is accordant with the SNP data. 

For successfull pedigree reconstruction, `sequoia` ideally requires data on a few hundred SNPs with high MAF, low missingness, and low genotyping error rate; and sex and birth year information of as many individuals as possible. It usually works adequately when these conditions are not met. 


## Installation
The latest stable and thoroughly bug and performance checked version of the package is available from CRAN, and can be downloaded installed using (in R)
```
install.packages("sequoia")
```

The latest not quite so thoroughly tested version can be downloaded and installed using
```
remotes::install_github("JiscaH/sequoia", ref="stable")
```

This package requires compilation, as the bulk of the algorithm is written in Fortran. A pre-compiled .zip binary file (for Windows and the current R version only) can be found in [sequoia_archives](https://github.com/JiscaH/sequoia_archives) , as well as .tar.gz source archives. You can download these to your hard drive, and then install using
```
install.packages("C:/file/to/path/sequoia_1.1.1.zip",  repos = NULL)`
```
(sometimes you need to turn R off & on before it works). 


## Running 
The function to perform pedigree reconstruction is `sequoia()`:
```
library(sequoia)   # load the package
data(SimGeno_example, LH_HSg5, package="sequoia")  # example data
SeqOUT <- sequoia(GenoM = SimGeno_example, LifeHistData = LH_HSg5, 
                  Err = 0.005, Module="ped", quiet="verbose", Plot=TRUE)
head(SeqOUT$Pedigree)   #  the reconstructed pedigree
```

Detailed instructions on how to use the package are available in the vignette (vignettes/vignette-main.pdf), available in R
via `vignette("sequoia")` or `??sequoia` (that's 2 question marks).


## Versions
The version here may be a beta version with new features, or contain bug fixes, or be the current CRAN version. Previous versions can be found in [sequoia_archives](https://github.com/JiscaH/sequoia_archives), including pre-compiled Windows .zip files (which are not archived by CRAN).  
 

## Not R 
A stand-alone Fortran version can be found in [sequoia_notR](https://github.com/JiscaH/sequoia_notR) . It does not require the genetic data to pass through R, which may be nearly impossible for very large datasets. Note that it does not offer all the features of the R package, as detailed in its manual. Its version will typically be functionally identical to the latest R package beta version. 
