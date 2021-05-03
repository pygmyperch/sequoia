# sequoia

## General
Sequoia reconstructs multi-generational pedigrees from SNP data, as described in `Pedigree reconstruction using SNP data: parentage assignment, sibship clustering, and beyond` ( http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12665/full ). It performs parentage assignment (parent genotyped) and sibship clustering (parent not genotyped), as well as grandparent assignment to link sibship clusters to the rest of the pedigree. 

The R package also includes a range of other functions, amongst others to check if a pedigree is accordant with the SNP data, and to compare two pedigrees.

For successfull pedigree reconstruction, `sequoia` ideally requires data on a few hundred SNPs with high MAF, low missingness, low genotyping error rate, and in low LD, but real-world data usually works well too. In addition it is useful to have sex and birth year information of as many individuals as possible.  


## Download & installation
The latest stable and thoroughly bug and performance checked version of the package is available from CRAN, and can be downloaded installed using (in R)
```
install.packages("sequoia")
```

The latest not quite so thoroughly tested version can be downloaded and installed using
```
remotes::install_github("JiscaH/sequoia", ref="stable")
```

Note that this package requires compilation, as the bulk of the algorithm is written in Fortran. A pre-compiled .zip binary file of the development version (for Windows and the current R version only) can be found in [sequoia_archives](https://github.com/JiscaH/sequoia_archives). You can download this to your hard drive, and then install using
```
install.packages("C:/file/to/path/sequoia_2.3.3.zip",  repos = NULL)`
```
Sometimes it seems necessary to first rename the file to 'sequoia.zip', if R can't find the package back after installation, and sometimes you need to turn R off & on before it works.  

The [archives](https://github.com/JiscaH/sequoia_archives) also has the source files and pre-compiled Windows and MacOS binaries of old versions that used to be on CRAN, as well as source files of old versions that never made it to CRAN. 


## Running 
The function to perform pedigree reconstruction is also called `sequoia()`:
```
# load the package
library(sequoia)  

# load example data 
data(SimGeno_example, LH_HSg5, package="sequoia")  

# run pedigree reconstruction
SeqOUT <- sequoia(GenoM = SimGeno_example, 
                  LifeHistData = LH_HSg5, 
                  Err = 0.005,   # genotyping error rate
                  Module="ped", 
                  quiet="verbose", 
                  Plot=TRUE)
# the result is a list with the pedigree, run parameters, 
# and various other elements.                 

# graphical summary of results
SummarySeq(SeqOUT)
```

Detailed instructions on how to use the package are available in the vignette (vignettes/vignette-main.pdf), available in R
via `vignette("sequoia")` or `help(package="sequoia")`.



## Not R 
A stand-alone Fortran version can be found in [sequoia_notR](https://github.com/JiscaH/sequoia_notR) . It does not require the genetic data to pass through R, which may be nearly impossible for very large datasets. Note that it does not offer all the features of the R package, as detailed in its manual. Its version will typically be functionally identical to the latest R package development version. 
