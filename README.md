# sequoia
Sequoia provides a method to reconstruct multi-generational pedigrees based on SNP data, as described in the manuscript `Pedigree reconstruction using SNP data: parentage assignment, sibship clustering, and beyond` ( http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12665/full ). The bulk of the algorithm is written in Fortran, to minimise computation times.

The package is available from CRAN, and can be installed using
`install.packages("sequoia")`

The version here may sometimes be newer, as it takes about a week for updates to turn into compiled packages on all CRAN servers. If you have a Fortran compiler on your computer, download the .tar.gz source, else if you have a windows machine and the current version of R, the .zip binary might work. You can install these using
`install.packages("C:/file/to/path/sequoia_0.9.3.zip",  repos = NULL)`
followed by turning R off & on. 

You can access detailed instructions on how to use the package with the command
`vignette("sequoia")`

Note that the Fortran stand-alone version is not as up to date as this R version!

For further information, questions or comments, please contact me at jisca.huisman@gmail.com
