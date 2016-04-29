# sequoia
Pedigree inference based on SNP data

Sequoia provides a method to reconstruct multi-generational pedigrees based on SNP data, as described in the manuscript ``Pedigree reconstruction using SNP data: parentage assignment, sibship clustering, and beyond''. The bulk of the algorithm is written in Fortran, to minimise computation times.

This is a beta version: the parentage assignment part works OK, but the sibship clustering may still contain a few bugs, and may not always give the maximum assignment rate and minimum error rate.

For further information, please contact  jisca.huisman@gmail.com
