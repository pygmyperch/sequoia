#============================================================================
#============================================================================
# preparations

#' Convert genotype file if needed, and save parameter values.
#'
#' Write the user-specified parameter values to the (temporary) file
#' `SequoiaSpecs.txt'. If GenoFormat="raw" or "col", call \code{\link{GenoConvert}}.
#'
#'  Please do not add or delete rows in the SequoiaSpecs text file. Changing
#'  values is possible, but increasing the number of SNPs or individuals beyond
#'  the numbers present in the files, or refering to non-existing files, will
#'  cause errors and may cause R to crash.
#'
#' @param GenoFile character string with name of genotype file
#' @param GenoFormat One of "raw", "ped", "col" or "seq", see ?GenoConvert.
#' @param LifeHistFile Life history data, with 3 columns:
#'  ID: max. 30 characters long,
#'  Sex: 1 = females, 2 = males, other numbers = unkown,
#'  BirthYear: Use negative numbers to denote missing values. If the species
#'  has multiple generations per year, use an integer coding such that the
#'  candidate parents' `Birth year' is at least one larger than their putative
#'  offspring.
#' @param Err Estimated genotyping error rate.
#' @param MaxMismatch Maximum number of loci at which candidate parent and
#'  offspring are allowed to be opposite homozygotes, or be excluded.
#' @param Tfilter Threshold log-likelihood ratio between a proposed
#'   relationship versus unrelated, to select candidate relatives. Typically a
#'   negative value, related to the fact that unconditional likelihoods are
#'   calculated during the filtering steps. More negative values may decrease
#'   non-assignment, but will increase computational time.
#' @param Tassign Minimum log-likelihood ratio required for acceptance of
#'  proposed relationship, relative to next most likely relationship. Higher
#'  values result in more conservative assignments.
#' @param MaxSibshipSize  Maximum number of offspring for a single individual
#'  (a generous safety margin is advised).
#' @param NumSibRounds Number of iterations of sibship clustering.
#' @param Complexity  Either "full" (default), "simp" (no explicit consideration
#'   of inbred relationships) or "mono" (monogamous breeding system)
#' @param ParentageFile Filename for assigned parents
#' @param AgePriorFile Filename for ageprior (in & out, if applicable)
#' @param PedigreeFile Filename for pedigree after sibship clustering
#'
#' @return A dataframe with the number of individuals found in the genotype file and in the lifehistory file, and the number of SNPs scored.

SeqPrep <- function(GenoFile = NULL,
                     GenoFormat = "seq",
                     LifeHistFile = NULL,
                     Err = 0.0001,
                     MaxMismatch = 3,
                     Tfilter = -2.0,
                     Tassign = 0.5,
            				 MaxSibshipSize = 100,
            				 NumSibRounds = 3,
            				 Complexity = "full",
            				 AgePriorFile = "AgePriors.txt",
            				 ParentageFile = "ParentsAssigned.txt",
            				 PedigreeFile = "PedSeq.txt")

{
  if (is.null(GenoFile)) stop("Please provide genotype file")
  if (is.null(LifeHistFile)) stop("Please provide lifehistory file")
  if (GenoFormat != "seq") {
    GenoConvert(InFile = GenoFile, InFormat = GenoFormat)
    GenoFile = "GenoForSequoia.txt"
  } else {  # check if genotype file is in valid format
    tmp <- strsplit(readLines(GenoFile, n=1), split=" ")
    if (!all(tmp[[1]][-1] %in% c(0,1,2,-9))) stop("GenoFile not in <seq> format")
  }

  # count number of individuals & number of SNPs in genotype file
  IDs_geno <- scan(GenoFile, what = "character", flush = TRUE, quiet = TRUE)
  nIndG <- length(IDs_geno)
  nSnp <- length(scan(GenoFile, what = "character", nlines = 1, quiet = TRUE)) - 1
  if (any(sapply(IDs_geno, nchar)>=29)) {
    warning("IDs of one or more individuals too long (max 29 characters)")
  }
  LifeHist_tmp <- utils::read.table(LifeHistFile, header=T)
  names(LifeHist_tmp) <- c("ID", "Sex", "BY")
  nIndLH <- nrow(LifeHist_tmp)
  if (nIndLH >1) {
    if (all(LifeHist_tmp$BY < 0)) {
      nAgeClasses <- 1
    } else {
      nAgeClasses <- with(LifeHist_tmp, diff(range(BY[BY >= 0 & ID %in% IDs_geno],
                                                 na.rm = TRUE))) + 1
    }
  } else {
    nAgeClasses <- 1
  }
  if (nAgeClasses > 100) {
    nAgeClasses <- length(table(LifeHist_tmp$BY[LifeHist_tmp$BY >= 0]))
  }
  rm(LifeHist_tmp)
  Complex <- switch(Complexity, full = 2, simp = 1, mono = 0)

  # write the specifications to a file, to be read in by Fortran
  Specs <- c(GenotypeFilename = GenoFile,
             LifehistoryFilename = LifeHistFile,
             NumberIndivGenotyped = nIndG,
             NumberIndivLifehist = nIndLH,
             NumberSnps = nSnp,
             GenotypingErrorRate = Err,
             MaxMismatch = MaxMismatch,
             Tfilter = Tfilter,
             Tassign = Tassign,
             ParentageFileName = ParentageFile,
      			 AgePriorFileName = AgePriorFile,
      			 nAgeClasses = nAgeClasses,   # max age difference
      			 MaxSibshipSize = MaxSibshipSize,
      			 NumSibRounds = NumSibRounds,
      			 PedigreeFileName = PedigreeFile,
      			 DummyPrefixFemale = "F",
      			 DummyPrefixMale = "M",
      			 Complexity = Complex)
  OPT <- options()
  options(scipen = 10)
  utils::write.table(as.data.frame(Specs), file="SequoiaSpecs.txt",
              sep = "\t,\t", quote = FALSE, col.names = FALSE)
  options(OPT)
  data.frame(n=c(nIndG = nIndG, nIndLH = nIndLH, nSnp = nSnp))
}



#============================================================================
#============================================================================

#' Check data for duplicates.
#'
#' Check the genotype and life history files for duplicate IDs (not permitted)
#' and duplicated genotypes (not advised), and count how many individuals in
#' the genotype file are included in the life history file (permitted).
#'
#' The filenames are taken from the SequioiaSpecs text file, generated by
#' \code{\link{SeqPrep}}. The order of IDs in the genotype and life history file
#' is not required to be identical, and not all individuals in the genotype file
#' are required to be in the life history file.
#'
#' @return A dataframe with the numbers of duplicated genotypes, duplicated IDs in the genotype file, number of duplicated IDs in the life history file, and number of individuals with unknown sex.
#'
#' @useDynLib sequoia

SeqDup <- function()
{
  DUP <- .Fortran("duplicates",
                  nDupGenoID = integer(1),
                  nDupLhID = integer(1),
                  nDupGenos = integer(1),
                  nSexless = integer(1),
                  PACKAGE = "sequoia")

  duplicates <- c("duplicated genotypes" = DUP$nDupGenos,
    "duplicated IDs (genotype file)" = DUP$nDupGenoID,
    "duplicated IDs (lifehist. file)" = DUP$nDupLhID,
    "no sex or not in lifehist. file" = DUP$nSexless)

  if (DUP[[1]]>0) message("Duplicate genotypes found, please remove before continuing.")
  if (DUP[[2]]>0) message("Duplicate IDs found in genotype file, please remove before continuing.")
  if (DUP[[3]]>0) message("Duplicate IDs found in lifehistory file, first one will be used.")
  data.frame(n = duplicates)
}


#=====================================================================
#=====================================================================

#' Parentage assignment
#'
#' Assign genotyped parents to genotyped individuals. Parameter values and
#' filenames are taken from the SequioiaSpecs text file, generated by
#' \code{\link{SeqPrep}}.
#'
#' The output file, by default `ParentsAssigned.txt', contains the following
#' columns:
#' \describe{
#'   \item{id}{Individual ID}
#'   \item{dam}{Assigned mother, or NA}
#'   \item{sire}{Assigned father, or NA}
#'   \item{LLR_dam}{Log-Likelihood Ratio (LLR) of this female being the mother,
#'      versus the next most likely relationship between the focal individual and
#'      this female (see ... for relationships considered)}
#'   \item{LLR_sire}{idem, for male parent}
#'   \item{LLR_pair}{trio LLR, versus the next most likely configuration
#'      (containing one or neither parent)}
#'   \item{OH_dam}{Number of loci at which the offspring and mother are
#'      opposite homozygotes}
#'   \item{OH_sire}{idem, for male parent}
#'   \item{RowO}{Offspring's row identifier in the genotype file. Used by
#'      \code{link{SeqSibs}}, please do not modify}
#'   \item{RowD}{Dam's row identifier}
#'   \item{RowS}{Sire's row identifier}
#' }
#'
#' @return A dataframe with the number and proportion of genotyped individuals
#'   with a female and male parent assigned. The assigned parents are written
#'   to a text file, by default `ParentsAssigned.txt'.
#'   In addition, it will create a text file with parent-offspring pairs where
#'   it is unclear which is the parent and which is the offspring, as either
#'   individual does not have a birthyear, and inferred sex for individuals
#'   with no sex provided, but which could be assigned as parent.
#'
#' @useDynLib sequoia

SeqParents <- function() {
  TMP <- .Fortran("parents",
                  nParents = integer(2),
                  nAmbiguous = as.integer(0),
                  TotLik = double(42),
                  PACKAGE = "sequoia")
  Specs <- ReadSpecs()
  Assigned <- c(Mothers = TMP$nParents[1],
              Fathers = TMP$nParents[2])
  nInd <- FacToNum(Specs["NumberIndivGenotyped",])
  Assigned <- data.frame(Assigned, Prop = round(Assigned / nInd,2))

  if (TMP$nAmbiguous>0) {
     message("There are  ", TMP$nAmbiguous, "  potential parent-offspring pairs \n which are
             not assigned in the pedigree, possibly due to unknown birth year(s).\n",
             "Please see 'UnassignedParents'.")
  }
  list(NumberParents = Assigned,
       ParentsTotLik = TMP$TotLik[1:sum(TMP$TotLik!=0)])
}


#=====================================================================
#=====================================================================

#' Sibship clustering
#'
#' Create clusters of full- and half- siblings with a dummy parent,
#' assign grandparents to such sibships, and assign parental relationships
#' between dummies. Parameter values and filenames are taken from the
#' SequioiaSpecs text file, generated by \code{\link{SeqPrep}}.
#'
#'
#' The output file, by default `PedSeq.txt', contains the following columns:
#' \describe{
#'   \item{id}{Individual ID}
#'   \item{dam}{Assigned mother, or NA}
#'   \item{sire}{Assigned father, or NA}
#'   \item{LLR_dam}{Log-Likelihood Ratio (LLR) of this female being the mother,
#'      or the individual belonging to this sibship, versus the next most
#'      likely relationship between the focal individual and this female or
#'      sibship (see ... for relationships considered)}
#'   \item{LLR_sire}{idem, for male (dummy) parent}
#'   \item{LLR_pair}{trio LLR, versus the next most likely configuration
#'      (containing one or neither (dummy) parent)}
#' }
#'
#' @return A dataframe with the number and proportion of individuals with a
#' female and male genotyped or dummy parent assigned. The pedigree is written
#' to a text file, by default `PedSeq.txt'.
#'
#' @useDynLib sequoia

SeqSibs <- function() {
  TMP <- .Fortran("Sibships",
                  TotLik = double(42),
                  PACKAGE = "sequoia")
  list(SibTotLik = TMP$TotLik[1:sum(TMP$TotLik!=0)])
}


#################################################################


.onUnload <- function (libpath) {
  library.dynam.unload("sequoia", libpath)
}
