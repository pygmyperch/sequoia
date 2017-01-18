# Various miscelaneous functions

#======================================================================
#' Wrapper for suppressMessages
#'
#' Provide additional on/off argument to suppressMessages()
#'
#' @param X a function
#' @param quiet  suppress messages / not

supprOF <- function(X, quiet=FALSE) if(quiet) suppressMessages(X) else X


#======================================================================
#' Convert factor to numeric
#'
#' Converts factors to their numeric values, via \code{as.character}.
#'
#' @param x a factor
#' @return A numeric vector of the same length as x.
FacToNum <- function(x) as.numeric(as.character(x))


#======================================================================
#' Data input
#'
#' Read data with header row.
#'
#' @param ... parameters for read.table
#' @param sep column seperator
ReadTable <- function(..., sep="\t") utils::read.table(...,
                                                header=TRUE,
                                                stringsAsFactors=FALSE,
                                                sep=sep,
                                                na.strings=c("", NA))


#======================================================================

#' Read specs file
#'
#' Read the sequoia specs text file.
#'
#' @param SpecsFile comma-separated text file with parameter settings,
#'   containing 1 column with labels and 1 column with values
ReadSpecs <- function(SpecsFile = "SequoiaSpecs.txt") {
  utils::read.table(SpecsFile,
             row.names = 1, header = FALSE, sep = ",",
             strip.white = TRUE, stringsAsFactors = FALSE)
}


#======================================================================
#' table
#'
#' sets UseNA to 'ifany'.
#'
#' @param ... one or more objects which can be interpreted as factors
#'  (including character strings), or a list (or data frame) whose components
#'   can be so interpreted.
Table <- function(...) table(..., useNA="ifany")


#======================================================================
#' Value Matching
#'
#' Like \code{\%in\%}, returns a logical vector indicating if there is a match
#'  or not for its left operand, but returns NA for NA's in the left operand.
#'
#' @param x vector: the values to be matched
#' @param y vector: the values to be matched against
#'
#' @return A logical vector of the same length as x.
#'
# #' @examples
# #' X <- c(1:5, NA, NA)
# #' Y <- c(3:10)
# #' X %in% Y
# #' table(X %in% Y, useNA="ifany")
# #' X %ina% Y
# #' table(X %ina% Y, useNA="ifany")
"%ina%" <- function(x, y) ifelse(!is.na(x), match(x, y, nomatch = 0) > 0, NA)


#======================================================================
#' paste directory and file name
#'
#' @param folder foldername
#' @param fileName filename

pasteD <- function(folder, fileName) paste(folder, fileName, sep="/")


#======================================================================
#' create a table, and ensure that the levels TRUE, FALSE and NA are always all
#'  represented
#'
#' @param x  a logical vector

tbl.logic <- function(x) table(factor(x, levels=c(TRUE, FALSE, NA)),
                               useNA="always")


#======================================================================
#' Comparison
#'
#' Test which elements in a vector are equal to x, returning FALSE at missing
#' values in V
#'
#' @param x  a value
#' @param V  a vector of the same type as x
#'
#' @return a logical vector of the same length as v
#'
eqv <- function(x, V) {
  if (!is.na(x)) {
    y <- ifelse(!is.na(V), x==V, FALSE)
  } else {
    y <- logical(length=length(V))
  }
  y
}


#======================================================================
#' function in Examples from integer {base}
#'
#' @param x a number
#' @param tol tolerance
is.wholenumber <-
  function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
