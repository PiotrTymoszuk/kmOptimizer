# Non-exported utils

#' Check if the string needs non-standard name escape.
#'
#' @description Checks if a string is a syntactically valid name.
#' @param x a string vector.
#' @return a logical value.

  isValidName <- function(x) {

    grepl("^((([[:alpha:]]|[.][._[:alpha:]])[._[:alnum:]]*)|[.])$", x)

  }

