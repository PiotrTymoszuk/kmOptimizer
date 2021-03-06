% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{find_cutoff}
\alias{find_cutoff}
\title{Find optimal survival cutoff.}
\usage{
find_cutoff(data, time, event, variable, min_n = 2, .parallel = FALSE, ...)
}
\arguments{
\item{data}{a data frame with the survival time, a numeric event index
variable and the variable to be cut.}

\item{time}{name of the survival time variables, needs to be numeric.}

\item{event}{name of the binary event variable: 1 is event, 0 no event.}

\item{variable}{name of the variable to be cut, needs to be numeric.}

\item{min_n}{the minimal strata size.}

\item{.parallel}{logical, should the search be run in parallel?}

\item{...}{additional arguments passed to \code{\link[survival]{survdiff}}.}
}
\value{
an object of the 'survcut'. Contains te optimal cutoff value
('cutoff'), a tibble with the numeric variable of interest and and it's
strata and a tibble with the cutoff values, test statistics
and p values for the survival differences.
}
\description{
Finds the optimal cutoff of a variable corresponding to the
greatest difference in survival. The difference in survival is estimated with
\code{\link[survival]{survdiff}}.
}
\details{
In case, more minuma are find, the stratified data table is
adjusted accordingly.
}
