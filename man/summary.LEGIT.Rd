% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LEGIT.R
\docType{data}
\name{summary.LEGIT}
\alias{summary.LEGIT}
\title{Summarizing LEGIT fits}
\format{An object of class \code{function} of length 1.}
\usage{
summary.LEGIT
}
\arguments{
\item{object}{An object of class "LEGIT", usually, a result of a call to LEGIT.}

\item{...}{Further arguments passed to or from other methods.}
}
\value{
Returns a list of objects of class "summary.glm" containing the summary of each parts (main, genetic, environmental) of the model.
}
\description{
Shows the summary for all parts (main, genetic, environmental) of the LEGIT model.
}
\examples{
	train = example_2way(250, 1, seed=777)
fit_default = LEGIT(train$data, train$G, train$E, y ~ G*E)
summary(fit_default)
}
\keyword{datasets}
