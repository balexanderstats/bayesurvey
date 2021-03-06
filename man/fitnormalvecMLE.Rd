% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation.R
\name{fitnormalvecMLE}
\alias{fitnormalvecMLE}
\title{MLE estimates of a Normal Distribution
Fits a normal}
\usage{
fitnormalvecMLE(data, logit = F)
}
\arguments{
\item{data}{- a vector of data}

\item{logit}{- if a logit transformation is needed}
}
\value{
a list with the elements: mean with the mean of the distribution and var with the variance of the distribution.
}
\description{
MLE estimates of a Normal Distribution
Fits a normal
}
\examples{
set.seed(1)
data1  = rnorm(100)
fitnormalvecMLE(data1)
data2 = rnorm(100, mean  = .5, sd = 0.03)
fitnormalvecMLE(data2, logit = TRUE)
}
