% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/g.R
\name{g}
\alias{g}
\title{moment function generator}
\usage{
g(type)
}
\arguments{
\item{type}{the type of distribution. It take a value from "Poisson", "power law", "Gamma", "Beta", mixture of "2Poissons", and "2Exponential"}
}
\description{
If we need to estimate n parameters, we can use Moment Generating Function to get from E[x] to E[X^n] in order to obtain n equations
The value meanb is our estimate in the currently iteration of the mean of the beta distribution,
while x is our data. The gmm() function, in calling our g(), will calculate the average of m1 (m2, m3) over
all our data, and then in setting that average to 0, we are equating our parametric estimate of the
population mean to our sample mean, exactly what MM is supposed to do.
}

