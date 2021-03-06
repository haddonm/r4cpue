% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpuefuncs.r
\name{coef.CEout}
\alias{coef.CEout}
\title{coef.CEout - S3 method to extract the parameter values from a CEout object}
\usage{
\method{coef}{CEout}(inout)
}
\arguments{
\item{inout}{is a CEout object produced by standLnCE}
}
\value{
a matrix containing the optimum model parameters, including the
  log-space coefficients, the lower and upper 95% cofidence intervals,
  the standard error of each parameter, and p the probability that the
  parameter is significantly different from zero. These coefficients
  would need to be back-transformed to return to the original scale. This
  would use x = exp(coef + (sterr*sterr)/2); note the bias correction
}
\description{
Extract the parameter values from a CEout object from standLnCE
  by extending the generic 'coef' function to apply to CEout objects, such
  as generated by standLnCE
}
\examples{
\dontrun{
data(sps)
splabel = "Species"
labelM <- c("Year","Vessel","Month")
mods <- makemodels(labelM)
sps1 <- makecategorical(labelM,sps)
out <- standLM(mods,sps1,splabel)
coef(out)
}
}
