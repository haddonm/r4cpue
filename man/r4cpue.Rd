% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpue_doc.R
\docType{package}
\name{r4cpue}
\alias{r4cpue}
\title{r4cpue a set of functions to assist with CPUE standardization}
\description{
The r4cpue package provides two categories of functions
    analytical functions that assist with the cpue standardization
    plotting functions that illustrate the results of the standardizations.
    In addition there are example data sets
}
\section{Analytical functions}{

\describe{
  \item{coef.CEout}{S3 function applied to the output of standLM}
  \item{dosingle}{conducts a standardization of indat using the
      inmodel, that can be generated using makeonemodel}
  \item{geomean}{Calculates the bias corrected geometric mean}
  \item{getmeans}{Calculates a set of central tendency measures for
      nominal CPUE}
  \item{makecategorical}{converts a series of variables in a data.frame into
      factors}
  \item{makemodels}{Given a list of factors this generates a list of
      formula for inclusion in glm or standLM}
  \item{makeonemodel}{generates a single model for use in dosingle}
  \item{scaleCE}{Rescales a vector of CPUE to a mean of 1.0 or of avCE}
  \item{standLM}{Uses lm to standardize log-transformed CPUE}
  \item{summary.CEout}{an S3 function to summarize a standardization output}
}
}

\section{Plotting functions}{

\describe{
  \item{addnorm}{fits a normal distribution curve to a given histogram}
  \item{diagnosticPlot}{plots some diagnostic details for a standardization}
  \item{impactplot}{Plots the influence of each factor}
  \item{inthist}{a replacement for the hist function for use with integers}
  \item{plotdata}{plots graphs of untransformed and log-transformed data}
  \item{plotprep}{defines a base graphics window for use in RStudio}
  \item{plotstand}{Plots the optimum model vs the year-only model}
  \item{plotstandFY}{The same as plotstand but for Fishing year species}
  \item{yearBubble}{Generates a bubbleplot of x against Year}
}
}

\section{Utility functions}{

\describe{
   \item{addcount}{adds a new column to input data.frame that is a count
       of the number of years in which the identified variable, default
       'Vessel' occurs each year for each level of the factor}
   \item{fishery}{generates vectors of year, catch, effort, and cpue}
   \item{getfact}{extracts a given factor from the analysis with its
       standard errors and rescaled to a mean of 1.0}
   \item{getStand}{extracts the main year paramters from a standLM object
       with its StErr and confidence itervals}
   \item{properties}{Checks a data.frame for NAs and counts; used for QC}
   \item{removeEmpty}{removes empty strings from a vector of strings}
   \item{selectdata}{simplifies the selection of data from a data.frame
       by depth, years, zones, method, and fishery}
   \item{toExcel}{copys the selected vector or matrix to the clipboard
       it can be pasted directly into Excel or other software}
   \item{yearNA}{counts the NAs in each numeric field in a data.frame}
   \item{yearZero}{counts the zeros and NAs in identified fields in a df}
}
}

