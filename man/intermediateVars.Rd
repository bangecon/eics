% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intermediateVars.R
\name{intermediateVars}
\alias{intermediateVars}
\title{Checking for collider bias using random forests.}
\usage{
intermediateVars(
  data,
  treatment,
  controls,
  ntree = 100,
  zthresh = 1,
  out = FALSE,
  parallel = FALSE
)
}
\arguments{
\item{data}{is the universe of all variables considered, which may include
both controls and excluded variables (e.g. variables eliminated by a model
selection algorithm).}

\item{treatment}{is the name of the treatment variable contained in
\code{data}.}

\item{controls}{is a list() of names of selected controls contained in
\code{data}.}

\item{ntree}{is the numerical number of trees to use in the random forest
models.}

\item{zthresh}{is the z-score used as a threshold for flagging a variable
for suspicion as an alternative treatment, feedback, or intermediate bias.}

\item{out}{is a logical indicating whether to create a file named
"alternateTreatments.csv" in the current working directory that contains
the variables flagged as alternate treatments.}

\item{parallel}{is a logical indicating whether to use parallel processing}
}
\value{
\code{intermediateVars} returns a list that contains the following
components:
\item{treatment}{name of the treatment variable.}
\item{controls}{the list of control variables.}
\item{intermediateVariables}{the list of possible intermediate variables.}
\item{zThresh}{the z-score threshold for suspicion.}
\item{zMatrix}{the scaled variable importance matrix between the treatment,
controls, and outcome.}
}
\description{
This function checks whether control variables that a process (e.g.
researcher "intuition", a model-selection algorithm) selected might create
collider bias with respect to the current treatment.
}
