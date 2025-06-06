% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zooroh_accessor.R
\name{update_zres}
\alias{update_zres}
\title{Update one main zres object with new results}
\usage{
update_zres(zres1, zres2)
}
\arguments{
\item{zres1}{The main zres object that will be modified.}

\item{zres2}{The new zres object, with the new results that will be inserted in the
main zres object.}
}
\value{
an updated zres object containing the results from zres1 updated by those from
zres2.
}
\description{
The function is used for example when the main analysis failed for one individual.
The analysis is repeated for that individual with other parameters. The function
can then be used to insert the new results in the main zres object. To avoid
generating to large files, the updated zres object can take the same name as zres1
doing zres1 <- update_zres(zres1,zres2).
}
