% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zplot_prophbd.R
\name{zooplot_prophbd}
\alias{zooplot_prophbd}
\title{Plot proportion of the genome associated with different HBD classes}
\usage{
zooplot_prophbd(input, cols = NULL, style = "barplot", cumulative = FALSE)
}
\arguments{
\item{input}{a named list with one or several zres objects obtained after
running zoorun. The zres objects are the output of the zoorun function. For
instance, putting list(name1 = zres1, name2 = zres2). The function will then
use the names in the plot (in case several zres objects are used).}

\item{cols}{a vector with the colors to be used for each population or zres
object.}

\item{style}{select "barplot", "lines" or "boxplot" for the graphic styles. Boxplot
can be used with a single zres file or population.}

\item{cumulative}{a logical indicating whether mean autozygosity is estimated
per class (FALSE) or summed over all HBD class with a rate smaller than a
value (these cumulated values are obtained for every rate defined in the
model). By default, this value is FALSE. When FALSE, the percentages correspond
to the mean individual genome-wide probabilities of belonging to each HBD-class
or to the fraction of the genome in an autozygosity class. When TRUE, we obtain
the mean probability of belonging to an HBD class with a rate smaller or equal than
a threshold (here we use the pre-defined rates of the model as thresholds), averaged
over the whole genome and all individuals. This corresponds to report mean genomic
inbreeding coefficients estimated with respect to different base populations obtained
by selecting different thresholds T that determine which HBD classes are considered
in the estimation of the genomic inbreeding coefficient (setting the base population
approximately 0.5 * T generations ago).}
}
\value{
The function plots either the average proportion of the genome associated with
different HBD classes or the average genomic inbreeding coefficient estimated with respect
to different base populations (from young to older).
}
\description{
Plot the mean percentage of the genome in different HBD classes or the
inbreeding coefficient obtained by summing autozygosity associated with HBD
classes with a rate lower or equal to a threshold (e.g., including all HBD
classes with longer and more recent HBD segments than a selected threshold).
}
