% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{BBB_NMP_pl_subset}
\alias{BBB_NMP_pl_subset}
\title{Example for "pl" format specification}
\format{
A data frame with 1,000 rows and 35 variables:
\describe{
  \item{chr}{The chromosome number}
  \item{marker_name}{The marker id}
  \item{pos}{The position of the marker}
  \item{allele1}{The name of the first marker allele}
  \item{allele2}{The name of the second marker allele}
  \item{id1_pl1}{The AA phred likelihood for the first individual}
  \item{id1_pl2}{The AB phred likelihood for the first individual}
  \item{id1_pl3}{The BB phred likelihood for the first individual}
  \item{id2_pl1}{The AA phred likelihood for the second individual}
  \item{id2_pl2}{The AB phred likelihood for the second individual}
  \item{id2_pl3}{The BB phred likelihood for the second individual}
  \item{id3_pl1, id3_pl2, id3_pl3, id4_pl1, id4_pl2, id4_pl3, id5_pl1, id5_pl2, id5_pl3, id6_pl1, id6_pl2, id6_pl3, id7_pl1, id7_pl2, id7_pl3, id8_pl1, id8_pl2, id8_pl3, id9_pl1, id9_pl2, id9_pl3, id10_pl1, id10_pl2, id10_pl3}{The phred likelihoods for the other individuals} }
}
\usage{
BBB_NMP_pl_subset
}
\description{
A dataset containing real genotype likelihoods in phred scores for 1,000
SNPs. The data is available for ten individuals and the first 1,000 SNPs on
chromosome 1. The data corresponds to a low-fold sequencing experiment. There
are three columns per individuals (for genotypes 00, 01 and 11).
}
\keyword{datasets}
