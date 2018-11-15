#' Example from a small simulated data set
#'
#' A dataset containing simulated genotypes for 20 individuals.
#' Genotypes are available for 10,000 SNPs on 10 chromosomes.
#' The ten last columns correspond to the genotypes.
#'
#' @format A data frame with 10,000 rows and 24 variables:
#' \describe{
#' \item{chr}{The chromosome number}
#' \item{pos}{The position of the marker}
#' \item{allele1}{The name of the first marker allele}
#' \item{allele2}{The name of the second marker allele}
#' \item{id1}{The genotype for the first individual}
#' \item{id2}{The genotype for the second individual}
#' \item{id3, id4, id5, id6, id7, id8, id9, id10, id11, id12, id13, id14, id15, id16, id17, id18, id19, id20}{The genotypes of
#' the remaining individuals}
#' }
"genosim"

#' Subset of a dataset with genotypes for 20 sheeps
#'
#' A dataset containing real genotypes for 20 individuals.
#' Genotypes are available for 14,831 SNPs from the first three chromosomes
#' were selected. The twenty last columns correspond to the genotypes.
#' Missing genotypes are set to 9.
#'
#' @format A data frame with 14,831 rows and 25 variables:
#' \describe{
#' \item{chr}{The chromosome number}
#' \item{marker_name}{The name of the marker}
#' \item{pos}{The position of the marker}
#' \item{allele1}{The name of the first marker allele}
#' \item{allele2}{The name of the second marker allele}
#' \item{id1}{The genotype for the first individual}
#' \item{id2}{The genotype for the second individual}
#' \item{id3, id4, id5, id6, id7, id8, id9, id10, id11, id12, id13, id14, id15, id16, id17, id18, id19, id20}{The genotypes of
#' the remaining individuals}
#' }
"genoex"

#' Subset of a dataset with genotypes for 6 individuals from a cattle population.
#'
#' A dataset containing real genotypes for 6 individuals.
#' Genotypes are available for a low density array with 6370 SNPs on 29 autosomes.
#' The six last columns correspond to the genotypes.
#' Missing genotypes are set to 9.
#'
#' @format A data frame with 6370 rows and 10 variables:
#' \describe{
#' \item{chr}{The chromosome number}
#' \item{pos}{The position of the marker}
#' \item{allele1}{The name of the first marker allele}
#' \item{allele2}{The name of the second marker allele}
#' \item{id1}{The genotypes for the first individuals (id1)}
#' \item{id2}{The genotypes for the second individuals (id2)}
#' \item{id3}{The genotypes for the third individuals (id3)}
#' \item{id4}{The genotypes for the fourth individuals (id4)}
#' \item{id5}{The genotypes for the fifth individuals (id5)}
#' \item{id6}{The genotypes for the sixth individuals (id6)}
#' }
"typs"

#' A file with marker allele frequencies for the cattle population.
#'
#' The allele frequencies of the first allele were computed in a larger sample.
#'
#' @format A data frame with 6370 rows and one variable.
#' \describe{
#' \item{frq}{The allele frequencies estimated for the first allele}
#' }
"typsfrq"

#' Example for "ad" format specification
#'
#' A dataset containing real allele depth information for 1,000 SNPs.
#' The data is available for ten individuals and the first 1,000 SNPs on chromosome 1.
#' The data corresponds to a low-fold sequencing experiment.
#' There are two columns per individuals (read counts for allele 1 and for allele 2).
#'
#' @format A data frame with 1,000 rows and 25 variables:
#' \describe{
#' \item{chr}{The chromosome number}
#' \item{marker_name}{The marker id}
#' \item{pos}{The position of the marker}
#' \item{allele1}{The name of the first marker allele}
#' \item{allele2}{The name of the second marker allele}
#' \item{id1_ad1}{The read count for the first individual at the first markrer}
#' \item{id1_ad2}{The read count for the first individual at the second markrer}
#' \item{id2_ad1}{The read count for the second individual at the first markrer}
#' \item{id2_ad2}{The read count for the second individual at the second markrer}
#' \item{id3_ad1, id3_ad2, id4_ad1, id4_ad2, id5_ad1, id5_ad2, id6_ad1, id6_ad2, id7_ad1, id7_ad2, id8_ad1, id8_ad2, id9_ad1, id9_ad2, id10_ad1, id10_ad2}{The
#' read counts for the other individuals}
#' }
"BBB_NMP_ad_subset"

#' Example for "gp" format specification
#'
#' A dataset containing real genotype probabilities for 1,000 SNPs.
#' The data is available for ten individuals and 1,000 SNPs on chromosome 10.
#' The data corresponds to a low-fold sequencing experiment.
#' There are three columns per individuals (for genotypes 00, 01 and 11).
#'
#' @format A data frame with 1,000 rows and 35 variables:
#' \describe{
#' \item{chr}{The chromosome number}
#' \item{marker_name}{The marker id}
#' \item{pos}{The position of the marker}
#' \item{allele1}{The name of the first marker allele}
#' \item{allele2}{The name of the second marker allele}
#' \item{id1_gp1}{The AA genotype probability for the first individual}
#' \item{id1_gp2}{The AB genotype probability for the first individual}
#' \item{id1_gp3}{The BB genotype probability for the first individual}
#' \item{id2_gp1}{The AA genotype probability for the second individual}
#' \item{id2_gp2}{The AB genotype probability for the second individual}
#' \item{id2_gp3}{The BB genotype probability for the second individual}
#' \item{id3_gp1, id3_gp2, id3_gp3, id4_gp1, id4_gp2, id4_gp3, id5_gp1, id5_gp2, id5_gp3, id6_gp1, id6_gp2, id6_gp3, id7_gp1, id7_gp2, id7_gp3, id8_gp1, id8_gp2, id8_gp3, id9_gp1, id9_gp2, id9_gp3, id10_gp1, id10_gp2, id10_gp3}{The genotype probabilities for
#' the other individuals}
#' }
"BBB_NMP_GP_subset"

#' Example for "pl" format specification
#'
#' A dataset containing real genotype likelihoods in phred scores for 1,000
#' SNPs. The data is available for ten individuals and the first 1,000 SNPs on
#' chromosome 1. The data corresponds to a low-fold sequencing experiment. There
#' are three columns per individuals (for genotypes 00, 01 and 11).
#'
#' @format A data frame with 1,000 rows and 35 variables: \describe{
#'   \item{chr}{The chromosome number} \item{marker_name}{The marker id}
#'   \item{pos}{The position of the marker} \item{allele1}{The name of the first
#'   marker allele} \item{allele2}{The name of the second marker allele}
#'   \item{id1_pl1}{The AA phred likelihood for the first individual}
#'   \item{id1_pl2}{The AB phred likelihood for the first individual}
#'   \item{id1_pl3}{The BB phred likelihood for the first individual}
#'   \item{id2_pl1}{The AA phred likelihood for the second individual}
#'   \item{id2_pl2}{The AB phred likelihood for the second individual}
#'   \item{id2_pl3}{The BB phred likelihood for the second individual}
#'   \item{id3_pl1, id3_pl2, id3_pl3, id4_pl1, id4_pl2, id4_pl3, id5_pl1, id5_pl2, id5_pl3, id6_pl1, id6_pl2, id6_pl3, id7_pl1, id7_pl2, id7_pl3, id8_pl1, id8_pl2, id8_pl3, id9_pl1, id9_pl2, id9_pl3, id10_pl1, id10_pl2, id10_pl3}{The phred likelihoods for the other individuals} }
"BBB_NMP_pl_subset"

#' Example for "gt" format specification
#'
#' A dataset containing real genotypes for 1,000 SNPs.
#' The data is available for ten individuals and the first 1,000 SNPs on chromosome 1.
#' The data corresponds to a WGS experiment.
#' There is one columns per individual.
#'
#' @format A data frame with 1,000 rows and 14 variables:
#' \describe{
#' \item{chr}{The chromosome number}
#' \item{pos}{The position of the marker}
#' \item{allele1}{The name of the first marker allele}
#' \item{allele2}{The name of the second marker allele}
#' \item{id1}{The genotypes for the first individuals (id1)}
#' \item{id2}{The genotypes for the second individuals (id2)}
#' \item{id3, id4, id5, id6, id7, id8, id9, id10}{The genotypes of
#' the remaining individuals}

#' }
"BBB_PE_gt_subset"

#' A file with names or IDs for ten samples.
#'
#' The names (or IDs) are provided for ten samples.
#'
#' @format A data frame with 10 rows and one variable.
#' \describe{
#' \item{IDs}{The ids for ten individuals}
#' }
"BBB_samples"

#' The result of an analysis on 110 sheeps from the Soay population.
#'
#' The results were obtained by running the default model (10 classes
#'  with pre-defined rates) on 110 individuals genotyped at 37465 SNPs.
#'
#' @format the results are a zres object.
"soay_mix10r"

#' The result of an analysis on 22 sheeps from Rasa Aragonesa population.
#'
#' The results were obtained by running the default model (10 classes
#'  with pre-defined rates) on 22 individuals genotyped at 37465 SNPs.
#'
#' @format the results are a zres object.
"rara_mix10r"

#' The result of an analysis on 23 sheeps from Wiltshire population.
#'
#' The results were obtained by running the default model (10 classes
#'  with pre-defined rates) on 23 individuals genotyped at 37465 SNPs.
#'
#' @format the results are a zres object.
"wilt_mix10r"


