#'Extracts the realized autozygosity from the zres object
#'
#'Extracts the realized autozygosity from the zres object. Extraction is
#'performed for the indicated classes (all by default) and names are added to
#'the columns. The function must be used with more than one individual is the
#'zres object.
#'
#'@param zres The name of the zres object created by the zoorun function.
#'
#'@param classNum An array with the number of the classes to extract. All
#'  classes are extracted by default.
#'
#'@return The function returns a data frame with one row per individual and one
#'  column per extracted classes. In addition, it gives names to the columns.
#'  For a pre-defined model, the names of HBD classes are "R_X" where X is the
#'  rate of the corresponding class. For a model with rate estimation, the names
#'  of the HBD classes are "HBDclassX" where X is the number of the HBD class.
#'  For non-HBD classes, we use "NonHBD".
#'
#'@export

realized <- function(zres, classNum = NULL) {
  K = length(zres@krates[1,])
  if(is.null(classNum)){
    classNum = seq(1,K)
  }
  matout <- as.data.frame(zres@realized[,classNum])
  K2 <- length(classNum)
  classNames <- array('',K2)
  for (i in 1:K2){
    if(sd(zres@krates[,K])!=0){
      classNames[i]=paste("HBDclass",classNum[i],sep="")
      if(classNum[i]==K){classNames[i]="NonHBD"}
    }
    if(sd(zres@krates[,K])==0){
      classNames[i]=paste("R_",zres@krates[1,classNum[i]],sep="")
      if(classNum[i]==K){classNames[i]="NonHBD"}
    }
  }
  colnames(matout) <- classNames
  return(matout)
}

#'Computes the realized inbreeding coefficient
#'
#'Computes the realized inbreeding coefficient with respect to a base population
#'by summing the autozygosity from all HBD class with a rate lower or equal to the threshold T.
#'This amounts to set the base population approximately 0.5*T generations ago.
#'HBD classes with a higher rate are then no longer considered as autozygous.
#'
#'@param zres The name of the zres object created by the zoorun function.
#'
#'@param T The value chosen to define the base population. When T is not provided, all HBD classes
#'are considered to estimate the inbreeding coefficient.
#'
#'@return An array with the compute inbreeding coefficients for all the individuals in the analysis.
#'
#'@export

cumhbd <- function(zres, T = NULL) {
  K = length(zres@krates[1,])
  if(sd(zres@krates[,K])!=0){
    print(c("cumhbd is not recommended for a model with different rates per individual!"))
  }
  nind = length(zres@ids)
  f <- array(0,nind)
  if(K > 2){x <- t(apply(zres@realized[,1:(K-1)],1,cumsum))}
  if(K == 2){x <- zres@realized[,1]}
  if(is.null(T)){
    if(K > 2){f = x[,(K-1)]}
    if(K == 2){f =x}
  }else{
    T2=array(0,nind)
    for (id in 1:nind){
      for (j in 1:(K-1)){
        rk = zres@krates[id,j]
        if(rk <= T){T2[id] = j}
      }
      if(T2[id] > 0){
        if(K > 2){f[id] = x[id,T2[id]]}
        if(K == 2){f[id] = x[id]}
        }
    }
  }

  return(f)
}

#'Extracts the HBD segments from the zres object
#'
#'Extracts the HBD segments (or RoH) from the zres object. Extraction is performed
#'for the indicated individuals and the selected region (all by default).
#'
#'@param zres The name of the zres object created by the zoorun function.
#'
#'@param ids An array with the ids of the individuals to extract. All individuals
#'are extracted by default.
#'
#'@param chrom the number of the chromosome where we are looking for HBD segments.
#'  This chromosome number refers to the position of the chromosome in the list
#'  of all chromsomes present in the input genotype data.
#'
#'@param startPos The starting position (on the chromosome) of the interval from which we will extract HBD segments
#'(1 by default).
#'
#'@param endPos The ending position (on the chromosome) of the interval from which we will extract HBD segments
#'(last position by default).
#'
#'@param inside A logical indicating whether we extract only segment within the interval (TRUE) or overlapping with the
#'interval (FALSE). By 'within the interval', we mean that both starting and end position of the HBD segment should be
#'in the interval. By 'overlapping', we mean that at least one part of the HBD segment should be located in the
#'interval.
#'
#'@return The function returns a data frame  with the HBD segments fitting the filtering rules (id and position).
#'The data frame has one line per identified HBD segment and nine columns: id is the number of the individual in
#'which the HBD segments is located, chrom is the chromosome of the HBD segments, start_snp is the number of the SNP
#' at which the HBD segment starts (the SNP number within the chromosome), start_end is the number of the SNP at
#' which the HBD segment ends (the SNP number within the chromosome), start_pos is the position at which the HBD
#' segment starts (within the chromosome), start_end is the position at which the HBD segment ends (within the
#' chromosome), number_snp is the number of consecutive SNPs in the HBD segment, length is the length of the HBD
#' segment (for instance in bp or in cM/1000000) and HBDclass is the HBD class associated with the HBD segment.
#'
#'@export

rohbd <- function(zres, ids = NULL, chrom = NULL, startPos = NULL, endPos = NULL, inside = TRUE) {
  if(is.null(ids)){ids <- zres@ids}
  if(is.null(chrom)){
    roh <- zres@hbdseg[zres@hbdseg$id %in% ids,]
  } else {
    if(is.null(startPos)){startPos = 0}
    if(is.null(endPos)){endPos = max(zres@hbdseg$end_pos)}
    if(inside){
      roh <- zres@hbdseg[zres@hbdseg$id %in% ids & zres@hbdseg$chrom == chrom
                         & zres@hbdseg$start_pos >= startPos & zres@hbdseg$end_pos <= endPos,]
    } else {
      roh <- zres@hbdseg[zres@hbdseg$id %in% ids & zres@hbdseg$chrom == chrom
                         & zres@hbdseg$start_pos <= endPos & zres@hbdseg$end_pos >= startPos,]
    }
  }
  return(roh)
}

#'Extracts the HBD probabilities from the zres object
#'
#'Extracts the locus specific HBD probabilities for an individual. A specific chromosomal region can be specified.
#'A threshold T can be used to determine which HBD classes are used in the computation of the HBD probability. This
#'function requires that the option "localhbd" was set to TRUE when creating the zres object.
#'
#'@param zres The name of the zres object created by the zoorun function.
#'
#'@param zooin The name of the zdata object created by the zoodata function. See
#'"zoodata" for more details.
#'
#'@param id The number of the individual to extract.
#'
#'@param chrom the number of the chromosome where we are looking for HBD segments.
#'  This chromosome number refers to the position of the chromosome in the list
#'  of all chromsomes present in the input genotype data.
#'
#'@param startPos The starting position (on the chromosome) of the interval from which we will extract HBD segments
#'(1 by default).
#'
#'@param endPos The ending position (on the chromosome) of the interval from which we will extract HBD segments
#'(last position by default).
#'
#'@param T The value chosen to define the base population (to determine which classes are used in estimated of HBD
#'probability, which classes are considered autozygous). When T is not provided, all HBD classes
#'are considered to estimate the local HBD probability.

#'
#'@return The function returns a vector of HBD probabilities for the specified individual and chromosomal region.
#'The HBD probabilities are computed as the sum of the probabilities for each HBD class with a rate smaller or equal than
#'the threshold (the sum from all the HBD classes when T is not specified).
#'
#'@export

probhbd <- function(zres, zooin, id, chrom = NULL, startPos = NULL, endPos = NULL, T = FALSE) {
  hbdp <- array(0,zooin@nsnps)
  K=length(zres@krates[id,]);K2=K-1
  if(T){
    K2=0
    for (k in 1:(K-1)){if(zres@krates[id,k] <= T){K2=k}}
  }
  if(is.null(chrom)){
    if( K2 < (K-1)){
      length(hbdp)
      length(zres@hbdp[[id]][1,])
      for (k in 1:K2){hbdp = hbdp + zres@hbdp[[id]][k,]}
    } else {
      hbdp = 1 - zres@hbdp[[id]][K,]
    }
  } else {
    if(is.null(startPos)){startPos = 0}
    if(is.null(endPos)){endPos = max(zres@hbdseg$end_pos)}
     ff <- (zooin@bp >= startPos & zooin@bp <= endPos &
              seq(1,zooin@nsnps) >= zooin@chrbound[chrom,1] &
              seq(1,zooin@nsnps) <= zooin@chrbound[chrom,2])
     hbdp = hbdp[ff]
     if( K2 < (K-1)){
        for (k in 1:K2){hbdp = hbdp + zres@hbdp[[id]][k,ff]}
     } else {
       hbdp = 1- zres@hbdp[[id]][K,ff]
     }
  }
  return(hbdp)
}

#'Returns number of HBD segments (or genome proportion) per size bins for DoRIS.
#'
#'The number of HBD segments in 1Mb bins is computed and the function returns a
#'table that can be used with DoRIS to estimate paste effective population size.
#'Alternatively, the proportion of the genome in the different 1 Mb bins is
#'computed. The user must specify the range of the bins (the smallest and largest
#'HBD segments considered).
#'#'
#'@param zres The name of the zres object created by the zoorun function.
#'
#'@param minv The minimum length of HBD segments (in Mb) used for the DoRIS analysis.
#'
#'@param maxv The maximum length of HBD segments (in Mb) used for the DoRIS analysis.
#'
#'@param glen The length of the genome in Mb (including only the autosomes, the space where
#'HBD segments have been searched for).
#'
#'@param method The argument is used to indicate that function will return counts of
#'HBD segments (method = "counts") or proportion of the genome (method ="sharing") per
#'1 Mb bins.
#'
#'@return A data frame with three columns and one row per 1Mb bins. The first column
#'indicates the the length of smallest HBD segments in the bin, the second column
#'indicates the length of the largest HBD segments in the bin and the third column
#'indicates the number of segments in the bin or the proportion of the genome
#'in the bin.
#'
#'@export

zoodoris <- function(zres, minv, maxv, glen, method='counts') {

win_st <- seq(minv,(maxv-1))
win_end <- seq((minv+1),maxv)

nseg = array(0,length(win_st))
sharing = array(0,length(win_st))
npairs = zres@nind

for (i in 1:length(nseg)){
sel <- zres@hbdseg[zres@hbdseg$length/1000000 > win_st[i] &
                     zres@hbdseg$length/1000000 <= win_end[i],]
nseg[i] <- dim(sel)[1]
sharing[i] <- sum(sel$length/1000000) /(glen*npairs)
}

out1 <- data.frame(cbind(win_st, win_end, nseg))
out2 <- data.frame(cbind(win_st, win_end, sharing))

if(method=='counts'){return(out1)}
if(method=='sharing'){return(out2)}

}

#'Realizes a simulation using a zooroh data set and a zooroh model
#'
#'Performs a simulation under a model similar too ZooRoH. It simulates the
#'genome as a mosaic of HBD and non-HBD segments. Several non-HBD classes can be
#'simulated. Classes with high rates (short segments from ancient ancestors) are
#'simulated first. Then, more recent classes are subsequently added. Recent HBD
#'segments will mask more ancient HBD segments. See details in the method
#'published in Molecular Ecology (Druet and Gautier, 2017). Genotypes are
#'simulated using provided allele frequencies (from the sample), their genetic
#'distances, the number of chromosomes and the number of SNPs. These simulations
#'do not take into account linkage disequilibrium information. This simulation
#'tool can for instance be used to check whether there is enough information in
#'the data set to estimate HBD segments and their partitioning in multiple
#'classes. Note that this simulation tool is not computationally efficient.
#'
#'@param simdata The name of the zdata object created by the zoodata function.
#'
#'@param simmodel The name of the zmodel object created by the zoomodel
#'  function. The simulation program uses only the number of classes, their rate
#'  and their mixing coefficient (it is the same for a pre-defined model or
#'  not).
#'
#'@param nsim The number of simulated individuals.
#'
#'@param fullout Indicates whether a more detailed output is requested.
#'
#'@return The function simulates genotypes with the properties of the data set
#'  (number of SNPs, genetic map, allele frequencies) and according to the
#'  specified model (number of HBD classes, rate of the classes and proportion
#'  of mixing). The simulation is created as in Druet and Gautier (2017) and new
#'  autozygosity masks more ancient autozygosity. The output is a new zoodata
#'  object that can be analyzed with the zoorun function.
#'
#'  If a more detailed output is requested with the fullout parameter, then the
#'  function returns list (for instance simres) containing the zoodata, a matrix
#'  with realized inbreeding per individual, estimated at SNP positions, in the
#'  different classes (with 1 being the most ancient class (non-HBD) and the
#'  highest number corresponds to the most recent class), and a matrix
#'  containing for each individual and at each marker position the simulated
#'  class (1 for non-HBD, 2 for most ancient HBD, etc). These elements can be
#'  accessed using the simres[[1]], simres[[2]] and simres[[3]], respectively.
#'@export

zoosimd <- function(simdata, simmodel, nsim, fullout = FALSE) {

  nlayers = length(simmodel@krates)-1
  runsim <- .Fortran("zoosim",as.integer(nlayers),as.integer(simdata@nchr),as.integer(simdata@nsnps),
                 as.double(simdata@freqs),as.integer(simdata@bp),
                 as.integer(simdata@chrbound),as.double(simmodel@krates[nlayers:1]),
                 as.double(simmodel@mix_coef[nlayers:1]),as.double(simmodel@err),as.integer(nsim),
                 matrix(as.integer(0),nsim,simdata@nsnps),
                 matrix(as.double(0),nsim,nlayers),matrix(as.integer(0),nsim,simdata@nsnps))

  simout <- new("zdata")
  simout@genos <- t(runsim[[11]])
  simout@bp <- simdata@bp
  simout@chrnames <- simdata@chrnames
  simout@chrbound <- simdata@chrbound
  simout@nind <- nsim
  simout@nsnps <- simdata@nsnps
  freqs <- apply(simout@genos, 1, getfreq1)
  simout@freqs <- freqs
  simout@nchr <- simdata@nchr
  simout@zformat = "gt"

  if(!fullout)return(simout)
  if(fullout){
    simres=list(simout,runsim[[12]],t(runsim[[13]]))
    return(simres)
  }
}

