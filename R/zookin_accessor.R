#'Computes the realized kinship
#'
#'Computes the realized kinship with respect to a base population by summing the
#'relatedness from all classes with a rate lower or equal to the threshold T.
#'This amounts to set the base population approximately 0.5*T generations ago.
#'Classes with a higher rate are then no longer considered as IBD.
#'
#'@param kres The name of the kres object created by the zoorun function.
#'
#'@param T The value chosen to define the base population. When T is not
#'  provided, all classes are considered to estimate the kinship.
#'
#'@return An array with the computed kinship for all the pairs of individuals in
#'  the analysis.
#'
#'@export

cumkin <- function(kres, T = NULL) {
  K = length(kres@krates)
  npairs = nrow(kres@pairs)
  kin <- array(0,npairs)
  if(K > 2){x <- t(apply(kres@realized[,1:(K-1)],1,cumsum))}
  if(K == 2){x <- kres@realized[,1]}
  if(is.null(T)){
    if(K > 2){kin = x[,(K-1)]}
    if(K == 2){kin =x}
  }else{
    T2=array(0,npairs)
    for (id in 1:npairs){
      for (j in 1:(K-1)){
        rk = kres@krates[j]
        if(rk <= T){T2[id] = j}
      }
      if(T2[id] > 0){
        if(K > 2){kin[id] = x[id,T2[id]]}
        if(K == 2){kin[id] = x[id]}
        }
    }
  }

  return(kin)
}

#'Extracts the IBD probabilities from the kres object
#'
#'Extracts the locus specific IBD probabilities for a pair of individuals. This
#'is the probability that one pair of chromosomes (one sampled in each
#'individual) are IBD at that position. A specific chromosomal region can be
#'specified. A threshold T can be used to determine which classes are used in
#'the computation of the IBD probability. This function requires that the option
#'"localhbd" was set to TRUE when creating the kres object.
#'
#'@param kres The name of the kres object created by the zookin function.
#'
#'@param zooin The name of the zdata object created by the zoodata function. See
#'  "zoodata" for more details.
#'
#'@param num The number of the pair of individuals to extract (correspond to the
#'  position in the list of pairs).
#'
#'@param chrom the number of the chromosome where we are looking for IBD
#'  segments. This chromosome number refers to the position of the chromosome in
#'  the list of all chromsomes present in the input genotype data.
#'
#'@param startPos The starting position (on the chromosome) of the interval from
#'  which we will extract IBD segments (1 by default).
#'
#'@param endPos The ending position (on the chromosome) of the interval from
#'  which we will extract IBD segments (last position by default).
#'
#'@param T The value chosen to define the base population (to determine which
#'  classes are used in estimated of IBD probability, which classes are
#'  considered autozygous). When T is not provided, all IBD classes are
#'  considered to estimate the local IBD probability.

#'
#'@return The function returns a vector of IBD probabilities averaged over the
#'  four combination of pairs of parental haplotypes for the specified pair of
#'  individuals and chromosomal region. The IBD probabilities are computed as
#'  the sum of the probabilities for each class with a rate smaller or equal
#'  than the threshold (the sum from all the IBD classes when T is not
#'  specified).
#'
#'@export

predhbd <- function(kres, zooin, num, chrom = NULL, startPos = NULL, endPos = NULL, T = FALSE) {
ibdp <- array(0,zooin@nsnps)
  K=length(kres@krates)+1;K2=K-1
  if(T){
    K2=0
    for (k in 1:(K-1)){if(kres@krates[k] <= T){K2=k}}
  }
  if(is.null(chrom)){
    if( K2 < (K-1)){
      length(ibdp)
      length(kres@ibdp[[num]][1,])
      for (k in 1:K2){ibdp = ibdp + kres@ibdp[[num]][k,]}
    } else {
      ibdp = 1 - kres@ibdp[[num]][K,]
    }
  } else {
    if(is.null(startPos)){startPos = 0}
    if(is.null(endPos)){endPos = zooin@bp[zooin@chrbound[chrom,2]]}
     ff <- (zooin@bp >= startPos & zooin@bp <= endPos &
              seq(1,zooin@nsnps) >= zooin@chrbound[chrom,1] &
              seq(1,zooin@nsnps) <= zooin@chrbound[chrom,2])
     ibdp = ibdp[ff]
     if( K2 < (K-1)){
        for (k in 1:K2){ibdp = ibdp + kres@ibdp[[num]][k,ff]}
     } else {
       ibdp = 1- kres@ibdp[[num]][K,ff]
     }
  }
  return(ibdp)
}



