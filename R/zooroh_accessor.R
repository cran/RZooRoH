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
  K = length(zres@krates[1,])+1
  if(is.null(classNum)){
    classNum = seq(1,K)
  }
  matout <- as.data.frame(zres@realized[,classNum])
  K2 <- length(classNum)
  classNames <- array('',K2)
  for (i in 1:K2){
    if(sd(zres@krates[,K-1])!=0){
      if(classNum[i] < K){classNames[i]=paste("HBDclass",classNum[i],sep="")}
      if(classNum[i]==K){classNames[i]="NonHBD"}
    }
    if(sd(zres@krates[,K-1])==0){
      if(classNum[i] < K){classNames[i]=paste("R_",zres@krates[1,classNum[i]],sep="")}
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
  K = length(zres@krates[1,])+1
  if(sd(zres@krates[,K-1])!=0){
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
  K=length(zres@krates[id,])+1;K2=K-1
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

#'Merge several zres objects generated by zoorun
#'
#'The function is used for example when the analysis has been split in several
#'tasks. All the zres must come from the same original
#'zoodata (zoorun is applied to the same zoodata with the same model). In
#'addition, the data should not overlap (the same individual should not be
#'present in multiple zres).
#'
#'@param list_zres A list with the name of the zres objects to be merged.
#'
#'@return a single zres object containing the results from the merged zres
#'  objects. Note that the HBD segments are not sorted by ID number.
#'
#'@export

merge_zres <- function(list_zres){

  allres <- new("zres")

  allres@nind <- list_zres[[1]]@nind+list_zres[[2]]@nind
  allres@ids <- c(list_zres[[1]]@ids,list_zres[[2]]@ids)
  allres@mixc <- rbind(list_zres[[1]]@mixc,list_zres[[2]]@mixc)
  allres@krates <- rbind(list_zres[[1]]@krates,list_zres[[2]]@krates)
  allres@realized <- rbind(list_zres[[1]]@realized,list_zres[[2]]@realized)
  allres@modlik <- c(list_zres[[1]]@modlik,list_zres[[2]]@modlik)
  allres@modbic <- c(list_zres[[1]]@modbic,list_zres[[2]]@modbic)
  allres@niter <- c(list_zres[[1]]@niter,list_zres[[2]]@niter)
  allres@optimerr <- c(list_zres[[1]]@optimerr,list_zres[[2]]@optimerr)
  allres@sampleids <- c(list_zres[[1]]@sampleids,list_zres[[2]]@sampleids)
  allres@hbdp <- c(list_zres[[1]]@hbdp,list_zres[[2]]@hbdp) ### stays a list
  allres@hbdseg <- rbind(list_zres[[1]]@hbdseg,list_zres[[2]]@hbdseg)

  if(length(list_zres)>2){
    for (i in 3:length(list_zres)){
      allres@nind <- allres@nind + list_zres[[i]]@nind
      allres@ids <- c(allres@ids,list_zres[[i]]@ids)
      allres@mixc <- rbind(allres@mixc,list_zres[[i]]@mixc)
      allres@krates <- rbind(allres@krates,list_zres[[i]]@krates)
      allres@realized <- rbind(allres@realized,list_zres[[i]]@realized)
      allres@modlik <- c(allres@modlik,list_zres[[i]]@modlik)
      allres@modbic <- c(allres@modbic,list_zres[[i]]@modbic)
      allres@niter <- c(allres@niter,list_zres[[i]]@niter)
      allres@optimerr <- c(allres@optimerr,list_zres[[i]]@optimerr)
      allres@sampleids <- c(allres@sampleids,list_zres[[i]]@sampleids)
      allres@hbdp <- c(allres@hbdp,list_zres[[i]]@hbdp) ### stays a list
      allres@hbdseg <- rbind(allres@hbdseg,list_zres[[i]]@hbdseg)
    }
  }

  x <- order(allres@ids)

  allres@ids <- allres@ids[x]
  allres@mixc <- allres@mixc[x,]
  allres@krates <- allres@krates[x,]
  allres@realized <- allres@realized[x,]
  allres@modlik <- allres@modlik[x]
  allres@modbic <- allres@modbic[x]
  allres@niter <- allres@niter[x]
  allres@optimerr <- allres@optimerr[x]
  allres@sampleids <- allres@sampleids[x]
  allres@hbdp <- allres@hbdp[x] ### stays a list
  ## allres@hbdseg

  n1 <- length(allres@ids)
  n2 <- length(unique(allres@ids))

  if(n2 < n1){
    print("Some duplicate IDs were present in the input files.")
    print("The merged zres files should contain distinct individuals")
  }

  return(allres)

}

#'Update one main zres object with new results
#'
#'The function is used for example when the main analysis failed for one individual.
#'The analysis is repeated for that individual with other parameters. The function
#'can then be used to insert the new results in the main zres object. To avoid
#'generating to large files, the updated zres object can take the same name as zres1
#'doing zres1 <- update_zres(zres1,zres2).
#'
#'@param zres1 The main zres object that will be modified.
#'
#'@param zres2 The new zres object, with the new results that will be inserted in the
#'main zres object.
#'
#'@return an updated zres object containing the results from zres1 updated by those from
#'zres2.
#'
#'@export

update_zres <- function(zres1,zres2){

  x <- match(zres2@ids,zres1@ids)

  zres1@mixc[x,] <- zres2@mixc
  zres1@krates[x,] <- zres2@krates
  zres1@realized[x,] <- zres2@realized
  zres1@modlik[x] <- zres2@modlik
  zres1@modbic[x] <- zres2@modbic
  zres1@niter[x] <- zres2@niter
  zres1@optimerr[x] <- zres2@optimerr
  zres1@hbdp[x] <- zres2@hbdp

  zres1@hbdseg <- zres1@hbdseg[!zres1@hbdseg$id %in% zres2@ids,]
  zres1@hbdseg <- rbind(zres1@hbdseg,zres2@hbdseg)

  return(zres1)

}

