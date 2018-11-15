
### FUNCTIONS FOR MIXKR models

### Estimate the likelihood for a set of parameters

lik_mixkr <- function(zpar,zooin,pem, zrs) {
  K <- length(zpar)+1
  fs=array(0,K)
  D=1
  for (j in 1:(K-1)){D=D+exp(zpar[j])}
  for (j in 1:(K-1)){fs[j]=exp(zpar[j])/D}
  fs[K]=1-sum(fs[1:(K-1)])
  .GlobalEnv$niter <- .GlobalEnv$niter + 1
#  oneiter <- .Fortran("zoolik",as.integer(K), as.integer(.GlobalEnv$zooin@nchr),as.integer(.GlobalEnv$zooin@nsnps),
#                      as.double(.GlobalEnv$pem),
#                      as.integer(.GlobalEnv$zooin@chrbound),as.double(.GlobalEnv$zrs),as.double(fs),
#                      as.integer(.GlobalEnv$zooin@bp),as.double(0))
  oneiter <- .Fortran("zoolik",as.integer(K), as.integer(zooin@nchr),as.integer(zooin@nsnps),
                      as.double(pem),
                      as.integer(zooin@chrbound),as.double(zrs),as.double(fs),
                      as.integer(zooin@bp),as.double(0))
  loglik <- oneiter[9][[1]][1]
  return(-loglik)
}

### Run Forward-Backward algorithm for a set of parameters (identical for KR and MIXKR)

fb <- function(zooin,id,pem,zrates,zmix,localhbd = FALSE) {
  K <- length(zmix)
  onefb <- .Fortran("zooFB",as.integer(K), as.integer(zooin@nchr),as.integer(zooin@nsnps),
                    as.double(pem),
                    as.integer(zooin@chrbound),as.double(zrates),as.double(zmix),
                    as.integer(zooin@bp),as.double(0),matrix(as.double(0),K,zooin@nsnps))
  realized <- apply(onefb[[10]],1,mean)
  if(!localhbd){return(realized)}
  if(localhbd){return(list(realized,onefb[[10]]))}
}

### Run Viterbi algorithm for a set of parameters (identical for KR and MIXKR)
viterbi <- function(zooin,id,pem,zrates,zmix) {
  K <- length(zmix)
  decout <- .Fortran("zooViterbi",as.integer(K), as.integer(zooin@nchr),as.integer(zooin@nsnps),
                     as.double(pem), as.integer(zooin@chrbound),as.double(zrates),as.double(zmix),
                     as.integer(zooin@bp),array(as.integer(0),zooin@nsnps))

  states<-decout[[9]]
  hbdsegments <- data.frame(id = numeric(), chrom = numeric(), hbdst = numeric(), hbdend = numeric(), bpst = numeric(),
                            bpend = numeric(), len1 = numeric(), len2 = numeric(), hbdstate = numeric())
  for (chr in 1:zooin@nchr){
    chrstates <- states[zooin@chrbound[chr,1]:zooin@chrbound[chr,2]]
    posbp <- zooin@bp[zooin@chrbound[chr,1]:zooin@chrbound[chr,2]]
    ns <- length(chrstates)
    v2 <- c(chrstates[2:ns],0)
    v0 <- c(0,chrstates[1:(ns-1)])
    seg1 <- which(chrstates != v0)
    seg2 <- which(chrstates != v2)
    nseg <- length(seg1)
    lensnp <- (seg2 - seg1) + 1
    lenbp <- posbp[seg2] - posbp[seg1] + 1
    allseg <-  cbind(rep(id,nseg),rep(chr,nseg),seg1,seg2,posbp[seg1],posbp[seg2],lensnp,lenbp,chrstates[seg1])
    allseg <- as.data.frame(allseg)
    hbdsegments <- rbind(hbdsegments, allseg[allseg[,9] < K,]) #### assumes that non-hbd is last (use isf else)
  }

  colnames(hbdsegments) <- c('id','chrom','start_snp','end_snp','start_pos','end_pos','number_snp','length','HBDclass')
  return(hbdsegments)
}

### Run EM algorithm
EM <- function(zooin,id, pem, zrates, zmix, estimr = 1, oner = 0, nrounds = 1000, minr = 1, maxr = 1024, convem = 1e-10, localhbd = FALSE) {
  K <- length(zmix)
  ems <- .Fortran("zooEM",as.integer(K), as.integer(zooin@nchr),as.integer(zooin@nsnps),
                  as.double(pem), as.integer(zooin@chrbound),as.double(zrates),as.double(zmix),
                  as.integer(zooin@bp),as.double(0),matrix(as.double(0),K,zooin@nsnps),
                  as.integer(estimr),as.integer(oner),as.integer(nrounds),as.double(minr),
                  as.double(maxr),as.double(convem),as.integer(0))
  realized <- apply(ems[[10]],1,mean)
  if(!localhbd){return(list(ems[[6]],ems[[7]],ems[[9]],realized,ems[[17]]))}
  if(localhbd){return(list(ems[[6]],ems[[7]],ems[[9]],realized,ems[[17]],ems[[10]]))}
}
