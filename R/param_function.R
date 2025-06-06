#### run mixkl model with step function

run_fun_mixkl <- function(zooin,id, zrates, zmix, XM, opti = TRUE, fb = FALSE, vit = FALSE,
                      gerr = 0.001, seqerr = 0.001, hemiprob = 0.000, maxr = 100000000, minmix = 1e-16,
                      maxiter = 1000, localhbd = FALSE, ibd = FALSE, ibdpair = NULL, haploid = FALSE,
                      optim_method = "L-BFGS-B",trim_ad = FALSE){
  K <- length(zrates)+1 ### number of states
  K2 <- length(zmix) ### number of mixing coefficients (without non-HBD)
  pem <- pemission(zooin, id, gerr, seqerr, hemiprob, zformat = zooin@zformat, ibd, ibdpair, haploid)
  skiperr=TRUE

  if(opti){
    zpar <- array(0,K2) ### number of parameters
    for (j in 1:K2){zpar[j]=log(zmix[j]/(1-zmix[j]))}
    if(minmix>0)myminpar=log(minmix/(1-minmix))
    niter=NULL;niter <<- 0
    checkoptim=NULL;checkoptim <<- 1

    if(zooin@zformat=="ad" & trim_ad){
      tad <- array(0,nrow(pem))
      for (i in 1:zooin@nchr){tad[zooin@chrbound[i,1]:zooin@chrbound[i,2]]=i}
      tad[(zooin@genos[,(2*id-1)]+zooin@genos[,(2*id)]) <= 1]=0
      pem2=pem[tad > 0,]
      nsnps2 = nrow(pem2)
      print(nrow(pem2))

      chrnum <- tad[tad > 0]
      snpnum <- seq(1:nsnps2)
      nchr2 <- zooin@nchr
      chrbound2 <- matrix(0, nchr2, 2)
      for (i in 1:nchr2){
        chrbound2[i,1]=min(snpnum[chrnum==i])
        chrbound2[i,2]=max(snpnum[chrnum==i])
      }
      bp2 <- zooin@bp[tad > 0]

      if(optim_method=="L-BFGS-B" & minmix > 0){
        tryCatch(optires <- optim(zpar, lik_fun_mixkl, nchr = nchr2 ,nsnps = nsnps2, chrbound = chrbound2,
                                  bp = bp2, pem = pem2, zrs = zrates, XM = XM, method=optim_method,
                                  lower=rep(myminpar,(K-1)),upper=rep(-myminpar,(K-1)),
                                  control = list(trace=TRUE, REPORT=10000,maxit=maxiter)),
                 error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
      }else{
        tryCatch(optires <- optim(zpar, lik_fun_mixkl, nchr = nchr2 ,nsnps = nsnps2, chrbound = chrbound2,
                                  bp = bp2, pem = pem2, zrs = zrates, XM = XM, method=optim_method,
                                  control = list(trace=TRUE, REPORT=10000,maxit=maxiter)),
                 error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
      }
    }else{

    if(optim_method=="L-BFGS-B" & minmix > 0){
    tryCatch(optires <- optim(zpar, lik_fun_mixkl, nchr = zooin@nchr, nsnps = zooin@nsnps, chrbound = zooin@chrbound,
                              bp = zooin@bp, pem = pem, zrs = zrates, XM = XM, method=optim_method,
                              control = list(trace=TRUE, REPORT=10000, maxit = maxiter),
                              lower=rep(myminpar,(K2)),upper=rep(-myminpar,(K2))),
             error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
    }else{
    tryCatch(optires <- optim(zpar, lik_fun_mixkl, nchr = zooin@nchr, nsnps = zooin@nsnps, chrbound = zooin@chrbound,
                              bp = zooin@bp, pem = pem, zrs = zrates, XM = XM, method=optim_method,
                              control = list(trace=TRUE, REPORT=10000, maxit = maxiter)),
              error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
    }}
    if(.GlobalEnv$checkoptim ==0){skiperr =FALSE} ### NEW
    if(.GlobalEnv$checkoptim == 1){
      zmix <- exp(optires$par[1:K2])/(1+exp(optires$par[1:K2]))
      zmix2 <- XM%*%zmix
    }
  }
  if(fb & skiperr){
    outfb <- fb_lay(zooin,id,pem,zrates,zmix2,localhbd)
  }
  if(vit & skiperr){
    outvit <- viterbi_lay(zooin,id,pem,zrates,zmix2)
  }

  ##### renvoie toujours les paramètres
  zresu <- new("oneres")
  loglik <- 0
  bicv <- 0
  zrates <- round(zrates,3)
  if(opti){
    if(.GlobalEnv$checkoptim == 1){
      loglik <- -optires$value ### optim miminizes -log(lik)
      bicv <- -2 * loglik + log(zooin@nsnps)*(K-1)}
  }

  zresu@mixc <- zmix
  zresu@krates <- zrates
  zresu@modlik <- loglik
  zresu@modbic <- bicv
  zresu@num <- id
  zresu@niter <- .GlobalEnv$niter
  zresu@optimerr <- -1
  if(opti & .GlobalEnv$checkoptim==0){zresu@optimerr <- 99}
  if(opti & .GlobalEnv$checkoptim==1){zresu@optimerr <- optires$convergence}

  if(!skiperr){zresu@mixc <- rep(0,K)}

  if(fb & skiperr){
    if(!localhbd){zresu@realized <- outfb}
    if(localhbd){
      zresu@realized <- outfb[[1]]
      zresu@hbdp <- outfb[[2]]
    }
  }
  if(fb & !skiperr){
    outfb1 <- matrix(as.double(0),K,zooin@nsnps)
    outfb0 <- apply(outfb1,1,mean)
    if(!localhbd){zresu@realized <- outfb0}
    if(localhbd){
      zresu@realized <- outfb0
      zresu@hbdp <- outfb1
    }
  }
  if(vit & skiperr){
    zresu@hbdseg <- outvit
  }
  return(zresu)
}

### Estimate the likelihood for a function of a set of parameters

lik_fun_mixkl <- function(zpar,nchr,nsnps,chrbound,bp,pem,zrs,XM) {
  zpar2=XM%*%zpar
  K <- length(zpar2)+1
  fs=array(0,(K-1))
  for (j in 1:(K-1)){fs[j]=exp(zpar2[j])/(1+exp(zpar2[j]))}
  .GlobalEnv$niter <- .GlobalEnv$niter + 1
  oneiter <- .Fortran("zoolayerlik",as.integer(K), as.integer(nchr),as.integer(nsnps),
                      as.double(pem),
                      as.integer(chrbound),as.double(zrs),as.double(fs),
                      as.integer(bp),as.double(0))
  loglik <- oneiter[9][[1]][1]
  return(-loglik)
}


