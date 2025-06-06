### FUNCTIONS FOR KL models

### Estimate the likelihood for a set of parameters

lik_kl <- function(zpar,nchr,nsnps,chrbound,bp,pem) {

  if(length(zpar) == 2){
    K <- 2 #### 1R model with two parameters
  } else {
    K <- length(zpar)/2 + 1 ## KL model with 2*(K-1) parameters (we don't need the non-HBD rate)
  }

  zrs=array(0,K-1)
  if(K > 2){
    zrs[1] <- 1 + exp(zpar[1])
    for (j in 2:(K-1)){zrs[j] <- zrs[j-1]+exp(zpar[j])}
  } else { ## 1 R model
    zrs[1] <- 1 + exp(zpar[1])
  }

  fs=array(0,(K-1))
  if(K > 2){
    for (j in 1:(K-1)){fs[j]=exp(zpar[K+j-1])/(1+exp(zpar[K+j-1]))}
  } else {
    D=1+exp(zpar[2])
    fs[1]=exp(zpar[2])/D
  }
  .GlobalEnv$niter <- .GlobalEnv$niter + 1
  oneiter <- .Fortran("zoolayerlik",as.integer(K), as.integer(nchr),as.integer(nsnps),
                      as.double(pem), as.integer(chrbound),as.double(zrs),as.double(fs),
                      as.integer(bp),as.double(0))
  loglik <- oneiter[9][[1]][1]
  return(-loglik)
}

