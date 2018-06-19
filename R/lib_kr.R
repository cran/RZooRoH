### FUNCTIONS FOR KR models

### Estimate the likelihood for a set of parameters

lik_kr <- function(zpar) {

  if(length(zpar) == 2){
    K <- 2 #### 1R model with two parameters
  } else {
    K <- (length(zpar)+1)/2 ## KR model with 2*K-1 parameters
  }

  zrs=array(0,K)
  if(K > 2){
    zrs[1] <- exp(zpar[1])
    for (j in 2:(K-1)){zrs[j] <- zrs[j-1]+exp(zpar[j])}
    zrs[K] <- exp(zpar[K])
  } else { ## 1 R model
    zrs[1] <- exp(zpar[1])
    zrs[2] <- zrs[1]
  }

  fs=array(0,K)
  if(K > 2){
    D=1
    for (j in 1:(K-1)){D=D+exp(zpar[K+j])}
    for (j in 1:(K-1)){fs[j]=exp(zpar[K+j])/D}
    fs[K]=1-sum(fs[1:(K-1)])
  } else {
    D=1+exp(zpar[2])
    fs[1]=exp(zpar[2])/D
    fs[2]=1-fs[1]
  }
  #  niter <<- niter + 1
  .GlobalEnv$niter <- .GlobalEnv$niter + 1
  oneiter <- .Fortran("zoolik",as.integer(K), as.integer(.GlobalEnv$zooin@nchr),as.integer(.GlobalEnv$zooin@nsnps),
                      as.double(.GlobalEnv$pem), as.integer(.GlobalEnv$zooin@chrbound),as.double(zrs),as.double(fs),
                      as.integer(.GlobalEnv$zooin@bp),as.double(0))
  loglik <- oneiter[9][[1]][1]
  return(-loglik)
}

