
#'RZooRoH: A package for estimating global and local individual autozygosity.
#'
#'Functions to identify Homozygous-by-Descent (HBD) segments associated with
#'runs of homozygosity (RoH) and to estimate individual autozygosity (or
#'inbreeding coefficient). HBD segments and autozygosity are assigned to
#'multiple HBD classes with a model-based approach relying on a mixture of
#'exponential distributions. The rate of the exponential distribution is
#'distinct for each HBD class and defines the expected length of the HBD
#'segments. These HBD classes are therefore related to the age of the segments
#'(longer segments and smaller rates for recent autozygosity / recent common
#'ancestor). The functions allow to estimate the parameters of the model (rates
#'of the exponential distributions, mixing proportions), to estimate global and
#'local autozygosity probabilities and to identify HBD segments with the Viterbi
#'decoding.
#'
#'@section Data pre-processing:  Note that the model is designed for autosomes.
#'  Other chromosomes and additional filtering (e.g. call rate, missing, HWE,
#'  etc.) should be performed prior to run RZooRoH with tools such as plink or
#'  bcftools. The model works on an ordered map and ignores SNPs with a null
#'  position.
#'
#'@section RZooRoH functions: The main functions included in the package are
#'  zoodata(), zoomodel() and zoorun(). There are also four function to plot
#'  the results: zooplot_partitioning(), zooplot_hbdseg(), zooplot_prophbd()
#'  and zooplot_individuals(). Finally, four accessors functions help to
#'  extract the results: realized(), cumhbd(), rohbd() and probhbd().
#'
#'  You can obtain individual help for each of the functions. By typing for
#'  instance: help(zoomodel) or ? zoomodel.
#'
#'  To run RZooRoH, you must first load your data with the zoodata() function.
#'  It will create a zooin object required for further analysis. Next, you need
#'   to define the model you want to run. You can define
#'  a default model by typing for instance, my.mod <- zoomodel(). Finally, you
#'  can run the model with the zoorun function. You can choose to estimate
#'  parameters with different procedures, estimate global and local
#'  homozygous-by-descent (HBD) probabilities with the Forward-Backward
#'  procedure or identify HBD segments with the Viterbi algorithm. The results
#'  are saved in a zres object.
#'
#'  The four plot functions zooplot_partitioning(), zooplot_hbdseg(),
#'  zooplot_prophbd() and zooplot_individuals() use zres objects to make different
#'  graphics. Similarly, the accessor functions help to extract information
#'  from the zres objects (see vignette for more details).
#'
#'  To get the list of data sets (for examples):
#'
#'  data(package="RZooRoH")
#'
#'  And to get the description of one data set, type ? name_data (with name_data
#'  being the name of the data set). For instance:
#'
#'  ? genosim
#'
#' @examples
#'
#' # Start with a small data set with six individuals and external frequencies.
#' freqfile <- (system.file("exdata","typsfrq.txt",package="RZooRoH"))
#' typfile <- (system.file("exdata","typs.txt",package="RZooRoH"))
#' frq <- read.table(freqfile,header=FALSE)
#' typdata <- zoodata(typfile,supcol=4,chrcol=1,poscol=2,allelefreq=frq$V1)
#' # Define a model with two HBD classes with rates equal to 10 and 100.
#' Mod3R <- zoomodel(K=3,base_rate=10)
#' # Run the model on all individuals.
#' typ.res <- zoorun(Mod3R,typdata)
#' # Observe some results: likelihood, realized autozygosity in different
#' # HBD classes and identified HBD segments.
#' typ.res@modlik
#' typ.res@realized
#' typ.res@hbdseg
#' # Define a model with one HBD and one non-HBD class and run it.
#' Mod1R <- zoomodel(K=2,predefined=FALSE)
#' typ2.res <- zoorun(Mod1R,typdata)
#' # Print the estimated rates and mixing coefficients.
#' typ2.res@krates
#' typ2.res@mixc
#'
#' # Get the name and location of a second example file.
#' myfile <- (system.file("exdata","genoex.txt",package="RZooRoH"))
#' # Load your data with default format:
#' example2 <- zoodata(myfile)
#' # Define the default model:
#' my.model <- zoomodel()
#' # Run RZooRoH on your data with the model (parameter estimation with optim). This can
#' # take a few minutes because it is a large model for 20 individuals:
#' \donttest{my.res <- zoorun(my.model,example2)}
#' # To estimate the parameters with the EM-algorithm, run the Forward-Backward
#' # algorithm to estimate realized autozygosity and the Viterbi algorithm to
#' # identify HBD segments (a few mintues too, see above).
#' \donttest{my.res2 <- zoorun(my.model, example2, fb=TRUE, vit=TRUE, method = "estem")}
#'# To run the model on a subset of individuals with 1 thread:
#' \donttest{my.res3 <- zoorun(my.model, example2, ids=c(7,12,16,18), nT = 1)}
#' # Define a smaller model and run it on two individuals.
#' my.mod2 <- zoomodel(K=4,base_rate=10)
#' \donttest{my.res4 <- zoorun(my.mod2, example2, ids=c(9,18))}
#'
#'@import graphics
#'@import stats
#'@import utils
#'@docType package
#'@name RZooRoH
NULL

setClass(Class = "zres",
         representation(nind = "numeric", ids = "vector", mixc = "matrix", krates = "matrix", realized = "matrix",
                        hbdp = "list", hbdseg ="data.frame", modlik ="vector", modbic ="vector", niter = "vector",
                        optimerr = "vector", sampleids = "vector")
)


is.zres <- function (x)
{
  res <- (is(x,"zres") & validObject(x))
  return(res)
}

setClass(Class = "oneres",
         representation(mixc = "vector", krates = "vector", hbdp = "matrix", hbdseg ="data.frame",
                        realized = "vector", modlik ="numeric", modbic ="numeric", num ="numeric",
                        niter = "numeric", optimerr ="numeric")
)

is.oneres <- function (x)
{
  res <- (is(x,"oneres") & validObject(x))
  return(res)
}

#### estimate parameters for one individual

#'Run the ZooRoH model
#'
#'Apply the defined model on a group of individuals: parameter estimation,
#'computation of realized autozygosity and homozygous-by-descent probabilities,
#' and identification of HBD segments (decoding).
#'
#'@param zoomodel A valid zmodel object as defined by the zoomodel function. The
#'  model indicates whether rates of exponential distributions are estimated or
#'  predefined, the number of classes, the starting values for mixing
#'  coefficients and rates, the error probabilities. See "zoomodel" for more
#'  details.
#'
#'@param zooin A valid zdata object as obtained by the zoodata function. See
#'"zoodata" for more details.
#'
#'@param ids An optional argument indicating the individual (its position in the
#'  data file) that must be proceeded. It can also be a vector containing the
#'  list of numbers that must be proceeded. By default, the model runs for all
#'  individuals.
#'
#'@param method Specifies whether the parameters are estimated by optimization
#'  with the L-BFGS-B method from the optim function (option method="opti",
#'  default value) or with the EM-algortihm (option method="estem"). When the EM
#'  algorithm is used, the Forward-Backward algorithm is run automatically (no
#'  need to select the fb option - see below). We recommend to set minr to 1
#'  when using the EM algorithm. If the user don't want to estimate the
#'  parameters he must set method="no".
#'
#'@param fb A logical indicating whether the forward-backward algorithm is run
#'  (optional argument - true by default). The Forward-Backward algorithm
#'  estimates the local probabilities to belong to each HBD or non-HBD class. By
#'  default, the function returns only the HBD probabilities for each class,
#'  averaged genome-wide, and corresponding to the realized autozygosity
#'  associated with each class. To obtain HBD probabilities at every marker
#'  position, the option localhbd must be set to true (this generates larger
#'  outputs).
#'
#'@param vit A logical indicating whether the Viterbi algorithm is run (optional
#'  argument - true by default). The Viterbi algorithm performs the decoding
#'  (determining the underlying class at every marker position). Whereas the
#'  Forward-Backward algorithms provide HBD probabilities (and how confident a
#'  region can be declared HBD), the Viterbi algorithm assigns every marker
#'  position to one of the defined classes (HBD or non-HBD). When informativity
#'  is high (many SNPs per HBD segments), results from the Forward-Backward and
#'  the Viterbi algorithm are very similar. The Viterbi algorithm is best suited
#'  to identify HBD segments. To estimate realized inbreeding and determine HBD
#'  status of a position, we recommend to use the Forward-Backward algorithm
#'  that better reflects uncertainty.
#'
#'@param localhbd A logical indicating whether the HBD probabilities for each
#'  individual at each marker are returned when using the EM or the
#'  Forward-Backward algorithm (fb option). This is an optional
#'  argument that is false by default.
#'
#'@param nT Indicates the number of threads used when running RZooRoH in
#'  parallel (optional argument - one thread by default).
#'
#'@param maxiter Indicates the maximum number of iterations with the EM
#'  algorithm (optional argument - 1000 by default).
#'
#'@param convem Indicates the convergence criteria for the EM algorithm
#'  (optional argument / 1e-10 by default).
#'
#'@param minr With optim and the reparametrized model (default), this indicates
#'  the minimum difference between rates of successive classes. With the EM
#'  algorithm, it indicates the minimum rate for a class. It is an optional
#'  argument set to 0. Adding such constraints might slow down the speed of
#'  convergence with optim and we recommend to run first optim without these
#'  constraints. With the EM algorithm, we recommend to use a value of 1.
#'
#'@param maxr With optim and the reparametrized model (default), this indicates
#'  the maximum difference between rates of successive classes. With the EM
#'  algorithm, it indicates the maximum rate for a class. It is an optional
#'  argument set to an arbitrarly large value (100000000). Adding such
#'  constraints might slow down the speed of convergence with optim and we
#'  recommend to run first optim without these constraints.
#'
#'@return The function return a zoores object with several slots accesses by
#' the "@" symbol. The three main results are zoores@@realized (the matrix with
#' partitioning of the genome in different HBD classes for each individual),
#' zoores@@hbdseg (a data frame with identified HBD segments) and zoores@@hbdp
#' (a list of matrices with HBD probabilities per SNP and per class).
#'
#' Here is a list with all the slots and their description:
#'
#' \enumerate{
#' \item zoores@@nind the number of individuals in the analysis,
#' \item zoores@@ids a vector containing the numbers of the analyzed individuals
#' (their position in the data file),
#' \item zoores@@mixc the (estimated) mixing coefficients per class for all
#' individuals,
#'  \item zoores@@krates the (estimated) rates for the exponential distributions
#'  associated with each HBD or non-HBD class for all individuals,
#'  \item zoores@@niter the number of iterations for estimating the
#'  parameters (per individual),
#'  \item zoores@@modlik a vector containing the likelihood of the model for
#'   each individual,
#'  \item zoores@@modbic a vector containing the value of the BIC for each
#'  individual,
#'  \item zoores@@realized a matrix with estimated realized autozygosity per HBD
#'  class (columns) for each individual (rows). These values are obtained with
#'  the Forward-Backward algorithm - fb option),
#'  \item zoores@@hbdp a list of matrices with the
#'  local probabilities to belong to an underlying hidden state (computed for
#'  every class and every individual). Each matrix has one row per snp and
#'  one column per class. To access the matrix from individual i, use the
#'  brackets "[[]]", for instance zoores@@hbdp[[i]],
#'  \item zoores@@hbdseg a data frame with the list
#'  of identified HBD segments with the Viterbi algorithm (the columns are the
#'  individual number, the chromosome number, the first and last SNP of the
#'  segment, the positions of the first and last SNP of the segment, the number
#'  of SNPs in the segment, the length of the segment, the HBD state of the
#'  segment),
#'  \item zoores@@optimerr a vector indicating whether optim ran with or
#'  without error (0/1),
#'  \item zoores@@sampleids is a vector with the names of the
#'  samples (when provided in the zooin object through the zoodata function).}
#'
#'
#'@examples
#'
#' # Start with a small data set with six individuals and external frequencies.
#' freqfile <- (system.file("exdata","typsfrq.txt",package="RZooRoH"))
#' typfile <- (system.file("exdata","typs.txt",package="RZooRoH"))
#' frq <- read.table(freqfile,header=FALSE)
#' typ <- zoodata(typfile,supcol=4,chrcol=1,poscol=2,allelefreq=frq$V1)
#' # Define a model with two HBD classes with rates equal to 10 and 100.
#' Mod3R <- zoomodel(K=3,base_rate=10)
#' # Run the model on all individuals.
#' typ.res <- zoorun(Mod3R, typ)
#' # Observe some results: likelihood, realized autozygosity in different
#' # HBD classes and identified HBD segments.
#' typ.res@modlik
#' typ.res@realized
#' typ.res@hbdseg
#' # Define a model with one HBD and one non-HBD class and run it.
#' Mod1R <- zoomodel(K=2,predefined=FALSE)
#' typ2.res <- zoorun(Mod1R, typ)
#' # Print the estimated rates and mixing coefficients.
#' typ2.res@krates
#' typ2.res@mixc

#' # Get the name and location of a second example file and load the data:
#' myfile <- (system.file("exdata","genoex.txt",package="RZooRoH"))
#' ex2 <- zoodata(myfile)
#' # Run RZooRoH to estimate parameters on your data with the 1 HBD and 1 non-HBD
#' # class (parameter estimation with optim).
#' my.mod1R <- zoomodel(predefined=FALSE,K=2,krates=c(10,10))
#' \donttest{my.res <- zoorun(my.mod1R, ex2, fb = FALSE, vit = FALSE)}
#' # The estimated rates and mixing coefficients:
#' \donttest{my.res@mixc}
#' \donttest{my.res@krates}
#' # Run the same model and run the Forward-Backward alogrithm to estimate
#' # realized autozygosity and the Viterbi algorithm to identify HBD segments:
#' \donttest{my.res2 <- zoorun(my.mod1R, ex2)}
#' # The table with estimated realized autozygosity:
#' \donttest{my.res2@realized}
#' # Run a model with 4 classes (3 HBD classes) and estimate the rates of HBD
#' # classes with one thread:
#' my.mod4R <- zoomodel(predefined=FALSE,K=4,krates=c(16,64,256,256))
#' \donttest{my.res3 <- zoorun(my.mod4R, ex2, fb = FALSE, vit = FALSE, nT =1)}
#' # The estimated rates for the 4 classes and the 20 individuals:
#' \donttest{my.res3@krates}
#' # Run a model with 5 classes (4 HBD classes) and predefined rates.
#' # The model is run only for a subset of four selected individuals.
#' # The parameters are estimated with the EM-algorithm, the Forward-Backward
#' # alogrithm is ued to estimate realized autozygosity and the Viterbi algorithm to
#' # identify HBD segments. One thread is used.
#' mix5R <- zoomodel(K=5,base=10)
#' \donttest{my.res4 <- zoorun(mix5R,ex2,ids=c(7,12,16,18), method = "estem", nT = 1)}
#' # The table with all identified HBD segments:
#' \donttest{my.res4@hbdseg}
#'
#'@useDynLib RZooRoH, .registration = TRUE
#'@export
#'@import iterators
#'@import foreach
#'@import parallel
#'@import doParallel
#'@import methods

zoorun <- function(zoomodel, zooin, ids = NULL, method = "opti",
                   fb = TRUE, vit = TRUE, minr = 0, maxr = 100000000,
                   maxiter = 1000, convem = 1e-10, localhbd = FALSE, nT = 1){

  if(is.null(ids)){ids=seq(from=1,to=zooin@nind)}
  if(is.zmodel(zoomodel) == FALSE){print("First element is not a valid model - see help for zoomodel function.")}
  if(is.zdata(zooin) == FALSE){print("First element is not a valid data set - see help for zoodata function.")}
  zrates <- zoomodel@krates
  zmix <- zoomodel@mix_coef
  gerr <- zoomodel@err
  seqerr <- zoomodel@seqerr

  if(.Platform$OS.type == "windows" & nT > 1){
    warning("Multithreading is not supported under windows. Only one thread will be used.\n")}
  if(.Platform$OS.type != "windows"){registerDoParallel(cores= nT)}

  if(zoomodel@typeModel != "mixkr" & zoomodel@typeModel != "MIXKR"
     & zoomodel@typeModel != "kr" & zoomodel@typeModel != "KR"){
    print("The model type must be kr or mixkr")
  }

  opti=FALSE;estem=FALSE
  if(method == "opti"){opti=TRUE}
  if(method == "estem"){estem=TRUE}
  if(method != "opti" & method != "estem" & method !="no"){
    print(c("Unrecognized method ::",method))
    print("Will use optim.")
    opti = TRUE
  }

  if(estem & minr < 1){
    print("We recommend to set minr to 1 with the EM algorithm")
  }

  if(estem){fb = FALSE}

  zoores <- new("zres")
  zoores@nind = length(ids)
  zoores@ids = ids
  zoores@mixc = matrix(0,zoores@nind,length(zmix))
  zoores@krates = matrix(0,zoores@nind,length(zrates))

  if(fb || estem){zoores@realized = matrix(0,zoores@nind,length(zmix))}

  i=NULL
  if(estem & minr < 1) minr = 1
  if(.Platform$OS.type == 'windows'){
    if(zoomodel@typeModel == "mixkr" || zoomodel@typeModel == "MIXKR"){
      xx <- foreach (i=1:length(ids), .combine=c) %do% {
        id <- ids[i]
        run_mixkr(zooin,id, zrates, zmix, opti, estem, fb, vit, gerr, seqerr, minr, maxr, maxiter, convem, localhbd)
      }
    }
    if(zoomodel@typeModel == "kr" || zoomodel@typeModel == "KR"){
      xx <- foreach (i=1:length(ids), .combine=c) %do% {
        id <- ids[i]
        run_kr(zooin,id, zrates, zmix, opti, estem, fb, vit, gerr, seqerr, minr, maxr, maxiter, convem, localhbd)
      }
    }
  }else{
    if(zoomodel@typeModel == "mixkr" || zoomodel@typeModel == "MIXKR"){
      xx <- foreach (i=1:length(ids), .combine=c) %dopar% {
        id <- ids[i]
        run_mixkr(zooin,id, zrates, zmix, opti, estem, fb, vit, gerr, seqerr, minr, maxr, maxiter, convem, localhbd)
      }
    }
    if(zoomodel@typeModel == "kr" || zoomodel@typeModel == "KR"){
      xx <- foreach (i=1:length(ids), .combine=c) %dopar% {
        id <- ids[i]
        run_kr(zooin,id, zrates, zmix, opti, estem, fb, vit, gerr, seqerr, minr, maxr, maxiter, convem, localhbd)
      }
    }
  }

  #### transfer list of oneres into zoores object

  if(length(ids) == 1){
    xx <- array(list(xx),1)
  }

  zoores@ids <- unlist(lapply(xx,slot,"num"))
  zoores@modlik <- unlist(lapply(xx,slot,"modlik"))
  zoores@modbic <- unlist(lapply(xx,slot,"modbic"))
  zoores@niter <- unlist(lapply(xx,slot,"niter"))
  zoores@mixc <- do.call("rbind",lapply(xx,slot,"mixc"))
  zoores@krates <- do.call("rbind",lapply(xx,slot,"krates"))
  zoores@optimerr <- unlist(lapply(xx,slot,"optimerr"))
  if(fb || estem){
    zoores@realized <- do.call("rbind",lapply(xx,slot,"realized"))
    if(localhbd){zoores@hbdp <- lapply(xx,slot,"hbdp")}
  }
  if(vit){zoores@hbdseg <- do.call("rbind", lapply(xx,slot,"hbdseg"))}
#  zoores@sampleids <- .GlobalEnv$zooin@sample_ids[zoores@ids]
  zoores@sampleids <- zooin@sample_ids[zoores@ids]
  return(zoores)
}

#### run mixkr model

run_mixkr <- function(zooin,id, zrates, zmix, opti = TRUE, estem = FALSE, fb = FALSE, vit = FALSE,
                      gerr = 0.001, seqerr = 0.001, minr = 0, maxr = 100000000,
                      maxiter = 1000, convem = 1e-10, localhbd = FALSE){
  K <- length(zmix)
#  zrs=NULL;pem=NULL
#  zrs <<- zrates
#  pem <<- pemission(id, gerr, seqerr, zformat = .GlobalEnv$zooin@zformat)
  pem <- pemission(zooin, id, gerr, seqerr, zformat = zooin@zformat)

  if(opti){
    zpar <- array(0,K-1)
    for (j in 1:(K-1)){zpar[j]=log(zmix[j]/zmix[K])}

    niter=NULL;niter <<- 0
    checkoptim=NULL;checkoptim <<- 1
    tryCatch(optires <- optim(zpar, lik_mixkr, zooin = zooin, pem = pem, zrs = zrates, method="L-BFGS-B",
                     control = list(trace=TRUE, REPORT=10000)),
             error=function(e){cat("Warning, error returned by optim - replaced by EM ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
    if(.GlobalEnv$checkoptim ==0){estem = TRUE; fb =FALSE}
    if(.GlobalEnv$checkoptim == 1){
      zmix <- 1/(1+exp(-optires$par[1:(K-1)]))
      zmix[K] <- 1 - sum(zmix[1:(K-1)])
    }
  }
  if(estem){
    estimrates <- 0
    onerate <- 0
    outem <- EM(zooin,id, pem, zrates, zmix, estimrates, onerate, maxiter, minr, maxr, convem, localhbd)
    print(paste("EM algorithm, iterations and loglik ::",outem[[5]],outem[[3]],sep=" "))
    zmix <- outem[[2]]
  }
  if(fb){
    outfb <- fb(zooin,id,pem,zrates,zmix,localhbd)
  }
  if(vit){
    outvit <- viterbi(zooin,id,pem,zrates,zmix)
  }

  ##### renvoie toujours les paramètres
  zresu <- new("oneres")
  loglik <- NA
  bicv <- NA
  zrates <- round(zrates,3)
  if(opti){
    if(.GlobalEnv$checkoptim == 1){
    loglik <- -optires$value ### optim miminizes -log(lik)
#    bicv <- -2 * loglik + log(.GlobalEnv$zooin@nsnps)*(K-1)}
    bicv <- -2 * loglik + log(zooin@nsnps)*(K-1)}
  }
  if(estem){
    loglik <- outem[[3]]
#    bicv <- -2 * loglik + log(.GlobalEnv$zooin@nsnps)*(K-1)
    bicv <- -2 * loglik + log(zooin@nsnps)*(K-1)
    zrates <- outem[[1]]
    zmix <- outem[[2]]
    .GlobalEnv$niter <- outem[[5]]
    zresu@realized <- outem[[4]]
    if(localhbd){zresu@hbdp <- outem[[6]]}
  }

  zresu@mixc <- zmix
  zresu@krates <- zrates
  zresu@modlik <- loglik
  zresu@modbic <- bicv
  zresu@num <- id
  zresu@niter <- .GlobalEnv$niter
  zresu@optimerr <- -1
  if(opti){zresu@optimerr <- .GlobalEnv$checkoptim}

  if(fb){
    if(!localhbd){zresu@realized <- outfb}
    if(localhbd){
      zresu@realized <- outfb[[1]]
      zresu@hbdp <- outfb[[2]]
      }
  }
  if(vit){
    zresu@hbdseg <- outvit
  }
  return(zresu)
}

#### run kr model

run_kr <- function(zooin,id, zrates, zmix, opti = TRUE, estem = FALSE, fb = FALSE, vit = FALSE,
                   gerr = 0.001, seqerr = 0.001, minr = 1, maxr = 100000000,
                   maxiter = 1000, convem = 1e-10, localhbd = FALSE){
#  pem=NULL
#  pem <<- pemission(id,  gerr, seqerr, zformat = .GlobalEnv$zooin@zformat)
  pem <- pemission(zooin,id,  gerr, seqerr, zformat = zooin@zformat)
  K <- length(zmix) #### vérifier que length(zmix) = length (zrates)
  if(opti){
    if(K > 2){
      zpar <- array(0,2*K-1)
      zpar[1] <- log(zrates[1])
      for(j in 2:(K-1)){zpar[j] <- log(zrates[j]-zrates[j-1])}
      zpar[K]=log(zrates[K])
      for (j in 1:(K-1)){zpar[K+j]=log(zmix[j]/zmix[K])}
    } else { #### 1R model
      zpar <- array(0,2)
      zpar[1] <- log(zrates[1]) ### the common rate
      zpar[2] <- log(zmix[1]/zmix[2])
    }

    niter=NULL;niter <<- 0
    checkoptim=NULL;checkoptim <<- 1

    if(maxr < 100000000 | minr > 0){
    if(K == 2){
      lowpar <- c(log(minr), -Inf)
      uppar <- c(log(maxr), Inf)
    }
    if(K > 2){
      lowpar <- c(rep(log(minr),K), rep(-Inf,(K-1)))
      uppar <- c(rep(log(maxr),K), rep(Inf,(K-1)))
    }
    tryCatch(optires <- optim(zpar, lik_kr, zooin = zooin, pem = pem, method="L-BFGS-B", lower = lowpar, upper = uppar,
                              control = list(trace=TRUE, REPORT=10000)),
             error=function(e){cat("Warning, error returned by optim - replaced by EM ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
    }else{
      tryCatch(optires <- optim(zpar, lik_kr, zooin = zooin, pem = pem, method="L-BFGS-B",
                                control = list(trace=TRUE, REPORT=10000)),
               error=function(e){cat("Warning, error returned by optim - replaced by EM ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
    }

    if(.GlobalEnv$checkoptim ==0){estem = TRUE; fb =FALSE}

    if(.GlobalEnv$checkoptim == 1){
      if(K > 2){
        zrates[1] <- exp(optires$par[1])
        for (j in 2:(K-1)){zrates[j] <- zrates[j-1] + exp(optires$par[j])}
        zrates[K] <- exp(optires$par[K])
        zmix <- 1/(1+exp(-optires$par[(K+1):(2*K-1)]))
        zmix[K] <- 1 - sum(zmix[1:(K-1)])
      } else {
        zrates <- exp(optires$par[1])
        zrates[2] <- zrates[1]
        zmix[1] <- 1/(1+exp(-optires$par[2]))
        zmix[2] <- 1 - zmix[1]
        }
    }
  }

  if(estem){
    estimrates <- 1
    onerate <- 0
    if (K ==2) onerate <- 1
    if(minr < 1){minr=1}
    outem <- EM(zooin,id, pem, zrates, zmix, estimrates, onerate, maxiter, minr, maxr, convem, localhbd)
    print(paste("EM algorithm, iterations and loglik ::",outem[[5]],outem[[3]],sep=" "))
    zrates <- outem[[1]]
    zmix <- outem[[2]]
  }
    if(fb){
    outfb <- fb(zooin,id,pem,zrates,zmix, localhbd)
  }
  if(vit){
    outvit <- viterbi(zooin,id,pem,zrates,zmix)
  }

  ##### renvoie toujours les paramètres
  zresu <- new("oneres")
  loglik <- NA
  bicv <- NA
  zrates <- round(zrates,3)
  if(opti){
    if (.GlobalEnv$checkoptim == 1){
    loglik <- -optires$value ### optim miminizes -log(lik)
#    if(K > 2)bicv <- -2 * loglik + log(.GlobalEnv$zooin@nsnps)*(2*K-1)
#    if(K == 2)bicv <- -2 * loglik + 2*log(.GlobalEnv$zooin@nsnps)
    if(K > 2)bicv <- -2 * loglik + log(zooin@nsnps)*(2*K-1)
    if(K == 2)bicv <- -2 * loglik + 2*log(zooin@nsnps)
    }}
  if(estem){
    loglik <- outem[[3]]
#    bicv <- -2 * loglik + log(.GlobalEnv$zooin@nsnps)*(K-1)
#    if(K > 2)bicv <- -2 * loglik + log(.GlobalEnv$zooin@nsnps)*(2*K-1)
#    if(K == 2)bicv <- -2 * loglik + 2*log(.GlobalEnv$zooin@nsnps)
    bicv <- -2 * loglik + log(zooin@nsnps)*(K-1)
    if(K > 2)bicv <- -2 * loglik + log(zooin@nsnps)*(2*K-1)
    if(K == 2)bicv <- -2 * loglik + 2*log(zooin@nsnps)
    zrates <- outem[[1]]
    zmix <- outem[[2]]
    .GlobalEnv$niter <- outem[[5]]
    zresu@realized <- outem[[4]]
    if(localhbd){zresu@hbdp <- outem[[6]]}
  }

  zresu@mixc <- zmix
  zresu@krates <- zrates
  zresu@modlik <- loglik
  zresu@modbic <- bicv
  zresu@num <- id
  zresu@niter <- .GlobalEnv$niter
  zresu@optimerr <- -1
  if(opti){zresu@optimerr <- .GlobalEnv$checkoptim}

  if(fb){
    if(!localhbd){zresu@realized <- outfb}
    if(localhbd){
      zresu@realized <- outfb[[1]]
      zresu@hbdp <- outfb[[2]]
    }
  }
  if(vit){
    zresu@hbdseg <- outvit
  }

  return(zresu)
}


