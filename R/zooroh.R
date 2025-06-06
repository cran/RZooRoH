
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
#'decoding. Functions also allow to compute kinship between individuals and to
#'predict inbreeding in the future progeny of a genotyped couple.
#'
#'@section Data pre-processing:  Note that the model is designed for autosomes.
#'  Other chromosomes and additional filtering (e.g. bi-allelic markers, markers
#'  and individuals filtering on call rate, coding of missing genotypes, HWE,
#'  etc.) should be performed prior to run RZooRoH with tools such as plink or
#'  bcftools. The model works on an ordered map and ignores SNPs with a null
#'  position.
#'
#'@section RZooRoH functions: The main functions included in the package are
#'  zoodata(), zoomodel() and zoorun(). The zoorun() function can also be
#'  applied to two (phased) haplotypes to obtain IBD probabilities with the same
#'  model. This is possible for haploid data (haploid organism or eventually
#'  specific cases with sex chromosomes) or for diploid individuals, but this
#'  requires then a prior phasing step. The zookin() functions estimates kinship
#'  between pairs of individuals with the ZooRoH model. To that end, it computes
#'  the IBD relationship between the four possible pairs of haplotypes between
#'  the two individuals (this is only possible with phased data). There are also
#'  four functions to plot the results: zooplot_partitioning(), zooplot_hbdseg(),
#'  zooplot_prophbd() and zooplot_individuals(). Eight accessors functions help to
#'  extract the results: realized(), cumhbd(), rohbd(), probhbd(), merge_zres and
#'  update_zres after HBD estimation and cumkin() and predhbd() after using zookin().
#'
#'  You can obtain individual help for each of the functions. By typing for
#'  instance: help(zoomodel) or ? zoomodel.
#'
#'  To run RZooRoH, you must first load your data with the zoodata() function.
#'  It will create a zooin object required for further analysis. Next, you need
#'  to define the model you want to run. You can define a default model by
#'  typing for instance, my.mod <- zoomodel(). Finally, you can run the model
#'  with the zoorun function. You can choose to estimate parameters with
#'  different procedures, estimate global and local homozygous-by-descent (HBD)
#'  probabilities with the Forward-Backward procedure or identify HBD segments
#'  with the Viterbi algorithm. The results are saved in a zres object.
#'
#'  The four plot functions zooplot_partitioning(), zooplot_hbdseg(),
#'  zooplot_prophbd() and zooplot_individuals() use zres objects to make
#'  different graphics. Similarly, the accessor functions help to extract
#'  information from the zres objects (see vignette for more details).
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
#' Mod2L <- zoomodel(K=2,base_rate=10)
#' # Run the model on all individuals.
#' typ.res <- zoorun(Mod2L,typdata)
#' # Observe some results: likelihood, realized autozygosity in different
#' # HBD classes and identified HBD segments.
#' typ.res@modlik
#' typ.res@realized
#' typ.res@hbdseg
#' # Define a model with one HBD and one non-HBD class and run it.
#' Mod1R <- zoomodel(K=1,predefined=FALSE)
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
#'# To run the model on a subset of individuals with 1 thread:
#' \donttest{my.res3 <- zoorun(my.model, example2, ids=c(7,12,16,18), nT = 1)}
#' # Define a smaller model and run it on two individuals.
#' my.mod2 <- zoomodel(K=3,base_rate=10)
#' \donttest{my.res4 <- zoorun(my.mod2, example2, ids=c(9,18))}
#'
#'@import graphics
#'@import stats
#'@import utils
#' @keywords internal
"_PACKAGE"
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
#'and identification of HBD segments (decoding). It is also possible to apply
#'the model to a pair of phased haplotypes with the 'ibd' option.
#'
#'@param zoomodel A valid zmodel object as defined by the zoomodel function. The
#'  model indicates whether rates of exponential distributions are estimated or
#'  predefined, the number of classes, the starting values for mixing
#'  coefficients and rates, the error probabilities. See "zoomodel" for more
#'  details.
#'
#'@param zooin A valid zdata object as obtained by the zoodata function. See
#'  "zoodata" for more details.
#'
#'@param ids An optional argument indicating the individual (its position in the
#'  data file) that must be proceeded. It can also be a vector containing the
#'  list of numbers that must be proceeded. By default, the model runs for all
#'  individuals.
#'
#'@param parameters Specifies whether the parameters are estimated by
#'  optimization with the L-BFGS-B method from the optim function (optional
#'  argument - true by default). If the user doesn't want to estimate the
#'  parameters he must set parameters=FALSE. In that case, the forward-backaward
#'  and Viterbi algorithms are run with the provided parameters.
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
#'  argument - false by default). The Viterbi algorithm performs the decoding
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
#'  individual at each marker are returned when using the Forward-Backward
#'  algorithm (fb option). This is an optional argument that is false by
#'  default.
#'
#'@param nT Indicates the number of threads used when running RZooRoH in
#'  parallel (optional argument - one thread by default).
#'
#'@param optim_method Indicates which method the optim R function will use to
#'  estimate the parameters of the model ("L-BFGS-B" by default). The possible
#'  methods are "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent".
#'  Type "? optim" to have more information. In our experience, the "L-BFGS-B"
#'  method works well but the method achieving the best likelihood is variable
#'  (according to the data sets, the model, the priors, the constraints). The
#'  same goes for the efficiency (speed). When the zoorun does not converge, you
#'  can test with another method. Note that the only method allowing to put
#'  constraints on parameters is "L-BFGS-B" (other methods are unconstrained).
#'
#'@param maxiter Indicates the maximum number of iterations when estimating the
#'  parameters with the R optim function (optional argument - 100 by default).
#'  Iterations are not defined identically across methods. For instance, in one
#'  iteration of the "L-BFGS-B" method, the likelihood of the model, estimated
#'  with the forward algorithm, is evaluated multiple times. So, a value of 100
#'  iterations is good for the "L-BFGS-B" method but larger values are required
#'  for some other algorithms.
#'
#'@param minmix This indicates the minimal value for the mixing coefficients. By
#'  default it is set to 0 with the classical mixkl and kl models
#'  (unconstrained). However, when using the step option or the "Interval"
#'  HBDclass, the values is set to 1e-16 to avoid numerical problems. Note that
#'  constraints are only allowed with the "L-BFGS-B" method from optim.
#'
#'@param maxr This indicates the maximum difference between rates of successive
#'  classes. It is an optional argument set to an arbitrarily large value
#'  (100000000). Adding such constraints might slow down the speed of
#'  convergence and we recommend to run first without this constraint
#'  (constraints are only allowed with the "L-BFGS-B" method from optim).
#'
#'@param ibd A logical indicating whether the function will be used to compute
#'  IBD between pairs of phased haplotypes instead of HBD within individuals. In
#'  that case, the user must provide a matrix with the pairs of haplotypes that
#'  will be analyzed. This option can only be used if phased data are provided
#'  as input with the zformat set to "vcf" or "haps". This is an optional
#'  parameter set to false by default.
#'
#'@param ibdpairs A matrix with four columns, indicating the pair of haplotypes
#'  being analyzed in an IBD analysis. Haplotypes are indicated by two columns,
#'  one column for the id of the individuals and a second column for the
#'  haplotype number within individual (1 or 2). The first and third columns
#'  indicate the id of the individuals carrying the first and second haplotype,
#'  respectively. The second and four columns indicates the haplotype numbers
#'  within the first and second individuals, respectively. With haploid data,
#'  the matrix must have only two columns indicating the number of the first and
#'  second haplotypes from the pairs. This is an optional parameter, the matrix
#'  must be provided only when the ibd option is true.
#'
#'@param haploid This is an optional parameter indicating whether haplotypes
#'  belong to an haploid organisms or chromosome (false by default). It can be
#'  used only in combination with the 'ibd' option and requires phased data as
#'  input with the zformat set to "haps". When haploid is true, then the
#'  ibdpairs matrix has only two columns indicating simply the haplotype
#'  numbers. When haploid is true the number of haplotypes can be uneven while
#'  even numbers are required when haploid is set to false.
#'
#'@param RecTable This is an optional parameter indicating whether a finite
#'  number of genetic distances are used (false by default). This function can
#'  be used only with the "Interval" HBDclass. The "Interval" option can be
#'  slow, in particular if large "intervals" of generations are defined. To
#'  speed up computations, some variables are precomputed for a finite set of
#'  genetic distances, select to cover a broad range of possible values. The
#'  real genetic distance between two genetic markers is then replaced by the
#'  closest value in the table (the difference between the true and used genetic
#'  distances being also lower than 10"\%").
#'
#'@param trim_ad This is an option still under evaluation (for testing only)
#'
#'@param hemiprob This is an option still under evaluation (for testing only)
#'
#'@return The function return a zoores object with several slots accesses by the
#'  "@" symbol. The three main results are zoores@@realized (the matrix with
#'  partitioning of the genome in different HBD classes for each individual),
#'  zoores@@hbdseg (a data frame with identified HBD segments) and zoores@@hbdp
#'  (a list of matrices with HBD probabilities per SNP and per class).
#'
#'  Here is a list with all the slots and their description:
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
#'  parameters - the number of calls of the forward algorithm (per individual),
#'  \item zoores@@modlik a vector containing the likelihood of the model for
#'   each individual,
#'  \item zoores@@modbic a vector containing the value of the BIC for each
#'  individual,
#'  \item zoores@@realized a matrix with estimated realized autozygosity per HBD
#'  class (columns) for each individual (rows). These values are obtained with
#'  the Forward-Backward algorithm - fb option),
#'  \item zoores@@hbdp a list of matrices with the
#'  local probabilities to belong to an underlying hidden state (computed for
#'  every class and every individual). Each matrix has one row per class and
#'  one column per snp. To access the matrix from individual i, use the
#'  brackets "[[]]", for instance zoores@@hbdp[[i]],
#'  \item zoores@@hbdseg a data frame with the list
#'  of identified HBD segments with the Viterbi algorithm (the columns are the
#'  individual number, the chromosome number, the first and last SNP of the
#'  segment, the positions of the first and last SNP of the segment, the number
#'  of SNPs in the segment, the length of the segment, the HBD state of the
#'  segment). In case of IBD, the first column indicates the number of the pair
#'  of haplotypes,
#'  \item zoores@@optimerr a vector indicating whether optim ran with or without
#'  error (0 indicates successful completion / 1 indicates that the iteration
#'  limit has been reached / 51 and 52 indicate warnings from the "L-BFGS-B" method
#'  / 99 indicates numerical problem). See the optim R function for more
#'  details.
#'  \item zoores@@sampleids is a vector with the names of the
#'  samples (when provided in the zooin object through the zoodata function).
#'  When an IBD analysis is run, it indicates the pairs of haplotypes separated
#'  by "_". For diploid individuals, it combines individual ID, haplotype ID (1
#'  or 2) for the two individuals.}
#'
#'@examples
#'
#' # Start with a small data set with six individuals and external frequencies.
#' freqfile <- (system.file("exdata","typsfrq.txt",package="RZooRoH"))
#' typfile <- (system.file("exdata","typs.txt",package="RZooRoH"))
#' frq <- read.table(freqfile,header=FALSE)
#' typ <- zoodata(typfile,supcol=4,chrcol=1,poscol=2,allelefreq=frq$V1)
#' # Define a model with two HBD classes with rates equal to 10 and 100.
#' Mod2L <- zoomodel(K=2,base_rate=10)
#' # Run the model on all individuals.
#' typ.res <- zoorun(Mod2L, typ)
#' # Observe some results: likelihood, realized autozygosity in different
#' # HBD classes and identified HBD segments.
#' typ.res@modlik
#' typ.res@realized
#' typ.res@hbdseg
#' # Define a model with one HBD and one non-HBD class and run it.
#' Mod1R <- zoomodel(K=1,predefined=FALSE)
#' typ2.res <- zoorun(Mod1R, typ)
#' # Print the estimated rates and mixing coefficients.
#' typ2.res@krates
#' typ2.res@mixc

#' # Get the name and location of a second example file and load the data:
#' myfile <- (system.file("exdata","genoex.txt",package="RZooRoH"))
#' ex2 <- zoodata(myfile)
#' # Run RZooRoH to estimate parameters on your data with the 1 HBD and 1 non-HBD
#' # class (parameter estimation with optim).
#' my.mod1R <- zoomodel(predefined=FALSE,K=1,krates=c(10))
#' \donttest{my.res <- zoorun(my.mod1R, ex2, fb = FALSE, vit = FALSE)}
#' # The estimated rates and mixing coefficients:
#' \donttest{my.res@mixc}
#' \donttest{my.res@krates}
#' # Run the same model and run the Forward-Backward alogrithm to estimate
#' # realized autozygosity and the Viterbi algorithm to identify HBD segments:
#' \donttest{my.res2 <- zoorun(my.mod1R, ex2)}
#' # The table with estimated realized autozygosity:
#' \donttest{my.res2@realized}
#' # Run a model with 3 layers (3 HBD classes / 1 non-HBD class) and estimate
#' # the rates of HBD classes with one thread:
#' my.mod3L <- zoomodel(predefined=FALSE,K=3,krates=c(16,64,256))
#' \donttest{my.res3 <- zoorun(my.mod3L, ex2, fb = FALSE, vit = FALSE, nT =1)}
#' # The estimated rates for the 3 classes and the 20 individuals:
#' \donttest{my.res3@krates}
#' # Run a model with 4 layers and predefined rates.
#' # The model is run only for a subset of four selected individuals.
#' # The parameters are estimated, the Forward-Backward alogrithm is ued to
#' # estimate realized autozygosity and the Viterbi algorithm to identify
#' # HBD segments. One thread is used.
#' mix4L <- zoomodel(K=4,base=10)
#' \donttest{my.res4 <- zoorun(mix4L,ex2,ids=c(7,12,16,18), nT = 1)}
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

zoorun <- function(zoomodel, zooin, ids = NULL, parameters = TRUE, fb = TRUE, vit = TRUE, localhbd = FALSE,
                   nT = 1, optim_method = "L-BFGS-B", maxiter = 100, minmix = 1, maxr = 100000000,
                   ibd = FALSE, ibdpairs = NULL, haploid = FALSE, RecTable = FALSE, trim_ad = FALSE,
                   hemiprob = 0){

  if(is.null(ids) & !ibd){ids=seq(from=1,to=zooin@nind)}
  if(ibd & is.null(ibdpairs)){
    print("Pairs of haplotypes were not provided - analayzing all pairs.")
    if(!haploid){
      nhap=2*zooin@nind
      npairs=0.5*nhap*(nhap - 1)
      ibdpairs=matrix(0,npairs,4)
      numpair=0
      for (i1 in 1:(zooin@nind-1)){
        numpair=numpair+1
        ibdpairs[numpair,1]=i1;ibdpairs[numpair,2]=1
        ibdpairs[numpair,3]=i1;ibdpairs[numpair,4]=2
        for (j1 in 1:2){
          for (i2 in (i1+1):zooin@nind){
            for (j2 in 1:2){
              numpair=numpair+1
              ibdpairs[numpair,1]=i1
              ibdpairs[numpair,2]=j1
              ibdpairs[numpair,3]=i2
              ibdpairs[numpair,4]=j2
            }
          }
        }
      }
      numpair=numpair+1
      ibdpairs[numpair,1]=zooin@nind;ibdpairs[numpair,2]=1
      ibdpairs[numpair,3]=zooin@nind;ibdpairs[numpair,4]=2
    }
    if(haploid){
      npairs=0.5*zooin@nind*(zooin@nind - 1)
      ibdpairs=matrix(0,npairs,2)
      numpair=0
      for (i1 in 1:(zooin@nind-1)){
        for (i2 in (i1+1):zooin@nind){
          numpair=numpair+1
          ibdpairs[numpair,1]=i1
          ibdpairs[numpair,2]=i2
        }
      }
    }
  }
  if(ibd){ids=seq(from=1,to=nrow(ibdpairs))}
  if(is.zmodel(zoomodel) == FALSE){print("First element is not a valid model - see help for zoomodel function.")}
  if(is.zdata(zooin) == FALSE){print("First element is not a valid data set - see help for zoodata function.")}
  if(haploid & !ibd){print("Haploid data can only be used with the IBD option.")}
  zrates <- zoomodel@krates
  zmix <- zoomodel@mix_coef
  gerr <- zoomodel@err
  seqerr <- zoomodel@seqerr
  if(minmix==1){ #### default, the user did not specify it's own value
    if(zoomodel@typeClass=="Interval"){minmix=1e-16}
    if(zoomodel@typeClass=="SingleRate"){minmix=0}
    if(zoomodel@typeModel == "step_mixkl" || zoomodel@typeModel == "STEP_MIXKL"){minmix=1e-16}
  }
  HBDclass <- zoomodel@typeClass
  XM <- zoomodel@XM

  if(.Platform$OS.type == "windows" & nT > 1){
    warning("Multithreading is not supported under windows. Only one thread will be used.\n")}
  if(.Platform$OS.type != "windows"){registerDoParallel(cores= nT)}

  if(zoomodel@typeModel != "mixkl" & zoomodel@typeModel != "MIXKL"
     & zoomodel@typeModel != "kl" & zoomodel@typeModel != "KL"
     & zoomodel@typeModel != "step_mixkl" & zoomodel@typeModel != "STEP_MIXKL"){
    print("The model type must be kl, mixkl or step_mixkl")
  }

  if(parameters){method ="opti";opti=TRUE}
  if(!parameters){method="no";opti=FALSE}

  zoores <- new("zres")
  zoores@nind = length(ids)
  zoores@ids = ids
  zoores@mixc = matrix(0,zoores@nind,length(zmix))
  zoores@krates = matrix(0,zoores@nind,length(zrates))

  if(fb){zoores@realized = matrix(0,zoores@nind,length(zmix))}

  i=NULL
  if(.Platform$OS.type == 'windows'){
    if(zoomodel@typeModel == "mixkl" || zoomodel@typeModel == "MIXKL"){
      xx <- foreach (i=1:length(ids), .combine=c) %do% {
        id <- ids[i]
        if(ibd){ibdpair=ibdpairs[i,]}
        run_mixkl(zooin,id, zrates, zmix, opti, fb, vit, gerr, seqerr, hemiprob, maxr, minmix, maxiter,
                  localhbd,ibd,ibdpair,haploid,HBDclass,RecTable,optim_method,trim_ad)
      }
    }
    if(zoomodel@typeModel == "kl" || zoomodel@typeModel == "KL"){
      xx <- foreach (i=1:length(ids), .combine=c) %do% {
        id <- ids[i]
        if(ibd){ibdpair=ibdpairs[i,]}
        run_kl(zooin,id, zrates, zmix, opti, fb, vit, gerr, seqerr, hemiprob, maxr, minmix, maxiter,
               localhbd,ibd,ibdpair,haploid,optim_method,trim_ad)
      }
    }
    if(zoomodel@typeModel == "step_mixkl" || zoomodel@typeModel == "STEP_MIXKL"){
      xx <- foreach (i=1:length(ids), .combine=c) %do% {
        id <- ids[i]
        if(ibd){ibdpair=ibdpairs[i,]}
        run_fun_mixkl(zooin,id, zrates, zmix, XM, opti, fb, vit, gerr, seqerr, hemiprob, maxr, minmix, maxiter,
                      localhbd, ibd, ibdpair, haploid, optim_method, trim_ad)
      }
    }
  }else{
    if(zoomodel@typeModel == "mixkl" || zoomodel@typeModel == "MIXKL"){
      xx <- foreach (i=1:length(ids), .combine=c) %dopar% {
        id <- ids[i]
        if(ibd){ibdpair=ibdpairs[i,]}
        run_mixkl(zooin,id, zrates, zmix, opti, fb, vit, gerr, seqerr, hemiprob, maxr, minmix, maxiter,
                  localhbd, ibd, ibdpair, haploid, HBDclass, RecTable, optim_method,trim_ad)
      }
    }
    if(zoomodel@typeModel == "kl" || zoomodel@typeModel == "KL"){
      xx <- foreach (i=1:length(ids), .combine=c) %dopar% {
        id <- ids[i]
        if(ibd){ibdpair=ibdpairs[i,]}
        run_kl(zooin,id, zrates, zmix, opti, fb, vit, gerr, seqerr, hemiprob, maxr, minmix, maxiter,
               localhbd,ibd,ibdpair,haploid,optim_method,trim_ad)
      }
    }
    if(zoomodel@typeModel == "step_mixkl" || zoomodel@typeModel == "STEP_MIXKL"){
      xx <- foreach (i=1:length(ids), .combine=c) %dopar% {
        id <- ids[i]
        if(ibd){ibdpair=ibdpairs[i,]}
        run_fun_mixkl(zooin,id, zrates, zmix, XM, opti, fb, vit, gerr, seqerr, hemiprob, maxr, minmix, maxiter, localhbd,
                      ibd,ibdpair,haploid,optim_method,trim_ad)
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
  if(fb){
    zoores@realized <- do.call("rbind",lapply(xx,slot,"realized"))
    if(localhbd){zoores@hbdp <- lapply(xx,slot,"hbdp")}
  }
  if(vit){zoores@hbdseg <- do.call("rbind", lapply(xx,slot,"hbdseg"))}
  if(!ibd){zoores@sampleids <- zooin@sample_ids[zoores@ids]}
  if(ibd & !haploid){zoores@sampleids <- paste(zooin@sample_ids[ibdpairs[zoores@ids,1]],
                                    ibdpairs[zoores@ids,2],
                                    zooin@sample_ids[ibdpairs[zoores@ids,3]],
                                    ibdpairs[zoores@ids,4],sep="_")}
  if(ibd & haploid){zoores@sampleids <- paste(zooin@sample_ids[ibdpairs[zoores@ids,1]],
                                               ibdpairs[zoores@ids,2],sep="_")}

  return(zoores)
}

#### run mixkl model

run_mixkl <- function(zooin,id, zrates, zmix, opti = TRUE, fb = FALSE, vit = FALSE,
                      gerr = 0.001, seqerr = 0.001, hemiprob = 0.000, maxr = 100000000, minmix = 1e-16,
                      maxiter = 1000, localhbd = FALSE, ibd = FALSE, ibdpair = NULL, haploid = FALSE,
                      hbdc_type = "SingleRate", RecTable = FALSE, optim_method = "L-BFGS-B",
                      trim_ad = FALSE){
  K <- length(zmix)+1 #### number of states
  pem <- pemission(zooin, id, gerr, seqerr, hemiprob, zformat = zooin@zformat, ibd, ibdpair, haploid)
  skiperr=TRUE

  if(opti){
    zpar <- array(0,K-1)
    for (j in 1:(K-1)){zpar[j]=log(zmix[j]/(1-zmix[j]))}
    if(minmix>0)myminpar=log(minmix/(1-minmix))
    if(minmix==0)myminpar=-10**100
    niter=NULL;niter <<- 0
    checkoptim=NULL;checkoptim <<- 1

    if(zooin@zformat=="ad" & trim_ad){
      tad <- array(0,nrow(pem))
      for (i in 1:zooin@nchr){tad[zooin@chrbound[i,1]:zooin@chrbound[i,2]]=i}
      tad[(zooin@genos[,(2*id-1)]+zooin@genos[,(2*id)]) <= 1]=0
      pem2=pem[tad > 0,]
      nsnps2 = nrow(pem2)
      print(c("Number of SNPs with cover > 1 :: ",nrow(pem2)))

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
        tryCatch(optires <- optim(zpar, lik_mixkl, nchr = nchr2 ,nsnps = nsnps2, chrbound = chrbound2,
                                  bp = bp2, pem = pem2, zrs = zrates,
                                  HBD_classes = hbdc_type, RecTable = RecTable, method=optim_method,
                                  lower=rep(myminpar,(K-1)),upper=rep(-myminpar,(K-1)),
                                  control = list(trace=TRUE, REPORT=10000,maxit=maxiter)),
                 error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
      }else{
        tryCatch(optires <- optim(zpar, lik_mixkl, nchr = nchr2 ,nsnps = nsnps2, chrbound = chrbound2,
                                  bp = bp2, pem = pem2, zrs = zrates,
                                  HBD_classes = hbdc_type, RecTable = RecTable, method=optim_method,
                                  control = list(trace=TRUE, REPORT=10000,maxit=maxiter)),
                 error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
      }
    }else{
    if(optim_method=="L-BFGS-B" & minmix > 0){
     tryCatch(optires <- optim(zpar, lik_mixkl, nchr = zooin@nchr ,nsnps = zooin@nsnps, chrbound = zooin@chrbound,
                               bp = zooin@bp, pem = pem, zrs = zrates,
                               HBD_classes = hbdc_type, RecTable = RecTable, method=optim_method,
                               lower=rep(myminpar,(K-1)),upper=rep(-myminpar,(K-1)),
                              control = list(trace=TRUE, REPORT=10000,maxit=maxiter)),
              error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
    }else{
      tryCatch(optires <- optim(zpar, lik_mixkl, nchr = zooin@nchr ,nsnps = zooin@nsnps, chrbound = zooin@chrbound,
                                bp = zooin@bp, pem = pem, zrs = zrates,
                                HBD_classes = hbdc_type, RecTable = RecTable, method=optim_method,
                                control = list(trace=TRUE, REPORT=10000,maxit=maxiter)),
             error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
    }}
    if(.GlobalEnv$checkoptim ==0){skiperr =FALSE} ### NEW
    if(.GlobalEnv$checkoptim == 1){
      zmix <- exp(optires$par[1:(K-1)])/(1+exp(optires$par[1:(K-1)]))
    }
  }
  if(fb & skiperr){
    outfb <- fb_lay(zooin,id,pem,zrates,zmix,localhbd,hbdc_type,RecTable)
  }
  if(vit & skiperr){
    outvit <- viterbi_lay(zooin,id,pem,zrates,zmix,hbdc_type,RecTable)
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

#### run kl model

run_kl <- function(zooin,id, zrates, zmix, opti = TRUE, fb = FALSE, vit = FALSE,
                   gerr = 0.001, seqerr = 0.001, hemiprob = 0.000, maxr = 100000000, minmix = 1e-16,
                   maxiter = 1000, localhbd = FALSE, ibd = FALSE, ibdpair = NULL, haploid = FALSE,
                   optim_method = "L-BFGS-B", trim_ad = FALSE){
  pem <- pemission(zooin, id, gerr, seqerr, hemiprob, zformat = zooin@zformat, ibd, ibdpair, haploid)
  K <- length(zmix)+1 #### K remains the number of states
  skiperr=TRUE

  if(opti){
    if(K > 2){
      zpar <- array(0,2*K-2) #### New, we don't need the rates for non-hbd segments
      zpar[1] <- log(zrates[1]-1)
      for(j in 2:(K-1)){zpar[j] <- log(zrates[j]-zrates[j-1])}
      for (j in 1:(K-1)){zpar[K+j-1]=log(zmix[j]/(1-zmix[j]))}
    } else { #### 1R model
      zpar <- array(0,2)
      zpar[1] <- log(zrates[1]-1) ### the common rate
      zpar[2] <- log(zmix[1]/(1-zmix[1]))
    }

    if(minmix>0)myminpar=log(minmix/(1-minmix))
    if(minmix==0)myminpar=-10**100
    niter=NULL;niter <<- 0
    checkoptim=NULL;checkoptim <<- 1

    if(K == 2){
      lowpar <- c(log(1e-16), myminpar)
      uppar <- c(log(maxr), -myminpar)
    }
    if(K > 2){
      lowpar <- c(rep(log(1e-16),(K-1)), rep(myminpar,(K-1)))
      uppar <- c(rep(log(maxr),(K-1)), rep(-myminpar,(K-1)))
    }

    if(zooin@zformat=="ad" & trim_ad){
      tad <- array(0,nrow(pem))
      for (i in 1:zooin@nchr){tad[zooin@chrbound[i,1]:zooin@chrbound[i,2]]=i}
      tad[(zooin@genos[,(2*id-1)]+zooin@genos[,(2*id)]) <= 1]=0
      pem2=pem[tad > 0,]
      nsnps2 = nrow(pem2)
      print(c("Number of SNPs with cover > 1 :: ",nrow(pem2)))

      chrnum <- tad[tad > 0]
      snpnum <- seq(1:nsnps2)
      nchr2 <- zooin@nchr
      chrbound2 <- matrix(0, nchr2, 2)
      for (i in 1:nchr2){
        chrbound2[i,1]=min(snpnum[chrnum==i])
        chrbound2[i,2]=max(snpnum[chrnum==i])
      }
      bp2 <- zooin@bp[tad > 0]

      if((maxr < 100000000 | minmix > 0) & optim_method=="L-BFGS-B"){
        tryCatch(optires <- optim(zpar, lik_kl, nchr =nchr2, nsnps= nsnps2, chrbound = chrbound2,
                                  bp = bp2, pem = pem, method="L-BFGS-B", lower = lowpar, upper = uppar,
                                  control = list(trace=TRUE, REPORT=10000, maxit=maxiter)),
                 error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
      }else{
        tryCatch(optires <- optim(zpar, lik_kl, nchr =nchr2, nsnps= nsnps2, chrbound = chrbound2,
                                  bp = bp2, pem = pem, method=optim_method,
                                  control = list(trace=TRUE, REPORT=10000, maxit=maxiter)),
                 error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
      }
    }else{
    if((maxr != 100000000 | minmix > 0) & optim_method=="L-BFGS-B"){
      tryCatch(optires <- optim(zpar, lik_kl, nchr =zooin@nchr, nsnps= zooin@nsnps, chrbound = zooin@chrbound,
                                bp = zooin@bp, pem = pem, method="L-BFGS-B", lower = lowpar, upper = uppar,
                                control = list(trace=TRUE, REPORT=10000, maxit=maxiter)),
               error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
    }else{
      tryCatch(optires <- optim(zpar, lik_kl, nchr =zooin@nchr, nsnps= zooin@nsnps, chrbound = zooin@chrbound,
                                bp = zooin@bp, pem = pem, method=optim_method,
                                control = list(trace=TRUE, REPORT=10000, maxit=maxiter)),
               error=function(e){cat("Warning, error returned by optim - try other algorithm ::",id,"\n"); .GlobalEnv$checkoptim <- 0})
    }}
    if(.GlobalEnv$checkoptim ==0){skiperr = FALSE} ### NEW

    if(.GlobalEnv$checkoptim == 1){
      if(K > 2){
        zrates[1] <- 1 + exp(optires$par[1])
        for (j in 2:(K-1)){zrates[j] <- zrates[j-1] + exp(optires$par[j])}
        zmix <- exp(optires$par[(K):(2*K-2)])/(1+exp(optires$par[(K):(2*K-2)])) ### 1:(K-1) rates, K:(2K-2) mixing
      } else {
        zrates <- 1 + exp(optires$par[1])
        zmix[1] <- 1/(1+exp(-optires$par[2]))
      }
    }
  }

  if(fb & skiperr){
    outfb <- fb_lay(zooin,id,pem,zrates,zmix,localhbd)
  }
  if(vit & skiperr){
    outvit <- viterbi_lay(zooin,id,pem,zrates,zmix)
  }

  ##### renvoie toujours les paramètres
  zresu <- new("oneres")
  loglik <- 0
  bicv <- 0
  zrates <- round(zrates,3)
  if(opti){
    if (.GlobalEnv$checkoptim == 1){
      loglik <- -optires$value ### optim miminizes -log(lik)
      if(K > 2)bicv <- -2 * loglik + log(zooin@nsnps)*(2*K-1)
      if(K == 2)bicv <- -2 * loglik + 2*log(zooin@nsnps)
    }}

  zresu@mixc <- zmix
  zresu@krates <- zrates
  zresu@modlik <- loglik
  zresu@modbic <- bicv
  zresu@num <- id
  zresu@niter <- .GlobalEnv$niter
  zresu@optimerr <- -1
  if(opti & .GlobalEnv$checkoptim==0){zresu@optimerr <- 99}
  if(opti & .GlobalEnv$checkoptim==1){zresu@optimerr <- optires$convergence}

  if(!skiperr){
    zresu@krates <- rep(0,K)
    zresu@mixc <- rep(0,K)
  }

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

