setClass(Class = "kres",
         representation(npairs = "numeric", pairs = "vector", realized = "matrix", krates = "vector",
                        ibdp = "list", ibdseg ="data.frame", sampleids = "vector", haplotype_ids = "vector")
)

is.kres <- function (x)
{
   res <- (is(x,"kres") & validObject(x))
   return(res)
}

#### estimate kinship for pairs of individuals

#'Use the ZooRoH model to estimate kinship between pairs of individuals
#'
#'Apply the defined model to estimate kinship for pairs of individuals by running
#'the ZooRoH model to estimate IBD probabilities between the four possible pairs
#'of chromosome of the two individuals (one chromosome from each individual). As
#'for zoorun this includes parameter estimation, computation of realized IBD,
#'IBD probabilities, and identification of IBD segments. Options are similar to
#'the zoorun function.
#'
#'@param zoomodel A valid zmodel object as defined by the zoomodel function. The
#'  model indicates whether rates of exponential distributions are estimated or
#'  predefined, the number of classes, the starting values for mixing
#'  coefficients and rates, the error probabilities. See "zoomodel" for more
#'  details. zookin can not be run with a KL model because four different rate
#'  parameters would be estimated for each pair of haplotypes, making
#'  interpretation difficult.
#'
#'@param zooin A valid zdata object as obtained by the zoodata function. See
#'  "zoodata" for more details.
#'
#'@param parameters Specifies whether the parameters are estimated by
#'  optimization with the L-BFGS-B method from the optim function (optional
#'  argument - true by default). If the user doesn't want to estimate the
#'  parameters he must set parameters=FALSE. In that case, the forward-backaward
#'  and Viterbi algorithms are run with the provided parameters.
#'
#'@param fb A logical indicating whether the forward-backward algorithm is run
#'  (optional argument - true by default). The Forward-Backward algorithm
#'  estimates the local probabilities to belong to each IBD or non-IBD class. By
#'  default, the function returns only the IBD probabilities for each class,
#'  averaged genome-wide, and corresponding to the realized autozygosity
#'  associated with each class. To obtain HBD probabilities at every marker
#'  position, the option localhbd must be set to true (this generates larger
#'  outputs).
#'
#'@param vit A logical indicating whether the Viterbi algorithm is run (optional
#'  argument - false by default). The Viterbi algorithm performs the decoding
#'  (determining the underlying class at every marker position). Whereas the
#'  Forward-Backward algorithms provide IBD probabilities (and how confident a
#'  region can be declared IBD), the Viterbi algorithm assigns every marker
#'  position to one of the defined classes (IBD or non-IBD). When informativity
#'  is high (many SNPs per IBD segments), results from the Forward-Backward and
#'  the Viterbi algorithm are very similar. The Viterbi algorithm is best suited
#'  to identify IBD segments. To estimate realized kinship and determine IBD
#'  status of a position, we recommend to use the Forward-Backward algorithm
#'  that better reflects uncertainty.
#'
#'@param localhbd A logical indicating whether the IBD probabilities for each
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
#'  iterations is good for the "L-BFGS-B" method but large values are required
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
#'@param kinpairs A matrix with two columns, indicating the pairs of individuals
#'  being analyzed. This information is required for a zookin analysis.
#'
#'@param RecTable This is an optimal parameter indicating whether a finite
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
#'@return The function return a kinres object with several slots accesses by the
#'  "@" symbol. It is very similar to the zoores object. The three main results
#'  are kinres@@realized (the matrix with partitioning of the genome in
#'  different IBD classes for each pair of individual), These values are scale
#'  like coancestry coefficient and represent the probability that two alleles
#'  sampled in the two individuals are IBD (and in that length class). The sum
#'  of the realized values in all IBD classes correspond to the kinship.
#'  kinres@@hbdseg (a data frame with identified IBD segments) and zoores@@hbdp
#'  (a list of matrices with IBD probabilities per SNP and per class). This gives
#'  at one marker position, the probability that two haplotypes, one sampled in
#'  each individual, are IBD (and in that length class) at that position. It
#'  corresponds also to the predicted HBD values in their future progeny.
#'
#'  Some slots present in zoores objects are not reported because for each pair
#'  of individuals we would then have four values. If you are interested in these
#'  values for each pair of haplotypes, you can run zoorun with the ibd option.
#'  These slots are @@mixc, @@modlik, @@modbic, @@niter, @@optimerr.
#'
#'  Here is a list with all the slots and their description:
#'
#' \enumerate{
#' \item kinres@@npairs the number of pairs of individuals in the analysis,
#' \item kinres@@ids a matrix containing the numbers of the analyzed pairs of
#' individuals (their position in the data file),
#'  \item kinres@@krates the rates for the exponential distributions
#'  associated with each IBD or non-IBD class for all individuals,
#'  \item kinres@@realized a matrix with estimated realized IBD per layer
#'   (columns) for each pair of individuals (rows). These values are obtained with
#'  the Forward-Backward algorithm - fb option - and averaging over the four
#'  possible haplotype pairs (between the two individuals),
#'  \item kinres@@ibdp a list of matrices with the local probabilities of IBD
#'  in different layers (computed for every class and every pair of individuals).
#'  The IBD probability is the probability that two alleles sampled at that locus
#'  in the two individuals (one in each) are IBD with respect to that layer (e.g.,
#'  the IBD segment must correspond to the layer). Each matrix has one row per
#'  layer / class and one column per snp. To access the matrix for pair (of individuals)
#'   i, use the brackets "[[]]", for instance kinres@@hbdp[[i]],
#'  \item zookin@@ibdseg a data frame with the list of identified IBD segments
#'   with the Viterbi algorithm (the columns are the haplotype pair number,
#'   the chromosome number, the first and last SNP of the segment, the positions
#'   of the first and last SNP of the segment, the number of SNPs in the segment,
#'   the length of the segment, the HBD state of the segment),
#'  \item kinres@@sampleids is a vector with the names of the pair of
#'  samples (when provided in the zooin object through the zoodata function).
#'  The ids from the two individuals are separated by a "_".
#'  \item kinres@@haplotype_ids is a vector with the information on
#'  haplotype pairs in the analysis and used in the ibdseg table.}

#'@useDynLib RZooRoH, .registration = TRUE
#'@export

zookin <- function(zoomodel, zooin, parameters = TRUE, fb = TRUE, vit = TRUE, localhbd = FALSE,
                   nT = 1, optim_method = "L-BFGS-B", maxiter = 1000, minmix = 1, maxr = 100000000,
                   kinpairs = NULL, RecTable = FALSE, trim_ad = FALSE, hemiprob = 0){

  if(zoomodel@typeModel=="kl"){stop(paste("zookin is not compatible with KL models."))}
  kinres <- new("kres")
  kinres@krates = zoomodel@krates

  if(is.null(kinpairs)){
    print("Pairs of individuals were not provided - analayzing all pairs.")
    npairs=0.5*zooin@nind*(zooin@nind - 1)
    kinpairs=matrix(0,npairs,2)
    numpair=0
    for (i1 in 1:(zooin@nind-1)){
      for (i2 in (i1+1):zooin@nind){
        numpair=numpair+1
        kinpairs[numpair,1]=i1
        kinpairs[numpair,2]=i2
      }
    }
  }
  kinres@npairs = nrow(kinpairs)
  kinres@pairs = kinpairs

  npairs=nrow(kinpairs)
  nibdpairs=nrow(kinpairs)*4
  allpairs <- matrix(0,nibdpairs,4)
  for(i in 1:npairs){
   mypairs <- matrix(0,4,4)
   mypairs[,1]=kinpairs[i,1]
   mypairs[,3]=kinpairs[i,2]
   mypairs[,2]=c(1,1,2,2)
   mypairs[,4]=c(1,2,1,2)
   i1=4*i-3
   i2=4*i
   allpairs[i1:i2,]=mypairs
  }

  kinres@sampleids <- paste(zooin@sample_ids[kinpairs[,1]],
                            zooin@sample_ids[kinpairs[,2]],
                            sep="_")

  tt <- zoorun(zoomodel=zoomodel,zooin=zooin,parameters=parameters,fb=fb,vit=vit,localhbd=localhbd,
               nT=nT,optim_method=optim_method,maxiter=maxiter,minmix=minmix,maxr=maxr,
               ibd=TRUE,ibdpairs=allpairs,haploid=FALSE,RecTable = RecTable, trim_ad=FALSE,hemiprob=hemiprob)

   NL=length(zoomodel@mix_coef)
   realizedkin <- matrix(0,npairs,(NL+1))
   for (i in 1:npairs){
      realizedkin[i,] <- 0.25*(tt@realized[4*i-3,]+tt@realized[4*i-2,]+tt@realized[4*i-1,]+tt@realized[4*i,])
   }

   kinres@realized=realizedkin[,(1:NL)]

   kinloc=list()
   if(localhbd){
      for (i in 1:npairs){
         kinloc[[i]]  = 0.25*(tt@hbdp[[4*i-3]] + tt@hbdp[[4*i-2]] + tt@hbdp[[4*i-1]] + tt@hbdp[[4*i]])
      }
      kinres@ibdp=kinloc
   }

   if(vit){kinres@ibdseg=tt@hbdseg}
   kinres@haplotype_ids = tt@sampleids
   return(kinres)
  }


