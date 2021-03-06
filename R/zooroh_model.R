
setClass(Class = "zmodel",
         representation(typeModel = "character", mix_coef = "vector", krates = "vector",
                        err = "numeric", seqerr = "numeric")
)

is.zmodel <- function (x)
{
  res <- (is(x,"zmodel") & validObject(x))
  return(res)
}

#'Define the model for the RZooRoH
#'
#'Help the user to create a model for RZooRoH, including default parameters. The
#'output is a zmodel object necessary to run RZooRoH.
#'
#'@param predefined Logical (TRUE or FALSE) to define whether rates of HBD and
#'  non-HBD-classes will be estimated by the model ("kr" model) or whether the
#'  rates of these classes are fixed and pre-defined by the user ("mixkr"
#'  model). The default value is "predefined = TRUE".
#'
#'@param K The number of HBD and non-HBD classes. There are always K-1 HBD
#'  classes and one non-HBD class. By default, K is set to 10 but this is not
#'  optimal for all data sets. Hence, we recommend to the users to select their
#'  own value. If K is set to 2 and rates are estimated, RZooRoH will use the
#'  same rate for the HBD and the non-HBD class (so-called 1R model).
#'
#'@param mix_coef The starting value for the mixing coefficients for all HBD and
#'  non-HBD classes. The mixing coefficients determine the frequency of the
#'  segments from different classes, they determine the probability to start a
#'  new segment in a given class when a segment ends. The mixing coefficients
#'  should sum to 1 and the function expects K mixing coefficients. The default
#'  values are 0.01 for HBD classes and the value for the non-HBD class is such
#'  that all mixing coefficients sum to 1. In case the parameters are not
#'  estimated (e.g. when running the forward-backward or the Viterbi algorithm
#'  alone), these are the mixing coefficients used by the RZooRoH model.
#'
#'@param base_rate is a integer used to define the rates of successive HBD
#'  classes (see krates below). This parameter is most useful when using a
#'  mixkr model with predefined rates. The rate of each HBD class will be equal
#'  to the base_rate raised to the exponent k (the class number). The non-HBD
#'  class will have the same rate as the last HBD class. For instance, with a
#'  base_rate of 2 and five classes, we have the following rates: 2, 4, 8, 16
#'  and 16. Similarly, with a base_rate of 10 and four classes, we have 10, 100,
#'  1000 and 1000. With this method, more HBD classes are defined for more
#'  recent ancestors (for which we have more information to estimate R) and less
#'  for ancient HBD classes (it doesn't make sense to try to distinguish R =
#'  1000 from R = 1010). In addition, since the expected length of HBD segments
#'  is expected to be approximately 1/R, the ratio between successive expected
#'  HBD lengths remains the same. This ratio also determines the ability of the
#'  model to distinguish segments from distinct classes. By keeping the ratio
#'  constant, the aptitude to discriminate between HBD classes is also constant.
#'  In addition, this method allows to cover a wide range of generations in the
#'  past, with more emphasis on recent ancestors. The default value for the
#'  base_rate is 2.
#'
#'@param krates Is an array with a rate for each HBD and non-HBD class. The
#'  function expects K positive rates. These rates are parameters of the
#'  exponential distribution that together with the distance in centimorgans
#'  defines the probability to end a HBD segments between two markers. Each HBD
#'  class has a distinct rate. Therefore, the expected length of HBD classes is
#'  defined by the rates. The expected length if equal to 1/R. These krates are
#'  associated with the age of the common ancestor of the HBD segment. The rate
#'  is approximately equal to the size of the inbreeding loop (twice the number
#'  of generations to the common ancestor) when the map is given in Morgans. By
#'  default, the rates are defined by the base_rate parameter (2, 4, 8, 16,
#'  ...).
#'
#'@param err Indicates the error term, the probability to observe an
#'  heterozygous genotype in a HBD segment. The genotype could be heterozygous
#'  due to a mutation occuring on the path to the common ancestor. It can also
#'  be associated with a genotype calling error or a technical error. In case GP
#'  or GL formats are used (with genotyped probabilities or phred scores) or
#'  when an AD format is used (based on read counts), this error term still
#'  represents the probability to observe an heterozygous genotype in a HBD
#'  segment. When an heterozygous genotype was called with a probability equal
#'  to 1.00, this heterozygosity in an HBD track might be associated to a
#'  mutation or to errors not accounted for by the model used to estimate the
#'  genotype probabilities (e.g., GATK). The emission probability to observe a
#'  heterozygous genotype in an HBD class will never go below the error term.
#'  The default value is 0.001.
#'
#'@param seqerr This parameters is used only with the AD format. In the AD
#'  format the user gives the number of reads for both alleles. A simple model
#'  is then used to estimate the genotype probabilities based on the read
#'  counts. In that model, the seqerr represents the probability to have a
#'  sequencing error in one read. The default value is 0.001.
#'
#'@param layers Logical (TRUE or FALSE) - When true, this parameter indicates
#'  that the data is modeled as mosaic of HBD and non-HBD classes at different
#'  levels. At each level, HBD and non-HBD classes have the same rate (the same
#'  expected length). Non-HBD classes are subsequently modelled as mosaic of
#'  non-HBD segments and HBD segments from more ancient generations (from
#'  smaller sizes). At each level, the mixing coefficients can be interpreted as
#'  the inbreeding coefficient at that level (TRUE by default). This model corresponds
#'  to the Nested 1R model (N1R).
#'
#'@return The function return an object that defines a model for RZooRoH and
#'  incuding the following elements: zmodel@@typeModel equal to "kr", "mixkr" or
#'  "mixkl" according to the selected model, zmodel@@mix_coef an array with
#'  mixing coefficients, zmodel@@krates an array with the rates of the HBD and
#'  non-HBD classes, zmodel@@err the parameter defining the probability to
#'  observe an heterozygous genotype in an HBD class, and zmodel@@seqerr the
#'  parameter defining the probability of sequencing error per read.
#'
#'@examples
#'
#'# To define a the default model, with 10 classes (9 HBD and 1 non-HBD class)
#'# and with pre-defined rates for HBD classes with a base of 2 (2, 4, 8, ...):
#'
#'mix10R <- zoomodel()
#'
#'# To see the parameters of the defined model, just type:
#'
#'mix10R
#'
#'# To define a model with pre-defined rates for 5 classes (4 HBD and 1 non-HBD
#'# class) and using a base of 10 to define rates (10, 100, 1000, ...):
#'
#'mix5R <- zoomodel(K=5,base=10)
#'
#'# To define a model with two classes, with estimation of rates for HBD classes
#'# and starting with a rate 10:
#'
#'my.mod1R <- zoomodel(predefined=FALSE,K=2,krates=c(10,10))
#'
#'# To define a model with four classes, with estimation of rates for HBD classes
#'# and choosing four initial rates:
#'
#'my.mod4R <- zoomodel(predefined=FALSE,K=4,krates=c(16,64,256,256))
#'
#'
#'@export

zoomodel <- function(predefined = TRUE, K = 10, mix_coef = rep(0, K), base_rate = 2,
                         krates = rep(0, K), err = 0.001, seqerr = 0.001, layers = TRUE) {
  mymod <- new("zmodel")
  mymod@typeModel = "kr"
  if(predefined){mymod@typeModel = "mixkr"}
  if(predefined & layers){mymod@typeModel = "mixkl"}
  if(!predefined & layers){mymod@typeModel = "kl"}

  if(sum(mix_coef) != 1){
    mymod@mix_coef <- c(rep(0.01,(K-1)),(1.01 - K*0.01))
  } else {
    if(length(mix_coef) != K){print("Length of mix_coeff does not match K !")}
    mymod@mix_coef = mix_coef
    }

  if(sum(krates == 0)){
    mymod@krates <- array(0,K)
    for (i in 1:(K-1)){mymod@krates[i] = base_rate**i}
    mymod@krates[K] = base_rate**(K-1)
  } else {
    mymod@krates = krates
    if(length(krates[krates < 0])){print("Problem: negative rates for HBD classes !!!")}
    if(length(krates) != K){print("Length of krates does not match K !")}
  }

  mymod@err = err
  mymod@seqerr = seqerr

  return(mymod)
}
