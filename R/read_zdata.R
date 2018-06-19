
setClass(Class = "zdata",
         representation(genos = "matrix", bp = "vector", chrbound = "matrix", nind = "numeric",
                        nsnps = "numeric", freqs ="vector", nchr ="numeric", zformat ="character",
                        sample_ids = "vector")
)

is.zdata <- function (x)
{
  res <- (is(x,"zdata") & validObject(x))
  return(res)
}

#### function to get allele frequencies from genotype data

allelefreq <- function(genoline) {
  n1=length(genoline)
  n2=length(genoline[genoline<3])
  if((n2/n1) > 0.00){
    f0 <- sum(genoline[genoline<3]) / (2*length(genoline[genoline<3]))
  }else{
    f0 <- 0.00
  }
  return(f0)
}

getfreq <- function(genoline, genformat){
  if(genformat == 'gp' | genformat== 'gl'){
    i1 <- seq(from=1,to=(length(genoline)),by=3)
    i2 <- seq(from=2,to=(length(genoline)),by=3)
    i3 <- seq(from=3,to=(length(genoline)),by=3)
    g11 <- genoline[i1]
    g12 <- genoline[i2]
    g22 <- genoline[i3]
    g0 <- g11 + g12 + g22
    if(genformat == 'gl'){
      g11 <- 10**(-g11/10)
      g12 <- 10**(-g12/10)
      g22 <- 10**(-g22/10)
     }
    gT <- g11 + g12 + g22
    if(genformat == 'gl'){
      g11[g0 > 0] <- g11[g0 > 0]/gT[g0 > 0]
      g12[g0 > 0] <- g12[g0 > 0]/gT[g0 > 0]
      g22[g0 > 0] <- g22[g0 > 0]/gT[g0 > 0]
    }
    n0 <- length(g0[g0 >0])
    f0 <- sum(2*g11[g0 > 0] + 1*g12[g0 > 0])
    if(n0 > 0)f0 <- f0/(2*n0)
    if(n0 == 0)f0 <- 0.00
  }

  if(genformat == 'ad'){
    i1 <- seq(from=1,to=(length(genoline)),by=2)
    i2 <- seq(from=2,to=(length(genoline)),by=2)
    ad1 <- genoline[i1]
    ad2 <- genoline[i2]
    ad <- ad1 + ad2
    f0 <- sum(ad1[ad > 0]/ad[ad > 0])
    n0 <- length(ad[ad >0])
    if(n0 > 0)f0 <- f0/(n0)
    if(n0 == 0)f0 <- 0.00
  }
  return(f0)
}

#'Read the genotype data file
#'
#'Read a data file and convert it to the RZooRoH format in an object called
#'"zooin".
#'
#'@param genofile The name of the input data file. Note that the model is
#'  designed for autosomes. Other chromosomes and additional filtering (e.g.
#'  call rate, missing, HWE, etc.) should be performed prior to run RZooRoH with
#'  tools such as plink or bcftools for instance. The model works on an ordered
#'  map and ignores SNPs with a null position.
#'
#'@param min_maf The minimum allele frequency to keep variants in the analysis
#'  (optional / set to 0.00 by default to keep all markers). Values such as 0.01
#'  allows to exclude monomorphic markers that are not informative and to reduce
#'  the size of the data and computational costs. There is no marker exclusion
#'  on call rate. However, we expect that data filtering is done prior to RZooRoH
#'  with tools such as plink or vcftools.
#'
#'@param zformat The code corresponding to the format of the data file ("gt" for
#'  genotypes, "gp" for genotype probabilities, "gl" for genotype likelihoods in
#'  Phred scores, "ad" for allelic depths). For all these formats, markers are
#'  ordered per rows and individuals per columns. Variants should be ordered by
#'  chromosome and position. By default, the first five columns are chromosome
#'  identification (e.g, "1", "chr1"), the name of the marker, the position of
#'  the marker in base pairs or better in cM divided 1000000 when genetic
#'  distances are known, the first marker allele and the second marker allele.
#'  Information per individual varies according to the format. With the "gt"
#'  format we have one column per individual with 0, 1 and 2 indicating the
#'  number of copies of the first allele (and 9 for missing). With the "gp" format we
#'  have three column per individual with the probabilities of genotype 11
#'  (homozygous for the first allele), genotype 12 and genotype 22 (this
#'  corresponds to the oxford gen format). Similarly, with the "gl" format, we
#'  have three column per individual with the likelihoods for genotypes 11, 12
#'  and 22 in Phred scores. Finally, with the "ad" format, we expect two columns
#'  per individuals: the number of reads for allele 1 and the number of reads
#'  for allele 2. For these three last formats, missing values must be indicated
#'  by setting all elements to 0. If one of the columns is non-null for one
#'  individual, the genotype will be considered non-missing. Note that the
#'  marker alleles specified in columns 4 and 5 are not used.
#'
#'  Conversion of a plink ped file or a VCF file to RZooRoH format can easily be
#'  performed using plink (version 1.9) or using bcftools.
#'
#'  For ped files, recode them to oxford gen format with plink --file myinput
#'  --recode oxford --autosome --out myoutput. The autosome option keeps only
#'  SNPs on autosomes as required by RZooRoH.
#'
#'  For vcf files, bcftools can be used to recode a vcf to the oxford gen format
#'  with the convert option: bcftools convert -t ^chrX,chrY,chrM -g outfile
#'  --chrom --tag GT myfile.vcf. The --chrom option is important to obtain
#'  chromosome number in the first column. The tag option allows to select which
#'  field from the vcf file (GT, PL, GL or GP) is used to generate the genotype
#'  probabilities exported in the oxford gen format. The -t option allows to
#'  exclude chromosomes (this is example and chromosome names must be adapted if
#'  necessary). The needed output data is then outfile.gen.
#'
#'  If some genotype probabilities are missing, with a value of "-nan", you must
#'  replace them with "0" (triple 0 is considered as missing). This can be done
#'  with this command:
#'
#'  sed -e 's/-nan/0/g' file.gen > newfile.gen
#'
#'@param supcol An optional argument that indicates the number of additional
#'  columns before the individual genotypes (five columns are expected by
#'  default as described in the zformat argument description). Note that the
#'  function requires at least two information columns: the chromosome number
#'  and the marker position.
#'
#'@param chrcol An optional argument that indicates the column number where the
#'  chromosome information is indicated (first column by default).
#'
#'@param poscol An optional argument that indicates the column number where the
#'  marker position is indicated (third column by default).
#'
#'@param allelefreq A vector with allele frequencies for the first marker allele
#'  (optional). By default, the allele frequencies are estimated from the data.
#'  The option allows to skip this computation or to provide external allele
#'  frequencies estimated by another method or on another data set.
#'
#'@param samplefile A file with names of the samples. It must match with the
#' number of genotypes. If none is provided, the position in the genofile is
#' used as ID.
#'
#'@return The function return a global object called zooin containing the
#'  following elements: zooin@@genos a matrix with the genotypes or genotype
#'  probabilities, zooin@@bp an array with marker positions, zooin@@chrbound a
#'  matrix with the first and last marker number for each chromosome,
#'  zooin@@nind the number of individuals, zooin@@nsnps the number of markers
#'  conserved after filtering for minor allele frequency, zooin@@freqs an array
#'  with the marker allele frequencies, zooin@@nchr the number of chromosomes,
#'  zooin@@zformat the format of the data ("gt","gp","gl","ad") and
#'  zooin@@sample_ids (the names of the samples).
#'
#' @examples
#'
#' # Get the name and location of example files
#'
#' myfile1 <- system.file("exdata","genoex.txt",package="RZooRoH")
#' myfile2 <- system.file("exdata","genosim.txt",package="RZooRoH")
#'
#' # Load your data with default format:
#'
#' zoodata(myfile1)
#'
#' # Load the first data file with default format and filtering out markers with MAF < 0.02:
#'
#' zoodata(myfile1, min_maf = 0.02)
#'
#' # Load the first data file with default format, with external allele frequencies (here a
#' # random set) and filtering out markers with MAF < 0.01:
#'
#' myrandomfreq <- runif(14831)
#' zoodata(myfile1, allelefreq = myrandomfreq, min_maf = 0.01)
#'
#' # Load the second data file and indicate your own format (chromosome number in column 1,
#' # map position in column 2, 4 columns before genotypes) and filtering out markers with
#' # MAF < 0.01:
#'
#' zoodata(myfile2, chrcol = 1, poscol =2, supcol = 4, min_maf = 0.01)
#'
#'@export
#'@import data.table

zoodata<-function(genofile, min_maf=0.00, zformat = "gt", supcol = 5,
                      chrcol = 1, poscol = 3, allelefreq = NULL, samplefile = NA){

  zooin=NULL;zooin <<- new("zdata")
  max_maf <- (1 - min_maf)
  freqs <- allelefreq
  genfile <- fread(genofile,data.table=FALSE)

  if(zformat != "gt" && zformat != "gp" && zformat != "gl" && zformat != "ad"){
    stop(paste("Unknown data format (zformat) ::",zformat,"\n",sep=""))
  }

  print(c('Number of positions in original file ::',length(genfile[,poscol])))
  freqs <- freqs[genfile[,poscol]>0]
  genfile <- genfile[genfile[,poscol]>0,]
  print(c('Number of positions after removing null positions ::',length(genfile[,poscol])))

  if (zformat == "gt") {
    nind <- length(genfile[1,])-supcol
    genos <- as.matrix(genfile[, (supcol+1):(nind+supcol)])
    if(is.null(allelefreq)){freqs <- apply(genos, 1, allelefreq)}
  }

  if (zformat == "gp") {
    nind <- (length(genfile[1,])-supcol)/3
    genos <- as.matrix(genfile[, (supcol+1):(3*nind+supcol)])
    if(is.null(allelefreq)){freqs <- apply(genos, 1, getfreq, genformat = "gp")}
  }

  if (zformat == "gl") {
    nind <- (length(genfile[1,])-supcol)/3
    genos <- as.matrix(genfile[, (supcol+1):(3*nind+supcol)])
    if(is.null(allelefreq)){freqs <- apply(genos, 1, getfreq, genformat = "gl")}
  }

  if (zformat == "ad") {
    nind <- (length(genfile[1,])-supcol)/2
    genos <- as.matrix(genfile[, (supcol+1):(2*nind+supcol)])
    if(is.null(allelefreq)){freqs <- apply(genos, 1, getfreq, genformat ="ad")}
  }

  chr <- as.character(genfile[(freqs >= min_maf & freqs <= max_maf), chrcol])
  bp <- genfile[(freqs >= min_maf & freqs <= max_maf), poscol]
  ff <- freqs[(freqs >= min_maf & freqs <= max_maf)]
  nchr <- length(unique(chr))

  ### solution avoiding loop
  chrbound <- matrix(0, nchr, 2)
  chrbp2 <- c(chr[2:(length(chr))],"chr+")
  chrbound[,2] <- which(chr != chrbp2)
  chrbound[1,1] <- 1
  if(nchr > 1){for (i in 2:nchr){chrbound[i,1] <- chrbound[(i-1),2] + 1 }}

  nsnps <- chrbound[nchr,2]

  .GlobalEnv$zooin@nind <- nind
  .GlobalEnv$zooin@genos <- genos[(freqs >= min_maf & freqs <= max_maf),]
  .GlobalEnv$zooin@nsnps <- nsnps
  .GlobalEnv$zooin@freqs <- ff
  .GlobalEnv$zooin@chrbound <- chrbound
  .GlobalEnv$zooin@nchr <- nchr
  .GlobalEnv$zooin@bp <- bp
  .GlobalEnv$zooin@zformat <- zformat
  if(!is.na(samplefile)){
    mysamples <- read.table(samplefile,header=FALSE)
    if(length(mysamples$V1) != nind){
      print("The number of sample IDs does not match the number of genotypes !")
    }
    .GlobalEnv$zooin@sample_ids <- as.character(mysamples$V1)
  } else {.GlobalEnv$zooin@sample_ids <- as.character(seq(1:nind))}
  print(c('Number of positions after MAF filtering  ::',.GlobalEnv$zooin@nsnps))

}

