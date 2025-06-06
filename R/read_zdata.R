
setClass(Class = "zdata",
         representation(genos = "matrix", bp = "vector", chrnames = "vector", chrbound = "matrix",
                        nind = "numeric", nsnps = "numeric", freqs ="vector", nchr ="numeric",
                        zformat ="character", sample_ids = "vector")
)

is.zdata <- function (x)
{
  res <- (is(x,"zdata") & validObject(x))
  return(res)
}

#### function to get allele frequencies from genotype data

getfreq1 <- function(genoline) {
  n1=length(genoline)
  n2=length(genoline[genoline<3])
  if((n2/n1) > 0.00){
    f0 <- sum(genoline[genoline<3]) / (2*length(genoline[genoline<3]))
  }else{
    f0 <- 0.00
  }
  return(f0)
}

getfreq2 <- function(genoline, genformat){
  if(genformat == 'gp' | genformat== 'gl'){
    i1 <- seq(from=1,to=(length(genoline)),by=3)
    i2 <- seq(from=2,to=(length(genoline)),by=3)
    i3 <- seq(from=3,to=(length(genoline)),by=3)
    g11 <- genoline[i1]
    g12 <- genoline[i2]
    g22 <- genoline[i3]
    g0 <- g11 + g12 + g22
    if(genformat == 'gp'){g0[g0 < 0.05]=0} #### considered missing
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
    g0[g11 >= 0.33 & g12 >= 0.33 & g22 >= 0.33]=0 #### when coding missing as 0.33 / 0.33 / 0.33
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

#### function to estimate the allele frequencies from haplotype data

getfreq3 <- function(genoline) {
  genoline <- genoline[which(!is.na(genoline))]
  n1=length(genoline)
  n2=length(genoline[genoline<2])
  if((n2/n1) > 0.00){
    f0 <- sum(genoline[genoline<2]) / (length(genoline[genoline<2]))
  }else{
    f0 <- 0.00
  }
  return(f0)
}

#### functions to estimate the allele frequencies with the EM algorithm

getfreqem1 <- function(genoline){
  n <- length(genoline)/3
  f_est <- .Fortran("freqem1",as.double(genoline),as.integer(n),as.double(0))
  return(f_est[3][[1]][1])
}

getfreqem2 <- function(genoline){
  n <- length(genoline)/3
  f_est <- .Fortran("freqem2",as.double(genoline),as.integer(n),as.double(0))
  return(f_est[3][[1]][1])
}

getfreqem3 <- function(genoline){
  n <- length(genoline)/2
  f_est <- .Fortran("freqem3",as.integer(genoline),as.integer(n),as.double(0))
  return(f_est[3][[1]][1])
}

#'Read the genotype data file
#'
#'Read a data file and convert it to the RZooRoH format in a 'zooin' object
#'required for further analysis.
#'
#'@param genofile The name of the input data file. Note that the model is
#'  designed for autosomes. Other chromosomes and additional filtering (e.g.
#'  call rate, missing, HWE, etc.) should be performed prior to run RZooRoH with
#'  tools such as PLINK or bcftools for instance. The model works on an ordered
#'  map and ignores SNPs with a null position.
#'
#'@param min_maf The minimum allele frequency to keep variants in the analysis
#'  (optional / set to 0.00 by default to keep all markers). Values such as 0.01
#'  allows to exclude monomorphic markers that are not informative and to reduce
#'  the size of the data and computational costs. There is no marker exclusion
#'  on call rate. However, we expect that data filtering is done prior to
#'  RZooRoH with tools such as PLINK or vcftools.
#'
#'@param zformat The code corresponding to the format of the data file ("gt" for
#'  genotypes, "gp" for genotype probabilities, "gl" for genotype likelihoods in
#'  Phred scores, "ad" for allelic depths). For all these formats, markers are
#'  ordered per rows and individuals per columns. Variants should be ordered by
#'  chromosome and position. By default, the format is inspired from the
#'  Oxford/GEN format, and the first five columns are chromosome identification
#'  (e.g, "1", "chr1"), the name of the marker, the position of the marker in
#'  base pairs or better in cM multiplied by 1,000,000 when genetic distances
#'  are known, the first marker allele and the second marker allele. Information
#'  per individual varies according to the format. With the "gt" format we have
#'  one column per individual with 0, 1 and 2 indicating the number of copies of
#'  the first allele (and 9 for missing). With the "gp" format we have three
#'  column per individual with the probabilities of genotype 11 (homozygous for
#'  the first allele), genotype 12 and genotype 22 (this corresponds to the
#'  oxford GEN format). Similarly, with the "gl" format, we have three column
#'  per individual with the likelihoods for genotypes 11, 12 and 22 in Phred
#'  scores. Finally, with the "ad" format, we expect two columns per individual:
#'  the number of reads for allele 1 and the number of reads for allele 2. For
#'  these three last formats, missing values must be indicated by setting all
#'  elements to 0. If one of the columns is non-null for one individual, the
#'  genotype will be considered non-missing. Note that the marker alleles
#'  specified in columns 4 and 5 are not used.
#'
#'  Conversion of a PLINK ped file or a VCF file to RZooRoH format can easily be
#'  performed using PLINK (version 1.9) or using bcftools.
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
#'  exclude chromosomes (this is an example and chromosome names must be adapted
#'  if necessary). The needed output data is then outfile.gen.
#'
#'  If some genotype probabilities are missing, with a value of "-nan", you must
#'  replace them with "0" (triple 0 is considered as missing). This can be done
#'  with this command:
#'
#'  sed -e 's/-nan/0/g' file.gen > newfile.gen
#'
#'  Note that some software recode missing as 0.333, 0.333 and 0.333. This was
#'  not previously considered as missing. This has no consequence of emission
#'  probabilities but impacts slightly estimation of allele frequencies (when
#'  the freqem option is false). Now, if the three genotype probabilities are
#'  above 0.33, the genotype is considered missing. Ideally, missing should be
#'  set to 0 0 0.
#'
#'  For applications related to identity-by-descent (IBD) estimation, the ZooRoH
#'  model is applied on two phased chromosomes or two haploid chromosomes. Two
#'  additional formats are therefore available to provide haplotypes (or haploid
#'  data). The first is "vcf" which refers to a phased VCF (for example, the
#'  output from Beagle5). In that case, the VCF must contain only the phased
#'  haplotype information (e.g 0|1, 0|0, 1|1, etc). The phased VCF can only be used
#'  with diploid individuals (after phasing). By default, we assume that
#'  the first and second column indicate chromosome and SNP position and that
#'  the individual haplotypes information starts in column 10. The second format
#'  is "haps" that is similar to the "GEN" format with five columns (chromosome
#'  identification (e.g, "1", "chr1"), the name of the marker, the position of
#'  the marker in base pairs or better in cM multiplied by 1,000,000 when
#'  genetic distances are known, the first marker allele and the second marker
#'  allele) followed by the haplotypes (two columns for diploid individuals and
#'  one in case of haploid chromosome). The alleles are coded as 0 and 1. The "haps"
#'  format is the only format that can be used with haploid data.
#'
#'  As for the genotype formats, it is possible to use tools such as bcftools to
#'  format file in the correct format.
#'
#'@param chrcol An optional argument that indicates the column number where the
#'  chromosome information is indicated (first column by default for all
#'  formats).
#'
#'@param poscol An optional argument that indicates the column number where the
#'  marker position is indicated (third column by default for all formats except
#'  phased vcf where it is the second column).
#'
#'@param supcol An optional argument that indicates the number of additional
#'  columns before the individuals genotypes or haplotypes (five columns are
#'  expected by default as described in the zformat argument description, except
#'  for phased VCF where this value is set to 9). Note that the function
#'  requires at least two information columns: the chromosome number and the
#'  marker position.
#'
#'@param haploid An optional argument that indicates whether we work on an
#'  haploid organism or chromosome (default = FALSE). This is only compatible
#'  with the 'IBD' option from the zoorun function (we don't estimate HBD in
#'  haploid data!). In the case you use haploid data, the only possible format
#'  is "haps".
#'
#'@param allelefreq A vector with allele frequencies for the first marker allele
#'  (optional). By default, the allele frequencies are estimated from the data.
#'  The option allows to skip this computation or to provide external allele
#'  frequencies estimated by another method or on another data set.
#'
#'@param freqem A logical indicating whether allele frequencies should be
#'  estimated with an EM algorithm. By default they are estimated with simpler
#'  approaches. The approach is ignored with the GT format. For high confidence
#'  genotypes (e.g., genotyping arrays, high-coverage sequencing data), it is
#'  not necessary to use this EM approach as genotypes are known. The approach
#'  is more relevant with low-fold sequencing for example, and more so with the
#'  PL or GP format (the approximation with the AD format being closer to the
#'  EM).
#'
#'@param samplefile A file with names of the samples (optional). It must match
#'  with the number of genotypes. If none is provided, the position in the
#'  genofile is used as ID.
#'
#'@return The function return a zooin object called containing the following
#'  elements: zooin@@genos a matrix with the genotypes, genotype probabilities
#'  or haplotypes, zooin@@bp an array with marker positions, zooin@@chrbound a
#'  matrix with the first and last marker number for each chromosome,
#'  zooin@@nind the number of individuals, zooin@@nsnps the number of markers
#'  conserved after filtering for minor allele frequency, zooin@@freqs an array
#'  with the marker allele frequencies, zooin@@nchr the number of chromosomes,
#'  zooin@@zformat the format of the data ("gt","gp","gl","ad","vcf","haps") and
#'  zooin@@sample_ids (the names of the samples).
#'
#' @examples
#'
#' # Get the name and location of example files
#'
#' myfile1 <- system.file("exdata","genoex.txt",package="RZooRoH")
#' myfile2 <- system.file("exdata","genosim.txt",package="RZooRoH")
#'
#' # Load your data with default format into a zooin object named "data1":
#'
#' data1 <- zoodata(myfile1)
#'
#' # Load the first data file with default format and filtering out markers with MAF < 0.02
#' # into a zooin object called "data1frq002":
#'
#' data1frq002 <- zoodata(myfile1, min_maf = 0.02)
#'
#' # Load the first data file with default format, with external allele frequencies
#' # (here a random set we create) and filtering out markers with MAF < 0.01:
#'
#' myrandomfreq <- runif(14831)
#' data1c <- zoodata(myfile1, allelefreq = myrandomfreq, min_maf = 0.01)
#'
#' # Load the second data file and indicate your own format (chromosome number in column 1,
#' # map position in column 2, 4 columns before genotypes) and filtering out markers with
#' # MAF < 0.01. The created zooin object is called "Sim5":
#'
#' Sim5 <- zoodata(myfile2, chrcol = 1, poscol =2, supcol = 4, min_maf = 0.01)
#'
#'@export
#'@import data.table

zoodata<-function(genofile, min_maf=0.00, zformat = "gt", chrcol = 1, poscol = 0, supcol = 0,
                  haploid = FALSE, allelefreq = NULL, freqem = FALSE, samplefile = NA){

  if(supcol==0){
    if(zformat=="gt" | zformat=="gp" | zformat=="gl" | zformat=="ad" | zformat=="haps"){supcol=5}
    else if(zformat=="vcf"){supcol=9}
  }

  if(poscol==0){
    if(zformat=="gt" | zformat=="gp" | zformat=="gl" | zformat=="ad" | zformat=="haps"){poscol=3}
    else if(zformat=="vcf"){poscol=2}
  }

  zooin <- new("zdata")
  max_maf <- (1 - min_maf)
  freqs <- allelefreq
  if(zformat == "vcf"){
    genfile <- fread(genofile,data.table=FALSE,header=TRUE)
  }else{
    genfile <- fread(genofile,data.table=FALSE)
  }

  if(zformat != "gt" && zformat != "gp" && zformat != "gl" && zformat != "ad" && zformat != "vcf" && zformat !="haps"){
    stop(paste("Unknown data format (zformat) ::",zformat,"\n",sep=""))
  }

  if(haploid & zformat != "haps"){
    stop("You must use the haps format with haploid data!")
  }

  print(c('Number of positions in original file ::',length(genfile[,poscol])))

  chr <- as.character(genfile[, chrcol])
  bp <- genfile[, poscol]

  if (zformat == "gt") {
    nind <- length(genfile[1,])-supcol
    genos <- as.matrix(genfile[, (supcol+1):(nind+supcol)])
    rm(genfile)
    gc()
    if(is.null(allelefreq)){freqs <- apply(genos, 1, getfreq1)}
  }

  if (zformat == "gp") {
    nind <- (length(genfile[1,])-supcol)/3
    genos <- as.matrix(genfile[, (supcol+1):(3*nind+supcol)])
    rm(genfile)
    gc()
    if(is.null(allelefreq) & !freqem){freqs <- apply(genos, 1, getfreq2, genformat = "gp")}
    if(is.null(allelefreq) & freqem){freqs <- apply(genos,1,getfreqem1)}
  }

  if (zformat == "gl") {
    nind <- (length(genfile[1,])-supcol)/3
    genos <- as.matrix(genfile[, (supcol+1):(3*nind+supcol)])
    rm(genfile)
    gc()
    if(is.null(allelefreq) & !freqem){freqs <- apply(genos, 1, getfreq2, genformat = "gl")}
    if(is.null(allelefreq) & freqem){freqs <- apply(genos,1,getfreqem2)}
  }

  if (zformat == "ad") {
    nind <- (length(genfile[1,])-supcol)/2
    genos <- as.matrix(genfile[, (supcol+1):(2*nind+supcol)])
    rm(genfile)
    gc()
    if(is.null(allelefreq) & !freqem){freqs <- apply(genos, 1, getfreq2, genformat ="ad")}
    if(is.null(allelefreq) & freqem){freqs <- apply(genos, 1, getfreqem3)}
  }

  if (zformat == "vcf") {
    #genfile[genfile=="."] <- NA
    pvcf <- as.data.frame(genfile[, (supcol+1):ncol(genfile)])
    pvcf[pvcf=="."] <- ".|."
    rm(genfile)
    haps1 <- lapply(pvcf, function(col) sapply(strsplit(as.character(col), "\\|"), "[[", 1))
    haps2 <- lapply(pvcf, function(col) sapply(strsplit(as.character(col), "\\|"), "[[", 2))
    haps1[haps1=="."] <- NA
    haps2[haps2=="."] <- NA
    rm(pvcf)
    hap1 <- as.data.frame(lapply(haps1, as.integer))
    hap2 <- as.data.frame(lapply(haps2, as.integer))
    nind <- length(hap1[1,])
    genos <- matrix(9,nrow(hap1),2*nind)
    genos[,seq(1,2*nind-1,by=2)]=as.matrix(hap1)
    genos[,seq(2,2*nind,by=2)]=as.matrix(hap2)
    rm(hap1,hap2)
    gc()
    if(is.null(allelefreq)){freqs <- apply(genos, 1, getfreq3)}
  }

  if (zformat == "haps") {
    nhap <- length(genfile[1,])-supcol
    if(haploid){nind <- nhap}
    if(!haploid){nind <- nind <- nhap/2}
    genfile[genfile=="."] <- NA
    hapin <- as.data.frame(lapply(genfile[, (supcol+1):(nhap+supcol)], as.integer))
    genos <- as.matrix(hapin)
    rm(genfile)
    gc()
    if(is.null(allelefreq)){freqs <- apply(genos, 1, getfreq3)}
  }

  gc()

  filterin <- (freqs >= min_maf & freqs <= max_maf & bp > 0)
  chr <- chr[filterin]
  bp <- bp[filterin]
  ff <- freqs[filterin]
  nchr <- length(unique(chr))

  ### solution avoiding loop
  chrbound <- matrix(0, nchr, 2)
  chrnames <- array("",nchr)
  chrbp2 <- c(chr[2:(length(chr))],"chr+")
  chrbound[,2] <- which(chr != chrbp2)
  chrnames <- chr[which(chr != chrbp2)]
  chrbound[1,1] <- 1
  if(nchr > 1){for (i in 2:nchr){chrbound[i,1] <- chrbound[(i-1),2] + 1 }}

  nsnps <- chrbound[nchr,2]

  zooin@nind <- nind
  zooin@genos <- genos[filterin,]
  rm(genos,freqs,filterin)
  gc()
  zooin@nsnps <- nsnps
  zooin@freqs <- ff
  zooin@chrnames <- chrnames
  zooin@chrbound <- chrbound
  zooin@nchr <- nchr
  zooin@bp <- bp
  zooin@zformat <- zformat
  if(!is.na(samplefile)){
    mysamples <- read.table(samplefile,header=FALSE)
    if(length(mysamples$V1) != nind){
      print("The number of sample IDs does not match the number of genotypes !")
    }
    zooin@sample_ids <- as.character(mysamples$V1)
  } else {zooin@sample_ids <- as.character(seq(1:nind))}
  print(c('Number of positions after MAF filtering  ::',zooin@nsnps))

  return(zooin)

}

