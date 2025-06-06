## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("RZooRoH")

## ----out.width="0.95\\linewidth", include=TRUE, fig.align="center", fig.cap=c("Distribution of length of HBD segments identified in Belgian Blue beef cattle with different models (results from Sole et al. (2017)). We observe that a model with 1 HBD class (1R) has more segments of intermediary size (10-100 kb). The model with 14 classes identifies more long segments (> 500 kb)."), echo=FALSE----
knitr::include_graphics("Figure1.pdf")

## ----out.width="0.95\\linewidth", include=TRUE, fig.align="center", fig.cap=c("Estimated level of autozygosity per HBD class in five humans (A), five dogs (B) and five sheeps (C). Each colour is associated with a distinct class (defined by its rate). The heights of each colour bar represent the estimated level of autozygosity associated with the class, and the total height represents the total estimated autozygosity (results from Druet and Gautier (2017)). Three dogs have total values close to 0.5 and the method indicates that 25 percent of the genome is associated with HBD classes with rates equal to 2 or 4. These values are compatible with parent-offspring matings combined with additional autozygosity from more distance ancestors. The most inbred sheep has autozygosity levels higher than 0.30 but mainly associated with more distant ancestors (approximately 8 generations in the past). These examples illustrate that the partitioning in different HBD class might help to understand the relationship between the parents or past demographic events at the individual level."), echo=FALSE----
knitr::include_graphics("Fig4hb_barplot_filteredRR.pdf")

## ----out.width="0.85\\linewidth", include=TRUE, fig.align="center", fig.cap=c("Comparison of individual autozygosity estimated with FEstim (1 HBD class) and a pruning strategy with recent autzoygosity classes estimated with a multiple HBD classes model (results from Druet and Gautier (2017))."), echo=FALSE----
knitr::include_graphics("Figure2.pdf")

## -----------------------------------------------------------------------------
file1 <- system.file("exdata","BBB_PE_gt_subset.txt",package="RZooRoH")
myfile1 <- read.table(file1,header=FALSE)
head(myfile1)

## -----------------------------------------------------------------------------
file2 <- system.file("exdata","BBB_NMP_pl_subset.txt",package="RZooRoH")
myfile2 <- read.table(file2,header=FALSE)
head(myfile2[,1:14])

## -----------------------------------------------------------------------------
file3 <- system.file("exdata","BBB_NMP_GP_subset.txt",package="RZooRoH")
myfile3 <- read.table(file3,header=FALSE)
head(myfile3[,1:14])

## -----------------------------------------------------------------------------
file4 <- system.file("exdata","BBB_NMP_ad_subset.txt",package="RZooRoH")
myfile4 <- read.table(file4,header=FALSE)
head(myfile4[,1:15])

## -----------------------------------------------------------------------------
file5 <- system.file("exdata","TAF_phased_ex.vcf",package="RZooRoH")
myfile5 <- read.table(file5,header=FALSE)
head(myfile5)

## -----------------------------------------------------------------------------
file6 <- system.file("exdata","Ref22_EX.haps",package="RZooRoH")
myfile6 <- read.table(file6,header=FALSE)
head(myfile6[,1:18])

## -----------------------------------------------------------------------------
library(RZooRoH)

## -----------------------------------------------------------------------------
file3 <- system.file("exdata","BBB_NMP_GP_subset.txt",package="RZooRoH")
BBB_GP <- zoodata(genofile = file3, zformat = "gp")

## -----------------------------------------------------------------------------
BBB_GP <- zoodata(file3, zformat = "gp")

## -----------------------------------------------------------------------------
file2 <- system.file("exdata","BBB_NMP_ad_subset.txt",package="RZooRoH")
BBB_AD <- zoodata(file2, zformat = "ad")

## -----------------------------------------------------------------------------
file1 <- system.file("exdata","BBB_PE_gt_subset.txt",package="RZooRoH")
BBB_GT <- zoodata(file1, supcol = 4, poscol = 2, chrcol = 1)

## -----------------------------------------------------------------------------
BBB_GT <- zoodata(file1, supcol = 4, poscol = 2, chrcol = 1, min_maf = 0.05)

## -----------------------------------------------------------------------------
BBB_AD2 <- zoodata(file2, zformat = "ad", freqem = TRUE)

## -----------------------------------------------------------------------------
mysamples <- system.file("exdata","BBB_samples.txt",package="RZooRoH")
BBB_GT <- zoodata(file1, supcol = 4, poscol = 2, chrcol = 1, samplefile = mysamples)

## -----------------------------------------------------------------------------
myfile6 <- system.file("exdata","Ref22_EX.haps",package="RZooRoH")
data6 <- zoodata(myfile6, zformat = "haps", haploid = TRUE)

## -----------------------------------------------------------------------------
head(BBB_GT@genos)

head(BBB_GP@genos[,1:6])

## -----------------------------------------------------------------------------
c(BBB_GT@nind, BBB_GT@nsnps, BBB_GT@nchr)

## -----------------------------------------------------------------------------
head(cbind(BBB_GT@bp, BBB_GT@freqs))

## -----------------------------------------------------------------------------
BBB_GT@sample_ids

## -----------------------------------------------------------------------------
MyMat = matrix(0,10,3)
MyMat[(1:2),1]=1;MyMat[(3:5),2]=1;MyMat[(6:10),3]=1
MyMat

## -----------------------------------------------------------------------------
StMix10L <- zoomodel(K=10,krates=c(2,4,6,8,10,12,14,16,18,20),mix_coef = c(0.01,0.01,0.01),
                     step = TRUE, XM = MyMat)

## -----------------------------------------------------------------------------
StMix10L@typeModel

## -----------------------------------------------------------------------------
Mix4I <- zoomodel(K=4,HBDclass ="Interval",krates=c(5,10,20,50),mix_coef=rep(0.001,4))

## -----------------------------------------------------------------------------
mix10L <- zoomodel()
mix10L

## -----------------------------------------------------------------------------
mix10L@krates

## -----------------------------------------------------------------------------
mix4L <- zoomodel(K=4,base=10)
mix4L

## -----------------------------------------------------------------------------
mix5R <- zoomodel(K=5,base=10,err=0.01,seqerr=0.005)
mix5R@err
mix5R@seqerr

## -----------------------------------------------------------------------------
my.mod3L <- zoomodel(predefined=FALSE,K=3,krates=c(16,64,256))
my.mod3L

## -----------------------------------------------------------------------------
my.mod3L <- zoomodel(predefined=FALSE,K=3,krates=c(16,64,256),
                     mix_coef=c(0.03,0.04,0.13))
my.mod3L@mix_coef

## -----------------------------------------------------------------------------
my.mod1L <- zoomodel(predefined=FALSE,K=1,krates=c(10))
my.mod1L

## -----------------------------------------------------------------------------
freqfile <- (system.file("exdata","typsfrq.txt",package="RZooRoH"))
typfile <- (system.file("exdata","typs.txt",package="RZooRoH"))
frq <- read.table(freqfile,header=FALSE)
bbb <- zoodata(typfile,supcol=4,chrcol=1,poscol=2,allelefreq=frq$V1)

## -----------------------------------------------------------------------------
bbb_results <- zoorun(my.mod1L, bbb, ids = c(1))

## -----------------------------------------------------------------------------
bbb_results1 <- zoorun(my.mod1L, bbb, ids = c(1,2), maxr = 100, minmix = 1e-8)

## -----------------------------------------------------------------------------
bbb_results2 <- zoorun(my.mod1L, bbb, ids = c(1,2), optim_method = "Nelder-Mead", 
                       maxiter = 1000)

## -----------------------------------------------------------------------------
bbb_results3 <- zoorun(my.mod1L, bbb, localhbd = TRUE)

## -----------------------------------------------------------------------------
Mod2L <- zoomodel(K=2,base_rate=10)
bbb_mod2l <- zoorun(Mod2L, bbb, localhbd = TRUE)

## -----------------------------------------------------------------------------
DTAF <- zoodata(file5,zformat="vcf")

## -----------------------------------------------------------------------------
TPairs <- as.matrix(cbind(c(3,3,1,2),c(1,2,2,2),c(4,4,2,4),c(1,1,1,1)))
TPairs

## -----------------------------------------------------------------------------
TIBD <- zoorun(Mod2L, DTAF, ibd = TRUE, ibdpairs = TPairs) 

## -----------------------------------------------------------------------------
D22 <- zoodata(file6,zformat="haps",haploid = TRUE)
Dpairs <- as.matrix(cbind(c(1,3,7,22,6),c(21,5,10,15,16)))
DIBD <- zoorun(my.mod1L,D22,ibd=TRUE,ibdpairs=Dpairs,haploid=TRUE)

## -----------------------------------------------------------------------------
bbb_results2@nind
bbb_results2@ids

## -----------------------------------------------------------------------------
TIBD@nind
TIBD@ids
TIBD@sampleids

## -----------------------------------------------------------------------------
bbb_mod2l@mixc

## -----------------------------------------------------------------------------
bbb_mod2l@krates

## -----------------------------------------------------------------------------
bbb_results3@krates

## -----------------------------------------------------------------------------
cbind(bbb_results3@modlik,bbb_results3@modbic)

## -----------------------------------------------------------------------------
bbb_mod2l@realized

## -----------------------------------------------------------------------------
bbb_results3@realized

## -----------------------------------------------------------------------------
cbind(bbb_results3@realized,bbb_results3@mixc)

## -----------------------------------------------------------------------------
mix10L <- zoomodel()
mix10L

## -----------------------------------------------------------------------------
round(soay_mix10l@realized[1:10,],3)

## ----out.width="80%",fig.align="center"---------------------------------------
x <- 1-soay_mix10l@realized[,11]
dpar <- par()
hist(x,nc=20,main="",xlab="Inbreeding coefficient",xlim=c(0.15,0.35),col='tomato')

## ----out.width="80%",fig.align="center"---------------------------------------
x <- t(apply(soay_mix10l@realized[,1:6],1,cumsum))
hist(x[,6],nc=20,main="",xlab="Inbreeding coefficient (T = 64)",
     xlim=c(0.15,0.35),col='tomato')

## -----------------------------------------------------------------------------
t(bbb_mod2l@hbdp[[3]][,1:15])

## -----------------------------------------------------------------------------
t(bbb_mod2l@hbdp[[6]][,4700:4720])

## -----------------------------------------------------------------------------
dim(bbb_mod2l@hbdseg)[1]
head(bbb_mod2l@hbdseg[,1:8])

## -----------------------------------------------------------------------------
summary(bbb_mod2l@hbdseg$length)
summary(bbb_mod2l@hbdseg$number_snp)

## -----------------------------------------------------------------------------
dim(soay_mix10l@hbdseg)[1]
head(soay_mix10l@hbdseg[,1:8])

## -----------------------------------------------------------------------------
summary(soay_mix10l@hbdseg$length)
summary(soay_mix10l@hbdseg$number_snp)

## -----------------------------------------------------------------------------
realized(bbb_mod2l)

## -----------------------------------------------------------------------------
head(round(realized(soay_mix10l,seq(1,6)),5))

## -----------------------------------------------------------------------------
F100 <- cumhbd(soay_mix10l, 100)
summary(F100)

## -----------------------------------------------------------------------------
F10 <- cumhbd(bbb_mod2l, 10)
F100 <- cumhbd(bbb_mod2l, 100)
cbind(F10,F100)

## -----------------------------------------------------------------------------
F20 <- cumhbd(bbb_results3, 20)
F50 <- cumhbd(bbb_results3, 50)
F100 <- cumhbd(bbb_results3, 100)
cbind(F20,F50,F100,bbb_results3@krates[,1])

## -----------------------------------------------------------------------------
roh25_10_20 <- rohbd(zres = soay_mix10l, chrom = 25, 
                     startPos= 10000000,endPos = 20000000, inside = FALSE)
dim(roh25_10_20)
head(roh25_10_20[,1:8])
summary(roh25_10_20$length)

## -----------------------------------------------------------------------------
roh25_10_20 <- rohbd(zres = soay_mix10l, chrom = 25, 
                     startPos= 10000000,endPos = 20000000, inside = TRUE)
dim(roh25_10_20)
head(roh25_10_20[,1:8])
summary(roh25_10_20$length)

## -----------------------------------------------------------------------------
roh25_10_20 <- rohbd(zres = soay_mix10l, ids = c(15,56,97,103,108))
dim(roh25_10_20)
head(roh25_10_20[,1:8])

## ----out.width="80%",fig.align="center"---------------------------------------
y6 <- probhbd(zres = bbb_mod2l, zooin = bbb, id = 6, chrom = 19, 
              startPos = 0, endPos = 50000000)
x <- bbb@bp[bbb@chrbound[19,1]:bbb@chrbound[19,2]]
x <- x[x >=0 & x<= 50000000]/1000000
plot(y6~x,type='b',ylim=c(0,1),ylab='HBD probability',col='red',
     xlab='Position on chr25 (Mb)')
y1 <- probhbd(zres = bbb_mod2l, zooin = bbb, id = 1, chrom = 19, 
              startPos = 0, endPos = 50000000)
y2 <- probhbd(zres = bbb_mod2l, zooin = bbb, id = 2, chrom = 19, 
              startPos = 0, endPos = 50000000)
par(new=TRUE)
plot(y1~x,type='b',ylim=c(0,1),ylab='',col='royalblue',xlab='',axes=FALSE)
par(new=TRUE)
plot(y2~x,type='b',ylim=c(0,1),ylab='',col='orange',xlab='',axes=FALSE)


## ----out.width="80%",fig.align="center"---------------------------------------
y6b <- probhbd(zres = bbb_mod2l, zooin = bbb, id = 6, chrom = 19, 
              startPos = 0, endPos = 50000000,T=20)
y1b <- probhbd(zres = bbb_mod2l, zooin = bbb, id = 1, chrom = 19, 
              startPos = 0, endPos = 50000000,T=20)
plot(y6b~x,type='l',ylim=c(0,1),ylab='HBD probability',col='red',
     xlab='Position on chr25 (Mb)')
par(new=TRUE)
plot(y1b~x,type='l',ylim=c(0,1),ylab='',col='royal blue',xlab='',axes=FALSE)

## -----------------------------------------------------------------------------
summary(realized(soay_mix10l))

## -----------------------------------------------------------------------------
summary(cumhbd(soay_mix10l,100))

## ----out.width="80%",fig.align="center"---------------------------------------
allseg <- rohbd(zres = soay_mix10l)
summary(allseg$length)
hist(allseg$length/1000000,xlab="Length of HBD segment (in cM)",main="",
     col='tomato',nc=100)

## ----fig.height=6,fig.width=10,out.width="80%",fig.align="center"-------------
zooplot_prophbd(list(Soay = soay_mix10l), cols = 'tomato', style = 'boxplot')

## ----fig.height=6,fig.width=10,out.width="80%",fig.align="center"-------------
zooplot_prophbd(list(Soay=soay_mix10l,Wiltshire=wilt_mix10l,
                     RasaAragonesa=rara_mix10l),style='barplot')

## ----fig.height=6,fig.width=12,out.width="80%",fig.align="center"-------------
zooplot_prophbd(list(Soay=soay_mix10l,Wiltshire=wilt_mix10l,
                     RasaAragonesa=rara_mix10l),style='lines', cumulative = TRUE)

## ----fig.height=6,fig.width=10,out.width="80%",fig.align="center"-------------
pop2 <- list(Soay=soay_mix10l,RasaAragonesa=rara_mix10l)
zooplot_individuals(pop2, cumulative = TRUE)

## ----fig.height=5,fig.width=10,show='hold',out.width="80%",fig.align="center"----
zooplot_partitioning(list(Wiltshire=wilt_mix10l), ylim = c(0,0.5), nonhbd = FALSE)

## ----fig.height=5,fig.width=10,show='hold'------------------------------------
pop3 <- list(Soay=soay_mix10l,RasaAragonesa=rara_mix10l,Wiltshire=wilt_mix10l)
zooplot_partitioning(pop3, randomids = TRUE, nrandom = c(20,20,20), plotids = FALSE,
                     ylim=c(0,0.5), nonhbd = FALSE)

## ----fig.height=10,fig.width=10,show='hold'-----------------------------------
zooplot_hbdseg(pop3, randomids = TRUE, nrandom = c(20,20,20),
               chr=5,coord=c(10000000,50000000))

## ----out.width="0.95\\linewidth", include=TRUE, fig.align="center", fig.cap=c("Impact on marker informativity on mean absolute error (MAE) for local HBD probabilities for all locus (left column) or for HBD locus only (right column) based on results from results from Druet and Gautier (2017). The figure shows the impact of average number of SNPs per HBD segment, site frequency spectrum (SFS), presence of genotyping errors or coverage with whole-genome sequencing (WGS) data."), echo=FALSE----
knitr::include_graphics("Effect_density_error_lowfoldII.pdf")

## ----out.width="0.72\\linewidth", include=TRUE, fig.align="center", fig.cap=c("HBD segments identified on chromosome 1 using different genotyping densities or using whole-genome sequencing at different coverages in 46 Belgian blue sires."), echo=FALSE----
knitr::include_graphics("BBB_hdbseg_densities.pdf")

## -----------------------------------------------------------------------------
targets <- as.matrix(cbind(c(1,1,1,2,2,3),c(2,3,4,3,4,4)))
TKIN <- zookin(Mod2L, DTAF, kinpairs = targets)

## ----eval=FALSE---------------------------------------------------------------
# taf3 <- zoodata("taf_phased.vcf",zformat="vcf", min_maf = 0.001)
# mix4L <- zoomodel(K=4,base_rate=5)
# kintaf_mix4l <- zookin(mix4L,taf3,kinpairs=mykinpairs,vit=FALSE,nT=4)

## ----out.width="75%",fig.align="center"---------------------------------------
par(mar=c(5.1,4.1,2.1,2.1),mfrow=c(1,1),mfcol=c(1,1))
barplot(t(kintaf_mix4l@realized),col=c('#003366','lightblue1',"lightgrey","whitesmoke"),
        xlab="Pairs of individuals",ylab="Kinship",cex.lab=1.2)

## ----eval=FALSE---------------------------------------------------------------
# rel <- kintaf_mix4L@realized[,1] > 0.05
# kintaf_local <- zookin(mix4l,taf3,kinpairs=mykinpairs[rel,],localhbd=TRUE)

## ----eval=FALSE---------------------------------------------------------------
# hbd1 <- predhbd(kintaf_local,zooin=taf3,chrom=1,num=1,startPos = 1,endPos = 140000000)
# map1 <- taf3@bp[1:length(hbd1)]/1000000
# hbd2 <- predhbd(kintaf_local,zooin=taf3,chrom=1,num=1,startPos=1,endPos=140000000,T=10)

## ----out.width="75%",fig.align="center"---------------------------------------
plot(hbd1~map1,type='l',ylim=c(0,1),ylab="IBD probability",
     xlab="Position on Chromosome 1")
par(new=T)
plot(hbd2~map1,type='l',ylim=c(0,1),ylab="",xlab="",axes=FALSE,col="brown1")

