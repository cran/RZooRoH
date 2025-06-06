#### function to compute emission prob. (should vary according to the format)
#### gerr is optional (genotyping error rate)
#### gerr represents the probability to observe of heterozygote in HBD state
#### for format gt, gerr accounts for genotyping error and mutation
#### for format gp & gl, the sequencing error is normally accounted in the probabilities gp and gl, and gerr accounts for mutation
#### for format ad: gerr accounts for incorrect genotype calling with our approach (function of the depth) and gerr
#### for ad we could add a term in our formula to account for the probability to observe a read incompatible with the genotype

pemission <- function(zooin, id , gerr = 0.001, seqerr = 0.001, FH = 0.000, zformat = "gt", ibd, ibdpair, haploid){

  pemission=matrix(1, zooin@nsnps, 2)
  f1 <- zooin@freqs
  f2 <- 1 - zooin@freqs

  if(ibd){
    if(zformat == "vcf" | zformat == "haps"){ #### gt = hap1 + hap2 of allele 1, 0 for 22, 1 for 12 and 2 for 11

    if(!haploid){
      if(ibdpair[2]==1){i1=2*ibdpair[1]-1}
      if(ibdpair[2]==2){i1=2*ibdpair[1]}
      if(ibdpair[4]==1){i2=2*ibdpair[3]-1}
      if(ibdpair[4]==2){i2=2*ibdpair[3]}
    }else{
      i1=ibdpair[1]
      i2=ibdpair[2]
    }

      geno1 = zooin@genos[,i1]+zooin@genos[,i2]
      geno1[is.na(geno1)]=9

      pemission[geno1==0, 1]=(1-FH)*(f2[geno1==0]**2)+FH*(f2[geno1==0])
      pemission[geno1==1, 1]=(1-FH)*(2*f1[geno1==1]*f2[geno1==1])
      pemission[geno1==2, 1]=(1-FH)*(f1[geno1==2]**2)+FH*(f1[geno1==2])

      pemission[geno1==0, 2]=f2[geno1==0]*(1-gerr)
      pemission[geno1==1, 2]=gerr
      pemission[geno1==2, 2]=f1[geno1==2]*(1-gerr)

      rm(geno1)

      return(pemission)
    }

  }else{

  if(zformat == "gt"){ #### gt = number of allele 1, 0 for 22, 1 for 12 and 2 for 11

    pemission[zooin@genos[,id]==0, 1]=(1-FH)*(f2[zooin@genos[,id]==0]**2)+FH*(f2[zooin@genos[,id]==0])
    pemission[zooin@genos[,id]==1, 1]=(1-FH)*(2*f1[zooin@genos[,id]==1]*f2[zooin@genos[,id]==1])
    pemission[zooin@genos[,id]==2, 1]=(1-FH)*(f1[zooin@genos[,id]==2]**2)+FH*(f1[zooin@genos[,id]==2])

    pemission[zooin@genos[,id]==0, 2]=f2[zooin@genos[,id]==0]*(1-gerr)
    pemission[zooin@genos[,id]==1, 2]=gerr
    pemission[zooin@genos[,id]==2, 2]=f1[zooin@genos[,id]==2]*(1-gerr)

    return(pemission)
  }

  if(zformat == "gp"){ #### prob genotype 11, genotype 12 and genotype 22

    pr11=zooin@genos[,(3*id-2)]
    pr12=zooin@genos[,(3*id-1)]
    pr22=zooin@genos[,(3*id)]

    pemission[,1]=(1-FH)*(pr11*f1**2+pr12*2*f1*f2+pr22*f2**2)+FH*(pr11*f1+pr22*f2)
    pemission[,2]=(pr11*f1+pr22*f2)*(1-gerr)+pr12*gerr

    pemission[((pr11+pr12+pr22) == 0),1]=1
    pemission[((pr11+pr12+pr22) == 0),2]=1

    return(pemission)
  }

  if(zformat == "gl"){ #### genotype likelihood for genotype 11, genotype 12 and genotype 22 in phred scores

    pr0 <- zooin@genos[,(3*id-2)] + zooin@genos[,(3*id-1)] + zooin@genos[,(3*id)]

    pr11=10**(-zooin@genos[,(3*id-2)]/10)
    pr12=10**(-zooin@genos[,(3*id-1)]/10)
    pr22=10**(-zooin@genos[,(3*id)]/10)

    prT=pr11+pr12+pr22
    pr11[prT > 0]=pr11[prT > 0]/prT[prT > 0]
    pr12[prT > 0]=pr12[prT > 0]/prT[prT > 0]
    pr22[prT > 0]=pr22[prT > 0]/prT[prT > 0]

    pemission[,1]=(1-FH)*(pr11*f1**2+pr12*2*f1*f2+pr22*f2**2)+FH*(pr11*f1+pr22*f2)
    pemission[,2]=(pr11*f1+pr22*f2)*(1-gerr)+pr12*gerr

    pemission[(pr0 == 0),1]=1
    pemission[(pr0 == 0),2]=1

    return(pemission)
  }

  if(zformat == "ad"){ ### AD, allele depth allele 1 / allele 2

    pr11=((1-seqerr)**zooin@genos[,(2*id-1)])*(seqerr**zooin@genos[,(2*id)])
    pr12=(0.5**zooin@genos[,(2*id-1)])*(0.5**zooin@genos[,(2*id)])
    pr22=(seqerr**zooin@genos[,(2*id-1)])*((1-seqerr)**zooin@genos[,(2*id)])

    prT=pr11+pr12+pr22
    pr11[prT > 0]=pr11[prT > 0]/prT[prT > 0]
    pr12[prT > 0]=pr12[prT > 0]/prT[prT > 0]
    pr22[prT > 0]=pr22[prT > 0]/prT[prT > 0]

    pemission[,1]=(1-FH)*(pr11*f1**2+pr12*2*f1*f2+pr22*f2**2)+FH*(pr11*f1+pr22*f2)
    pemission[,2]=(pr11*f1+pr22*f2)*(1-gerr)+pr12*gerr

    pemission[(zooin@genos[,(2*id-1)]+zooin@genos[,(2*id)]) == 0,1]=1
    pemission[(zooin@genos[,(2*id-1)]+zooin@genos[,(2*id)]) == 0,2]=1

    return(pemission)
  }

  if(zformat == "vcf" | zformat == "haps"){ #### gt = hap1 + hap2 of allele 1, 0 for 22, 1 for 12 and 2 for 11

    geno1 = zooin@genos[,(2*id-1)]+zooin@genos[,(2*id)]
    geno1[is.na(geno1)]=9

    pemission[geno1==0, 1]=(1-FH)*(f2[geno1==0]**2)+FH*(f2[geno1==0])
    pemission[geno1==1, 1]=(1-FH)*(2*f1[geno1==1]*f2[geno1==1])
    pemission[geno1==2, 1]=(1-FH)*(f1[geno1==2]**2)+FH*(f1[geno1==2])

    pemission[geno1==0, 2]=f2[geno1==0]*(1-gerr)
    pemission[geno1==1, 2]=gerr
    pemission[geno1==2, 2]=f1[geno1==2]*(1-gerr)

    rm(geno1)

    return(pemission)
  }
  } ## if not ibd

}

