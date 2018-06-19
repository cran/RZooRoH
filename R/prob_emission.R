#### function to compute emission prob. (should vary according to the format)
#### gerr is optional (genotyping error rate)
#### gerr represents the probability to observe of heterozygote in HBD state
#### for format gt, gerr accounts for genotyping error and mutation
#### for format gp & gl, the sequencing error is normally accounted in the probabilities gp and gl, and gerr accounts for mutation
#### for format ad: gerr accounts for incorrect genotype calling with our approach (function of the depth) and gerr
#### for ad we could add a term in our formula to account for the probability to observe a read incompatible with the genotype

pemission <- function(id , gerr = 0.001, seqerr = 0.001, zformat = "gt"){

  pemission=matrix(1, .GlobalEnv$zooin@nsnps, 2)
  f1 <- .GlobalEnv$zooin@freqs
  f2 <- 1 - .GlobalEnv$zooin@freqs

  if(zformat == "gt"){ #### gt = number of allele 1, 0 for 22, 1 for 12 and 2 for 11

    pemission[.GlobalEnv$zooin@genos[,id]==0, 1]=f2[.GlobalEnv$zooin@genos[,id]==0]**2
    pemission[.GlobalEnv$zooin@genos[,id]==1, 1]=2*f1[.GlobalEnv$zooin@genos[,id]==1]*f2[.GlobalEnv$zooin@genos[,id]==1]
    pemission[.GlobalEnv$zooin@genos[,id]==2, 1]=f1[.GlobalEnv$zooin@genos[,id]==2]**2

    pemission[.GlobalEnv$zooin@genos[,id]==0, 2]=f2[.GlobalEnv$zooin@genos[,id]==0]*(1-gerr)
    pemission[.GlobalEnv$zooin@genos[,id]==1, 2]=gerr
    pemission[.GlobalEnv$zooin@genos[,id]==2, 2]=f1[.GlobalEnv$zooin@genos[,id]==2]*(1-gerr)

    return(pemission)
  }

  if(zformat == "gp"){ #### prob genotype 11, genotype 12 and genotype 22

    pr11=.GlobalEnv$zooin@genos[,(3*id-2)]
    pr12=.GlobalEnv$zooin@genos[,(3*id-1)]
    pr22=.GlobalEnv$zooin@genos[,(3*id)]

    pemission[,1]=pr11*f1**2+pr12*2*f1*f2+pr22*f2**2
    pemission[,2]=(pr11*f1+pr22*f2)*(1-gerr)+pr12*gerr

    pemission[((pr11+pr12+pr22) == 0),1]=1
    pemission[((pr11+pr12+pr22) == 0),2]=1

    return(pemission)
  }

  if(zformat == "gl"){ #### genotype likelihood for genotype 11, genotype 12 and genotype 22 in phred scores

    pr0 <- .GlobalEnv$zooin@genos[,(3*id-2)] + .GlobalEnv$zooin@genos[,(3*id-1)] + .GlobalEnv$zooin@genos[,(3*id)]

    pr11=10**(-.GlobalEnv$zooin@genos[,(3*id-2)]/10)
    pr12=10**(-.GlobalEnv$zooin@genos[,(3*id-1)]/10)
    pr22=10**(-.GlobalEnv$zooin@genos[,(3*id)]/10)

    prT=pr11+pr12+pr22
    pr11[prT > 0]=pr11[prT > 0]/prT[prT > 0]
    pr12[prT > 0]=pr12[prT > 0]/prT[prT > 0]
    pr22[prT > 0]=pr22[prT > 0]/prT[prT > 0]

    pemission[,1]=pr11*f1**2+pr12*2*f1*f2+pr22*f2**2
    pemission[,2]=(pr11*f1+pr22*f2)*(1-gerr)+pr12*gerr

    pemission[(pr0 == 0),1]=1
    pemission[(pr0 == 0),2]=1

    return(pemission)
  }

  if(zformat == "ad"){ ### AD, allele depth allele 1 / allele 2

    pr11=((1-seqerr)**.GlobalEnv$zooin@genos[,(2*id-1)])*(seqerr**.GlobalEnv$zooin@genos[,(2*id)])
    pr12=(0.5**.GlobalEnv$zooin@genos[,(2*id-1)])*(0.5**.GlobalEnv$zooin@genos[,(2*id)])
    pr22=(seqerr**.GlobalEnv$zooin@genos[,(2*id-1)])*((1-seqerr)**.GlobalEnv$zooin@genos[,(2*id)])

    prT=pr11+pr12+pr22
    pr11[prT > 0]=pr11[prT > 0]/prT[prT > 0]
    pr12[prT > 0]=pr12[prT > 0]/prT[prT > 0]
    pr22[prT > 0]=pr22[prT > 0]/prT[prT > 0]

    pemission[,1]=pr11*f1**2+pr12*2*f1*f2+pr22*f2**2
    pemission[,2]=(pr11*f1+pr22*f2)*(1-gerr)+pr12*gerr

    pemission[(.GlobalEnv$zooin@genos[,(2*id-1)]+.GlobalEnv$zooin@genos[,(2*id)]) == 0,1]=1
    pemission[(.GlobalEnv$zooin@genos[,(2*id-1)]+.GlobalEnv$zooin@genos[,(2*id)]) == 0,2]=1

    return(pemission)
  }

}

