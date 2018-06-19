#'Plot individual curves with proportion of the genome in each HBD class or cumulated
#'proportion in HBD classes with rates smaller than a threshold.
#'
#'For each individual, the function plots the mean percentage of the genome in different
#' HBD classes or the inbreeding coefficient obtained by summing autozygosity associated
#' with HBD classes with a rate lower or equal to a threshold (e.g., including all HBD
#'classes with longer and more recent HBD segments than a selected threshold).
#'
#'@param input a named list with one or several zres objects obtained after
#'  running zoorun. The zres objects are the output of the zoorun function. For
#'  instance, putting list(name1 = zres1, name2 = zres2). The function will then
#'  use the names in the plot (in case several zres objects are used).
#'
#'@param cumulative a logical indicating whether individual autozygosity is plotted
#'  per class (FALSE) or summed over all HBD class with a rate smaller than a
#'  value (these cumulated values are obtained for every rate defined in the
#'  model). By default, this value is TRUE. When FALSE, the percentages correspond
#'  to the individual genome-wide probabilities of belonging to each HBD-class
#'  or to the fraction of the genome in an autozygosity class. When TRUE, we obtain
#'  the probability of belonging to an HBD class with a rate smaller or equal than
#'  a threshold (here we use the pre-defined rates of the model as thresholds), averaged
#'  over the whole genome for each individual. This corresponds to report individual genomic
#'  inbreeding coefficients estimated with respect to different base populations obtained
#'  by selecting different thresholds T that determine which HBD classes are considered
#'  in the estimation of the genomic inbreeding coefficient (setting the base population
#'  approximately 0.5 * T generations ago).
#'
#'@param toplot A list of vectors indicating the zres@@ids to be plotted. This
#'  option can be used to select the individuals to plot. The list must contain
#'  one vector per population or zres object. By default, all individuals are
#'  plotted.
#'
#'@param ncols when several populations are plotted, ncols determines how many results (graphs)
#'are plotted per row.
#'
#'@return The function plots either the individual proportions of the genome associated with
#'different HBD classes or individual genomic inbreeding coefficients estimated with respect
#'to different base populations (from young to older). With both option, the average values are
#'plotted in red.
#'
#'@export

zooplot_individuals <- function (input,  cumulative=TRUE, toplot=NULL, ncols=2) {
  layout ( matrix (1:2, nrow=1), widths=c(0.90, 0.1))
  par (mar =c(1, 0, 4, 3))
  par (oma=c(6, 6, 0, 0))

  if (is (input, "list")) {
    if(any (lapply (input, class) != "zres")) {
      stop ("Some objects are NOT of class  \"zres\"\n")
    }
  }else {
    if(is (input,"zres")){
      input <- list(input)
    }
    else {
      stop ("input should be a list\n")
    }}

  if(length(names(input))==0 & length(input) > 1){
    warning("No names were provided for the input list!\n We will use capital letters.\n")
    names(input)=LETTERS[1:length(input)]}


  ks <- c()
  for (i in 1:length (input)) {
    myres <- input[[i]]@realized
    ks [i] <- ncol (myres)
    myres <- myres [, -c(ncol(myres))]
    myres <- data.frame (id=input[[i]]@ids, fullid=input[[i]]@sampleids, myres)


    myres$pop <- i
    if (i ==1) {
      allres <- myres
    }else {
      allres <- rbind(allres, myres)
    }
  }

  if (length (unique (ks)) >1) {
    stop  ("different models used for the data\n")
  }else {
    k <- unique (ks) -1
  }


  if (cumulative ==TRUE) {
    xlab=expression("Value of the threshold T used to estimate F"['G-T'] * " (using HBD classes with R"['K']<=" T)")
    ylab =expression ("Genomic inbreeding coefficient" ~ ( F [G-T]  ))
    ylim <- c(0, max(apply (allres [, 3:(k+2)], 1, cumsum))) * 1.06

  }else {
    xlab="Rate of the HBD class"
    ylab="Proportion of the genome in HBD class"
    ylim <- c(0,  max(allres [, 3:(k+2)])) *1.06
  }

  npop = length(input)
  if(npop>1 & npop%%2==1){npop=npop+1}
  mymat <- matrix (1:npop, ncol=ncols, byrow=TRUE)
  layout (mymat)

  for (j in 1:length (input)) {
    myres <- allres [allres$pop ==j, ]
    mymean <- apply (myres[, 3:(k+2)] , 2, mean)
    if (!is.null (toplot)) {
      if (length (toplot) != length (input)) {
        stop ('length of topot should match the length of input\n')
      }
      if (sum (myres$id %in% toplot [[j]]) == length (toplot [[j]])  ) {
        myres <- myres [myres$id %in% toplot [[j]], ]
        cat ('population ',  j,  ": ", nrow (myres),  "\n")
      }else {
        warning ("some ids for data ",  j,  " was not found\n")
      }
    }else {
      if (j ==1) {
        warning ("All individuals are plotted; use toplot to select individuals\n")
      }
    }


    ##mymean <- apply (myres [, 3:(2+k)] , 2, mean)
    plot (1,  type="n",xlim=c(1, k), ylim=ylim,  xaxt="n", cex.lab=2, cex.axis=1.5, bty="l", xlab="", ylab="", main=names(input)[j])
    axis (side=1, at=1:k, labels=input[[1]]@krates [1, 1:k], cex.axis=1.5, cex.lab=2)

    for (n in 1:nrow (myres)){
      if (cumulative==TRUE) {
        lines (1:k, cumsum(  as.numeric(myres[n,  3:(k+2)])  )  ,  type="l", col="gray")
        if (n == nrow (myres)) {
          lines (1:k, cumsum(mymean) ,  type="l", col="indianred2", lwd=2)
        }

      }else {
        lines (1:k, myres[n,  3:(k+2)]  ,  type="l", col="gray")
        if (n == nrow (myres)) {
          lines (1:k, mymean ,  type="l", col="indianred2", lwd=2)
        }
      }
    }
  }
  mtext(side=1, xlab, line=3, cex=1, outer=TRUE)
  mtext(side=2, ylab, line=3, cex=1, outer=TRUE)
}
