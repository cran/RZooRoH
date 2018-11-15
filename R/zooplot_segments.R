#'Plot HBD segments identified with the ZooROH model
#'
#'Plot HBD segments identified with the ZooRoH model for one or several
#'populations.
#'
#'@param input a named list with one or several zres objects obtained after
#'  running zoorun. The zres objects are the output of the zoorun function. For
#'  instance, putting list(name1 = zres1, name2 = zres2). The function will then
#'  use the names in the plot (in case several zres objects are used).
#'
#'@param chr the number of the chromosome where we are looking for HBD segments.
#'  This chromosome number refers to the position of the chromosome in the list
#'  of all chromsomes present in the input genotype data.
#'
#'@param coord a vector with the start and end position (in bp) of the region to
#'  plot.
#'
#'@param minlen the minimal length (in cM or Mb) of HBD segments to be plotted
#'  (set to 0 by default).
#'
#'@param cols a vector with the colors to be used for each population or zres
#'  object.
#'
#'@param plotids a logical indicating whether the IDs of the individuals are
#'  plotted on the graph (TRUE by default).
#'
#'@param toplot a list of vectors indicating the zres@@ids to be plotted. This
#'  option can be used to select the individuals to plot. The list must contain
#'  one vector per population or zres object. By default, all individuals are
#'  plotted.
#'
#'@param randomids a logical indicating whether a randomset of individuals is
#'  plotted. This option allows to reduce the number of individuals in the plot.
#'  The option can not be used simultaneously with the toplot option. By
#'  default, randomids is FALSE.
#'
#'@param nrandom a vector indicating the number of individuals to be randomly
#'  sampled per population or per zres object when randomids is TRUE. By
#'  default, we select 10 individuals per zres object. This vector must have the
#'  same length as the input list.
#'
#'@param seed a value for the random seed used to sample individuals to plot
#'  (when the randomids option is TRUE).
#'
#'@return The function plots the HBD segments identified in the region, using
#'  different colors for different zres object. Each line represents a different
#'  individual.
#'
#'@export
#'@import RColorBrewer

zooplot_hbdseg <- function (input,  chr = NULL,  coord = NULL, minlen = 0, cols = NULL,
                    plotids = TRUE, toplot = NULL, randomids = FALSE,
                    nrandom = (rep(10, length (input))), seed = 100){

  minlen = minlen * 1000000
  if (is.null (chr) | is.null (coord)) {
    stop ('Please provide the chromosome and coordinates\n')
  }

  if (is.null (cols)) {
    cols <- brewer.pal (8, "Dark2")
  }

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

  set.seed=seed
  ns <- c()
  yats <- c()
  popats <- c() #to plot the population labels
  allids <- c()
  ranids <- list()

  ##missing hbds
  for (i in 1:length (input)) {
    myres <- input[[i]]@hbdseg [input[[i]]@hbdseg[, 2]==chr, ]
    missing <- input[[i]]@ids [! input[[i]]@ids %in% unique (myres [, 1])]
    if (length (missing) >0) {
      toadd <- data.frame (missing, chr, 0, 0, 0, 0, 0, 0, 0)
      colnames (toadd) <- colnames(myres)
      myres <- rbind (myres, toadd)
    }


    if(randomids == FALSE & is.null(toplot)) {
     myres <- myres
     myids <- input[[i]]@sampleids
     allids <- c(allids, myids)
    }

    if ( !is.null (toplot)  &  randomids == TRUE) {
      stop ("Choose only one option, either toplot, either randomids.\n")
    }

    if (randomids==TRUE) {
      if (i ==1) {
        cat ("random seed ",  seed,  " is used to sample ids\n")
      }

      if(!is.null (nrandom)  & (length(nrandom) == length(input)   )  ){
        if (length (input[[i]]@ids) >= nrandom [i]) {
        myids <- sample (input[[i]]@ids, nrandom [i])
        ranids[[i]] <- myids
        allids <- c (allids,  input[[i]]@sampleids [myids])
          myres <- myres [myres [, 1] %in% myids, ]
        }else {
          stop ("nrandom for populaltion ",  i,  " larger than the number of samples\n")
        }
        }else {
          stop ("When randomids is TRUE, a VECTOR of length equal to the input list
                must be provided with the option nrandom.\n")
        }
      }

      if ( !is.null (toplot)  &  (length(toplot) != length(input))) {
        stop ("toplot must be a list with the same length as the input list.\n")
        }

      if ( !is.null (toplot)  &  (length(toplot)==length(input))) { #may be i should check that toplot is a list of length n?
          myids <- input [[i]]@ids [ input[[i]]@ids  %in% toplot [[i]] ]
        if (sum (myids %in% myres [, 1]   ) == length (myids)  ){
          myres <- myres [myres [, 1] %in% myids, ]  ##the list will be sorted [order not same as the input]
          myids <- toplot[[i]]
          allids <- c (allids,  input [[i]]@sampleids[myids])
        }else {
          stop ("Results not found for ids provided for data ",  i,  "\n")
        }
      }

    ns [i] <- length (unique (myres[, 1]))
    myres$pop <- i
    if (i==1) {
      allres <- myres
    }else {
      allres <- rbind (allres, myres)
    }
  }
  filres <- allres [allres$length >= minlen & allres[, 6]>=coord[1] & allres [, 5]<= coord[2] ,  ]
  filres [, 5] <-   filres [, 5] /1000000
  filres [, 6] <-   filres [, 6] /1000000

  ##here should i get the limits from the data ?
  xlim <- coord/1000000


  if (plotids==TRUE){
    par (mar =c(6, 6, 2, 4))
  }else {
    par (mar=c(6, 2, 2, 4))
  }
  if (length (input) >1) {
    par (oma =c(0, 6, 0, 0))
  }

  plot (1,  type="n", xlim=xlim,  ylim=c(0,  (sum(ns)*4) + (length (input)*4)   ), ylab="", yaxt="n", xlab=paste("Postion on Chromosome-", chr, " (Mb) ",    sep=""), xaxs="i", yaxs="i")
  yat <- 0
  for (i in 1:length (input)) {
    if(randomids==FALSE & is.null(toplot)){ids <- input[[i]]@ids}
    if(!is.null(toplot)){ids <- toplot[[i]]}
    if(randomids==TRUE){ids <- ranids[[i]]}
#    ids <- unique(allres [allres$pop ==i, 1]) ##take the ids from the fullres
    for (j in 1:length (ids)) {
      yat <- yat + 4
      myind <- filres [filres$pop==i & filres[, 1]==ids [j], ]
      lines (c(xlim[1], xlim[2]), c(yat, yat), col=cols[i])
      yats <- c(yats, yat)
      if (nrow (myind) >0) {
        for (k in 1:nrow (myind)) {
          rect (myind [k, 5],  yat-1,  myind[k, 6],  yat+1, col=cols[i], border=cols [i])
        }
      }
    }

    yat <- yat + 4 ##space after plotting individuals from one population [if not used chanage the ylim !]
    popats [i] <- yat
  }

  if (plotids==TRUE) {
    axis (side=2, at=yats,  labels = allids, tick=FALSE, las=2, cex.axis=0.7)
  }

  if (length (input)>1) {
    labats <-    0.5*(c(0,   popats [-length (popats)]) + popats)
    for(l in  1:length (input)) {
      axis (side=2, at=labats[l], labels=names(input)[l], outer=TRUE, tick=FALSE, cex.axis=1, col.axis=cols[l], font.axis=2, padj=0, las=3) #pop label colors ..font.axis=2 for bold
    }
  }
}
