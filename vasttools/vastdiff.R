opt = list()

dat = read.table("test", header=1, sep="\t")

#-replicatesA=name1@name2@name3 -replicatesB=name4@name5
firstRepSet <- c(7,9,11,13,15,19,31,33)
secondRepSet <- c(17,21,23,25,27,29,35)

# Set number of replicates
firstRepN <- length(firstRepSet)
secondRepN <- length(secondRepSet)

# Make sure there are sample names
if(is.null( opt$sampleNameA ) ) {
  opt$sampleNameA <- firstRepSet[1]
}
if(is.null( opt$sampleNameB ) ) {
  opt$sampleNameB <- secondRepSet[1]
}
# Set output sample names for plot
sampOneName <- substr(opt$sampleNameA, 1, 9)
sampTwoName <- substr(opt$sampleNameB, 1, 9)


## INITIALIZE LISTS ##
shapeFirst <- vector("list", firstRepN)
shapeSecond <- vector("list", secondRepN)

psiFirst <- vector("list", firstRepN)
psiSecond <- vector("list", secondRepN)


# check if header is correct..  TODO

# Indexes of samples of interest
repAind <-firstRepSet
repBind <-secondRepSet


# Indexes of Quals
repA.qualInd <- repAind + 1
repB.qualInd <- repBind + 1

# make sure this succeeded  TODO

# CONST
alphaList <- seq(0,1,0.01)


### BEGIN READ INPUT ###
source("~/repos/vast-tools/R/Rlib/include_diff.R")
source("~/repos/vast-tools/R/Rlib/include.R")
opt$alpha = 1
opt$beta = 1
# Iterate through input, 'nLines' at a time to reduce overhead/memory

      tabLine <- dat[1,]
	 #writeLines(paste(tabLine[repA.qualInd], collapse="\t"), stderr());

      # Posterior parameters... Prior given from command line --alpha, --beta
      shapeFirst <- lapply( tabLine[repA.qualInd], function(x) {
								parseQual(x, opt$alpha, opt$beta)
							   } )
      shapeSecond <- lapply( tabLine[repB.qualInd], function(x) {
								parseQual(x, opt$alpha, opt$beta)
							    } )

      totalFirst <- unlist(lapply( shapeFirst, function(x) { x[1] + x[2] }))
      totalSecond <- unlist(lapply( shapeSecond, function(x) { x[1] + x[2] }))

      firstShapeMat <- do.call(rbind, shapeFirst)
      secondShapeMat <- do.call(rbind, shapeSecond)

      firstShapeAve <- c( mean(firstShapeMat[,1]), mean(firstShapeMat[,2]) )
      secondShapeAve <- c( mean(secondShapeMat[,1]), mean(secondShapeMat[,2]) )

      # Sample Posterior Distributions
      opt$minReads = 10
      opt$size = 500
      psiFirst <- lapply( shapeFirst, function(x) {
        #sample here from rbeta(N, alpha, beta) if > -e
        if(x[1]+x[2] < opt$minReads) { return(NULL) }
        rbeta(opt$size, shape1=x[1], shape2=x[2])
      })

      psiSecond <- lapply( shapeSecond, function(x) {
        #sample here from rbeta(N, alpha, beta)
        if(x[1]+x[2] < opt$minReads) { return(NULL) }
        rbeta(opt$size, shape1=x[1], shape2=x[2])
      })

      # calculate expected value of psi for each replicate
      expFirst <- unlist(lapply(shapeFirst, function(x) {
         if(x[1]+x[2] < opt$minReads) { return(NULL) }
         x[1] / (x[1] + x[2])
      }))
      expSecond <- unlist(lapply(shapeSecond, function(x) {
         if(x[1]+x[2] < opt$minReads) { return(NULL) }
         x[1] / (x[1] + x[2])
      }))


      # Create non-parametric Joint Distributions
      psiFirstComb <- do.call(c, psiFirst)
      psiSecondComb <- do.call(c, psiSecond)


      #    print(length(psiFirstComb))

      # if they aren't paired, then shuffle the joint distributions...
      library("psiplot")
      loadPackages(c("getopt", "optparse", "RColorBrewer", "reshape2", "ggplot2", "grid", "parallel", "devtools","MASS"), local.lib=paste(c(scriptPath,"/R/Rlib"), collapse=""))

        paramFirst <- try (suppressWarnings(
				fitdistr(psiFirstComb,
					"beta",
					list( shape1=firstShapeAve[1], shape2=firstShapeAve[2])
				)$estimate ), TRUE )
        paramSecond <- try (suppressWarnings(
				fitdistr(psiSecondComb,
					"beta",
					list( shape1=secondShapeAve[1], shape2=secondShapeAve[2])
				)$estimate ), TRUE )
        # if optimization fails its because the distribution is too narrow
        # in which case our starting shapes should already be good enough
	if(class(paramFirst) != "try-error") {
          psiFirstComb <- rbeta(opt$size, shape1=paramFirst[1], shape2=paramFirst[2])

        if(class(paramSecond) != "try-error") {
          psiSecondComb <- rbeta(opt$size, shape1=paramSecond[1], shape2=paramSecond[2])
        }
      }

      # get emperical posterior median of psi
      medOne <- median(psiFirstComb)
      medTwo <- median(psiSecondComb)

      # look for a max difference given prob cutoff...
      opt$prob=0.95
      # psiFirstComb = psiFirstComb[1:3500]
      if(medOne > medTwo) {
        max <- maxDiff(psiFirstComb, psiSecondComb, opt$prob)
      } else {
        max <- maxDiff(psiSecondComb, psiFirstComb, opt$prob)
      }
      #    writeLines(lines[i], stderr()) ### DEBUGGING

      # SIGNIFICANT from here on out:
      if( opt$filter ) {
        write(sprintf("%s\t%s\t%f\t%f\t%f\t%s", tabLine[1], tabLine[2], medOne, medTwo, medOne - medTwo, round(max,2)), stdout())
      }

      # check for significant difference
      if(max < opt$minDiff) { return(NULL) } # or continue...

      writeLines(lines[i], sighandle)

      eventTitle <- paste(c("Gene: ", tabLine[1], "  Event: ", tabLine[2]), collapse="")
      eventCoord <- paste(c("Coordinates: ", tabLine[3]), collapse="")
      #    eventTitleListed[[i]] <- paste(c("Gene: ", tabLine[1], "     ", "Event: ", tabLine[2]), collapse="")

      # Print visual output to pdf;
      if( medOne > medTwo ) {
        retPlot <- plotDiff(psiFirstComb, psiSecondComb, expFirst, expSecond, max, medOne, medTwo, sampOneName, sampTwoName , FALSE)
      } else {
        retPlot <- plotDiff(psiSecondComb, psiFirstComb, expFirst, expSecond, max, medTwo, medOne, sampTwoName, sampOneName , TRUE)
      }

      return(list(retPlot, eventTitle, eventCoord))  #return of mclapply function
  }, mc.cores=opt$cores, mc.preschedule=TRUE, mc.cleanup=TRUE) #End For

  for(it in 1:length(lines)) {
  # PRINT LIST OF PLOTS.
    if(is.null(plotListed[[it]]) || is.null(plotListed[[it]][[1]])) { next; }
    plotPrint(plotListed[[it]][[2]], plotListed[[it]][[3]], plotListed[[it]][[1]])
  }

} #End While

garbage <- dev.off()

close(sighandle)

q(status=0)
