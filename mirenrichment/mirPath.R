# The preparation code comes from the package pMim from https://github.com/lpantano/sydSeq.git
# I changed many parts of the code and stats, but the inspiration comes from that

require(limma)

mirPath = function (DEmirs=DEmirs, dtMi, dtG, classes, targets, pathways, Zmi = NULL,
                 Zg = NULL, cor=0.7, corP = "Fisher", nInt = 2, randomTarg = FALSE,
                 randomPath = FALSE, fdr="fdr", verbose = TRUE)
{
  if (length(DEmirs)==0){
      message("No deregulated miRNAs")
      return(data.frame())
  }
  if (verbose == TRUE) {
    if (sum(colnames(dtG) != names(classes)) > 0 | is.null(names(classes))) {
      cat("colnames of Gene Data do not match names of classes")
    }
  }
  if (verbose == TRUE) {
    if (sum(colnames(dtG) != names(classes)) > 0 | is.null(names(classes))) {
      cat("colnames of miRNA Data do not match names of classes")
    }
  }
  if (!is.null(cor)){
    if (is.matrix(cor)){
      if (ncol(cor) != nrow(dtG) | nrow(row) != nrow(dtMi)){
        stop("The correlation matrix should be same rows than miRNAs and same columns than genes.")
      }
      TR = cor
    } else{
      if (is.numeric(cor)){
        cor = abs(cor)
      }
    }
  } else{
    stop("Not supported correlation value. Should be NULL, matrix, or numeric.")
  }


  if (ncol(dtG) == ncol(dtMi)){
    y = names(classes)
    X = t(dtG)
    Y = t(dtMi)
    corYX = cor(Y[y, ], X[y, ])
    n = length(y)
    TR = corYX^2
  }else{
    warning("Not paired samples, all miRNA-mRNA will be considered equally")
    TR = matrix(ncol=ncol(dtG), nrow=nrow(dtMi))
    TR[] = cor
  }

  if (is.null(Zmi)) {
    require(limma)
    design = model.matrix(~classes)
    Dat = dtMi[, names(classes)]
    fit = lmFit(Dat, design)
    #ordinary.t <- fit$coef/fit$stdev.unscaled/pmax(pmax(fit$sigma,
    #                                                    sqrt(rowMeans(Dat))), 1)
    efit =ebayes(fit)
    #Tmi = ordinary.t[, 2]
    Zmi = qnorm(pt(efit$t[,2], fit$df.res))
  }

  if (is.null(Zg)) {
    require(limma)
    design = model.matrix(~classes)
    Dat = dtG[, names(classes)]
    fit = lmFit(Dat, design)
    #ordinary.t <- fit$coef/fit$stdev.unscaled/pmax(pmax(fit$sigma,
    #                                                    sqrt(rowMeans(Dat))), 1)
    efit = ebayes(fit)
    #Tg = ordinary.t[, 2]
    Zg = qnorm(pt(efit$t[,2], fit$df.res))
  }

  # mi = colnames(Y)
  mi = DEmirs
  msgr = intersect(unlist(pathways), unlist(targets))
  mapMat = matrix(0, length(mi), length(msgr))
  rownames(mapMat) = mi
  colnames(mapMat) = msgr
  for (i in rownames(mapMat)) {
    if (randomTarg == TRUE) {
      mapMat[i, sample(colnames(mapMat), sum(msgr %in%
                                               targets[[i]]))] = 1
    }
    else {
      mapMat[i, msgr %in% targets[[i]]] = 1
    }
  }
  # mapMat row=mirna, col=genes
  use = names(which(rowSums(mapMat) > 0))
  mapMat = mapMat[use, ]
  mi = use
  pathName = names(pathways)
  pathMat = matrix(0, length(pathName), length(msgr))
  rownames(pathMat) = pathName
  colnames(pathMat) = msgr
  # pathMat row=paths, col=genes
  for (j in 1:length(pathName)) {
    if (randomPath == TRUE) {
      pathMat[j, sample(colnames(mapMat), sum(msgr %in%
                                                pathways[[j]]))] = 1
    }
    else {
      pathMat[j, msgr %in% pathways[[j]]] = 1
    }
  }
  N = mapMat %*% t(pathMat)
  test = which(N >= nInt, 2)
  testPath = unique(test[,2])
  G = rep(NA,length=length(testPath))
  names(G) = pathName[testPath]
  #G = G[1:10]
  N = Nw = G
  #colnames(G) = pathname
  for (path_item in testPath){
    mir_items = rownames(test[test[,2]==path_item,,drop=F])
    name_path = rownames(pathMat)[path_item]
    pathway = pathMat[name_path, ]
    binding = round((colSums(mapMat[mir_items, ,drop=F]) / (colSums(mapMat[mir_items, ,drop=F ]) + 1)))
    keep = names(which(binding * pathway == 1))
    Zg_keep = Zg[keep]
    Zm_keep = Zmi[mir_items]

    if (length(keep)>0){
        t = TR[mir_items, keep ,drop=F]
        tmax = sapply(keep, function(x){
            is_target = names(which(Zm_keep*Zg_keep[x]<0))
            if (length(is_target)==0){return(NA)}
            idx = is_target[which.max(t[is_target, x] * Zg_keep[x])]
            t[idx, x]
        })
        Zg_keep = Zg_keep[!is.na(tmax)]
        tmax = tmax[!is.na(tmax)]

        tcorr = tmax
        tcorr[tcorr<0.5] = 0.5

        q <- pnorm(sum(tmax * abs(Zg_keep))/sqrt(sum(tcorr)));
        q <- 2*min(q, 1-q)

        G[name_path] = q
        N[name_path] = length(tmax)
        Nw[name_path] = sum(tmax)
    }
  }

  G = p.adjust(G, method = fdr)
  # Score = pchisq(-2 * (log(1 - G) + log(1 - Gnorm)), 4)


  data.frame(path=pathName[testPath],FDR=G,genes=N, weight=Nw)
  #list(Results = res, Scores = Score, Zg = Zg, Zmi = Zmi, cor = corYX, corTransform = TR,Direction = direction,targets = targets, pathways = pathways)
}
