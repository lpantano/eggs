
## ----setup---------------------------------------------------------------
library(knitr)
library(rmarkdown)
library(knitrBootstrap)
library(ggplot2)
library(reshape)
library(dplyr)
library(edgeR)
library(DESeq2)
library(genefilter)
library(CHBUtils)
library(gtools)
theme_set(theme_bw())
opts_chunk$set(tidy=TRUE, highlight=T, figalign="center",
               fig.height=6, fig.width=6, message=F, error=F, warning=F, bootstrap.show.code=FALSE)


## ----counts--------------------------------------------------------------
counts <- read.table("~/projects/ibb/inputs/TSI/counts.Genes.tab", header=T, row.names=1)
counts <- counts[rowSums(counts>0)>2,]
phe <- read.table("~/projects/geuvadis/pop_data_withuniqueid.txt",header=T)
tsi <- phe[match(names(counts),phe$hapmap_id),]
tsi$fake <- sample(c("g1","g2"),replace = T,size = ncol(counts))
#dse <- DESeqDataSetFromMatrix(countData = counts, colData = tsi,
#                       design = ~fake)
#dse.rlog <- rlog(dse,blind=TRUE)



## ----samples-------------------------------------------------------------
n <- nrow(tsi)
counts<-counts[rowSums(counts>0)>10,]

do_limma <- function(g1,g2){
    de <- data.frame(con=c(rep("g1",length(g1)),
                            rep("g2",length(g2))),
                     row.names=c(g1,g2))
    de <- de[names(counts),,drop=F]
    dge <- DGEList(counts)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~de$con)
    v <- voom(dge,design)
    lm <- lmFit(v,design)
    lmb <- eBayes(lm)
    res <- topTable(lmb,number = nrow(counts))
    res$adj.P.Val
    #list(sum(res$adj.P.Val < 0.1),min(res$adj.P.Val,na.rm = T))
}

get_estimation <- function(iter,x){
    print(x)
    g1 <- as.character(sample(tsi$hapmap_id,x))
    g2 <- as.character(tsi$hapmap_id[-match(g1,tsi$hapmap_id)])
    do_limma(g1,g2)
}

do_de <- function(size){
    print(size)   
    do.call(cbind,lapply(1:400,get_estimation,size))
}

#tab <- lapply(c(2,5,10,45),do_de)
summary_de <- function(ma){
    ma <- as.data.frame(ma[1:5000,]) 
    rank <- 1:nrow(ma)
    low <- apply(ma,1,quantile, .001) 
    up <- apply(ma,1,quantile, .099) 
    q = data.frame(rank=rank, low=low, up=up)
    l <- predict(loess(low ~ rank, data=q), se =T)
    low <- l$fit - qt(0.975,l$df)*l$se
    up <- l$fit + qt(0.975,l$df)*l$se
    data.frame(rank = q$rank, low=low, up=up)
}
tab.sum <- lapply(tab,summary_de)
#tab <- lapply(c(2),do_de)
#tab.join <- do.call(rbind,tab)


