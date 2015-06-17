pMim - Integrating RNA-Seq, miRNA-Seq and pathway information
========================================================

```{r include=FALSE}
opts_chunk$set(message=FALSE, warning=FALSE)
```


The following provides an example of how to use pMim. The markdown code that generated this html file, <a href="https://github.com/ellispatrick/sydSeq/blob/master/Examples/pMim/pMimExample.Rmd" target="_blank"> pMimExample.Rmd </a>,  along with the associated data can be found on my github account at <a href="https://github.com/ellispatrick/sydSeq/blob/master/Examples/pMim/" target="_blank"> https://github.com/ellispatrick/sydSeq/blob/master/Examples/pMim/ </a>. Start by reading in the RNA-Seq and miRNA-Seq data.


```{r}
### Load packages
library(devtools)
load_all("~/repos/sydSeq/sydSeq/") #For pMimCor
load_all("~/soft/src/multiMiR/") #For getting miRNA target genes easily
library(goseq) #For getting KEGG pathways easily
library(KEGG.db)
# load_all("~/repos/S") 
library(SimSeq)
library(limma)
### Load Notch2 Knockout counts for use as an example.
n_mirs = 4
fc_multiplier = 1.3
top_kegg=20

load('counts.RData')
#This loaded the gene counts "counts" and the miRNA counts "countsMi".

### For simplicity, while obviously not ideal, for the following we will work with TMM normalized data. We will also restrict to genes with average count over 20.

Data = counts[rowSums(counts)>=20,]
DataMi = countsMi[rowSums(countsMi)>=20,]
classes = substring(colnames(Data),1,2)
names(classes) = colnames(Data)

# We do not have mature miRNA counts, so do the following.
miR = rownames(DataMi)
miR = sub("mir", "miR", miR)
miR = c(miR, paste(miR, "-3p", sep = ""), paste(miR, "-5p", sep = ""))
miR = split(miR, c(rownames(DataMi), rownames(DataMi), rownames(DataMi)))

mirdb = read.table("~/repos/pcbc_c4_analysis/data/targets/mmu_miRDB_v5.0_prediction_result.txt")
map = read.table("~/repos/pcbc_c4_analysis/data/targets/mmu_refseq_mirdb_ensembl")
idx = match(mirdb$V2, map$V1)
mirdb = cbind(mirdb, gene=map[idx,"V2"])


kegg = getgo(rownames(Data), fetch.cats = "KEGG", genome = "mm10", id = "ensGene")
kegg2ens = Biobase::reverseSplit(kegg)
kname = mget(names(kegg2ens), env = KEGGPATHID2NAME)
names(kegg2ens) = kname

get_DE_genes = function(top_kegg){
    ###KEGG can be loaded as follows
    DEgenes = intersect(kegg2ens[[1]], rownames(Data))[1:top_kegg]
    
    
    kegg2ens_no = 
        lapply(kegg2ens, function(x){
            if (length(intersect(x,DEgenes))<1){
                return(x)
            }
        })
    
    kegg2ens_no = c(Filter(length, kegg2ens_no), positive=kegg2ens[[1]])
    others = sapply(kegg2ens,function(x){ length(intersect(x, DEgenes))})
    others = others[others>0]
    list(DEgenes=DEgenes, others=others)
}
#kegg2ens = c(kegg2ens[[1]],)
```

```{r}

targets = lapply(miR, function(z) as.character(unique(mirdb[mirdb[, "V1"] %in% z, "gene"])))
#targets is our list of miRNA targets

```

```{r sim-function}
getSim = function(DEgenes, n_mirs=4, fc_multiplier=1.3){
    simRna = lapply(rowMeans(Data),function(x){
        rnbinom(6, size=x, prob=0.4)
    })
    simRna = do.call(rbind,simRna)
    simRna[DEgenes,1:3] = simRna[DEgenes,1:3] * fc_multiplier
    
    y = voom(counts = simRna, design = model.matrix(~classes))
    simRnanorm = y$E
    fit = lmFit(y, design = model.matrix(~classes))
    fit = eBayes(fit)
    res = topTable(fit,coef=2,number = Inf)
    dim(res[res$adj.P.Val<0.1,])
    cat("% of real DE genes detected", sum(rownames(res[res$adj.P.Val<0.1,]) %in% DEgenes)/length(DEgenes), "\n")

    mirs = mirdb$gene %in% DEgenes
    DEmirs = names( summary(as.factor(mirdb[mirs,1])) )
    
    set.seed(42)
    DEmirs = unlist(sapply(miR, function(z){ 
        if (sum(DEmirs %in% z)>0) {return(z[1])}
    }))
    n_mirs = min(n_mirs, length(DEmirs))
    cat("mirs selected: ", n_mirs, "\n")
    DEmirs = sample(names(DEmirs), n_mirs)
    
    real_targets = unlist(sapply(DEmirs,function(x){
        intersect(targets[[x]], DEgenes)
    }))
    n_real_targets =  length(unique(unlist(real_targets)))
    
    #simMirs = SimData(DataMi, treatment=classes, genes.select = 1:nrow(DataMi), genes.diff=which(rownames(DataMi) %in% DEmirs), power=0.8, sort.method="unpaired", k.ind=3,switch.trt = TRUE)
    simMi = lapply(rowMeans(DataMi),function(x){
        rnbinom(6, size=x, prob=0.4)
    })
    simMi = do.call(rbind,simMi)
    simMi[DEmirs,4:6] = simMi[DEmirs,4:6] * fc_multiplier
    
    y = voom(counts = simMi, design = model.matrix(~classes))
    simMinorm = y$E
    fit = lmFit(y, design = model.matrix(~classes))
    fit = eBayes(fit)
    res = topTable(fit,coef=2,number = Inf)
    dim(res[res$adj.P.Val<0.1,])
    cat("% of real DE miRNAs detected", sum(rownames(res[res$adj.P.Val<0.1,]) %in% DEmirs)/length(DEmirs),"\n")
    simDEmirs = rownames(res[res$adj.P.Val<0.1,])
    
    colnames(simRnanorm) = names(classes)
    colnames(simMinorm) = names(classes)
    
    list(simrna=simRnanorm, simmirna=simMinorm, n_real_targets=n_real_targets, DEgenes=DEgenes, DEmirna=simDEmirs)
    
}
```


The pathway and miRNA target information can then be read in. Here we use KEGG for the pathway information and TargetScan for the miRNA target predictions.


We are now in a position to run pMim and view the top ten mir-pathways that are identified.

```{r}
options = expand.grid(c(2,5,10,20), c(4,10,20,30), c(0.5,1,1.3,1.5))
meta_sim = apply(options, 1, function(params){
    nde = params[1]
    nmir = params[2]
    multiplier = params[3]
    DEgenes = get_DE_genes(nde)
    sim = getSim(DEgenes$DEgenes,nmir,multiplier)
    cat("real DE got it:", sim$n_real_targets, "\n" )
    
    output = pMim(sim$DEmirna, dtMi = sim$simmirna, dtG = sim$simrna, classes = classes, targets = targets, pathways = kegg2ens, fdr="bonferroni")
    found=NA
    if (nrow(output)>0){
        rownames(output) = output$path
        output$true = 0
        output[names(others),"true"] = others
        output$score = output$genes/output$weight
        
        resEnrichment = output
        idx = which(resEnrichment$path=="Glycolysis / Gluconeogenesis")
        found=resEnrichment[idx,"FDR"]
    }
    data.frame(init_DE=nde, fc=multiplier, real_DE=sim$n_real_targets, ini_mir = nmir, n_DE_mirs=length(sim$DEmirna), found=found)
})
# output %>% filter(true>0) %>% arrange(score)
do.call(rbind,meta_sim)
```


```{r real data}

y = voom(counts = DataMi, design = model.matrix(~classes))
sRnanorm = y$E
fit = lmFit(y, design = model.matrix(~classes))
fit = eBayes(fit)
res = topTable(fit,coef=2,number = Inf)
real_DE_mirs = rownames(res[res$adj.P.Val<0.3,])

y = voom(counts = Data, design = model.matrix(~classes))
Rnanorm = y$E

colnames(sRnanorm) = names(classes)
colnames(Rnanorm) = names(classes)

real_output = pMim(real_DE_mirs, dtMi = sRnanorm, dtG = Rnanorm, classes = classes, targets = targets, pathways = kegg2ens, fdr="bonferroni")

```

We can also output the results into table which can be embedded in html code or their own html files.
An example of the output is give in <a href="pMimResults.html" target="_blank"> pMimResults.html </a>.

Clicking on the "genes" links will take you to the genes that lie in the selected mir-pathway. These tables generated by the googleVis package and can be sorted by clicking on the column headers. 


```{r results='asis'}
load('GeneSymbol.RData')

dir = 'mirpathways'
cutoff = 0.05
filename = 'pMimResults.html'
GeneSymbol = Symbol

htmlOutput = pMimHTML(output,dir = 'mirpathways',cutoff = 0.05,filename = 'pMimResults.html',GeneSymbol = Symbol,outputHTML=TRUE)
print(htmlOutput,type = 'html') 



```

