#set path to source directory, change this if you want to run the code
setwd("reproducibility/peer4rnaseq/")

#Plot expression values distribution before and after correction
ori<-t(read.table(bzfile("tsiexp.txt.bz2"),sep=","))
def<-t(read.table(bzfile("tsiexp.def.txt.bz2"),sep=","))

n<-19176

all<-data.frame(values=c(rowMeans(ori),rowMeans(def)),
                method=c(rep("original",n),rep("corrected",n)))


ggplot(all,aes(x=values))+
  geom_histogram()+
  theme_bw(base_size = 14)+
  facet_wrap(~method,nrow=2)

#plot alpha values: how important are factors in the model
ori.a<-t(read.table("tsiexp.alpha.def.txt",sep=","))
alpha<-1.0 / ori.a
barplot(alpha,names.arg=1:8,xlab="Factors",ylab="Importance")

#For futures post
# a1<-read.table("yrinew.alpha_0.01_1.txt")
# a2<-read.table("yrinew.alpha_0.1_10.txt")
# a3<-read.table("yrinew.alpha_1_100.txt")
# 
# 
# all<-data.frame(values=c(rowMeans(ori),rowMeans(def),rowMeans(a1),rowMeans(a2),rowMeans(a3)),
#                 method=c(rep("ori",n),rep("def",n),rep("a1",n),rep("a2",n),rep("a3",n)))
# 
# ggplot(all,aes(x=values))+
#   geom_histogram()+
#   facet_wrap(~method,nrow=2)
# 
# 
# b1<-read.table("yrinew.alpha_0.1_10Eps_0.01_1.txt")
# b2<-read.table("yrinew.alpha_0.1_10Eps_1_100.txt")
# b3<-read.table("yrinew.alpha_0.1_10Eps_10_1000.txt")
# 
# all2<-data.frame(values=c(rowMeans(ori),rowMeans(def),rowMeans(b1),rowMeans(b2),rowMeans(b3)),
#                 method=c(rep("ori",n),rep("def",n),rep("b1",n),rep("b2",n),rep("b3",n)))
# 
# ggplot(all2,aes(x=values))+
#   geom_histogram()+
#   facet_wrap(~method,nrow=2)


