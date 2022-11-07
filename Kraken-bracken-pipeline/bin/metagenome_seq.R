#!/usr/bin/env Rscript

library(metagenomeSeq)
library(ggfortify)
library(sva)

# Load data
arg <- commandArgs(trailingOnly = TRUE)

counts<-loadMeta(arg[1],sep="\t")

pdata = loadPhenoData(arg[3])
ord = match(colnames(counts$counts), rownames(pdata)) 
pdata = pdata[ord, ]
phenotypeData = AnnotatedDataFrame(pdata)

taxa = read.delim(arg[2] ,stringsAsFactors = FALSE)
OTUdata = AnnotatedDataFrame(taxa)
data = newMRexperiment(counts$counts,phenoData=phenotypeData,featureData = OTUdata)

# FILTERING
nrow(MRcounts(data))
rareFeatures = which(rowSums(MRcounts(data) <= 0) >= ncol(MRcounts(data)))
if(length(rareFeatures) > 0){
  data = data[-rareFeatures, ]
}

# NORMALIZATION
pdf("CumNormStat.pdf",paper="A4")
datap = cumNormStat(data, pFlag = TRUE, main = "Difference of the median deviance")
dev.off()
data = cumNorm(data, p = datap)
Type<-as.factor(pData(data)[,2])
normFactor = normFactors(data)
normFactor = log2(normFactor/median(normFactor) + 1)

if(ncol(pData(data))>2){
  Covariable<-as.factor(pData(data)[,3])
  mod = model.matrix(~0+Type+Covariable+normFactor)
  
} else {
  mod = model.matrix(~0+Type+normFactor)
}

colnames(mod)=gsub("Type|Covariable","",colnames(mod))
colnames(mod)=gsub("-","_",colnames(mod))

res <- tryCatch( {
  settings = zigControl(maxit = 10, verbose = F)
  res=fitZig(obj = data, mod = mod, useCSSoffset = FALSE, control = settings)
  return(res)
},
error= function(cond) {
  message("Using maximum number of iterations = 1")
  settings = zigControl(maxit = 1, verbose = F)
  res=fitZig(obj = data, mod = mod, useCSSoffset = FALSE, control = settings)
  return(res)
}
)

zigFit = slot(res, "fit")
finalMod = slot(res, "fit")$design

# MAKE CONTRAST
contrastconds<-read.delim(file=arg[4], col.names=NA)

for (c in 1:length(rownames(contrastconds))){
  #Obtaining the name of the contrast and the conditions to make the contrast
  contrast<-gsub("[)( ]","",strsplit(as.character(contrastconds[c,1]), "=")[[1]])

  #Performing the DE analysis for each pair
  contrast.matrix = makeContrasts( contrasts= contrast[2],
                                   levels = finalMod)
  fit = contrasts.fit(zigFit, contrast.matrix)
  fit = eBayes(fit)
  top<-topTable(fit, number = nrow(counts$counts))
  
  Comp<- merge(top, taxa, by=0, all=F)
  colnames(Comp)[1]<-"Taxid"
  write.table(Comp[order(Comp$adj.P.Val),],file=paste0(contrast[1],".xls"),
              quote = F,sep="\t",row.names =F)
  
}


# CLUSTERING PLOTS
NormObj<-MRcounts(data,norm = T)
LogNormObj<-MRcounts(data,norm = T, log=T)

pdf("Clustering_plots.pdf",paper="A4")

# Top 100 entities
rsd <- rowSds(as.matrix(NormObj))
sel <- order(rsd, decreasing=TRUE)[1:100]
heatmap(na.omit(as.matrix(NormObj[sel,])),margins=c(10,8),
        main="Heatmap 100 most DE entities",cexRow=0.01,cexCol=0.5,
        labCol=rownames(pdata))

# Dendrogram
pr.hc.c<- hclust(na.omit(dist(t(NormObj),method = "euclidean")))
plot(pr.hc.c, xlab="Sample Distance",main=("Hierarchical Clustering of Normalized Samples"), 
     labels=pdata[,1], cex = 0.7)

#PCoA log

col.group <- Type
levels(col.group) <-  palette.colors(nlevels(col.group))
col.group <- as.character(col.group)

z<-plotMDS(LogNormObj, labels=pdata[,1], col=col.group, gene.selection = "pairwise", plot=F)
edge<-sd(z$x)
plotMDS(LogNormObj, labels=pdata[,1], col=col.group, gene.selection = "pairwise",
        xlim=c(min(z$x)-edge,max(z$x) + edge)) 
title(main="MDS-PCoA log")

#PCA
data_pca<-as.matrix(LogNormObj)
data_pca<-as.data.frame(t(data_pca))
rownames(data_pca)<-rownames(pdata)
data_pca.PC = prcomp(data_pca)
data_pca$Filename<-colnames(pdata)
data_pca$Name<-pdata[,1]
data_pca$Type<-pdata[,2]

lim=max(abs(data_pca.PC$rotation[,c(1,2)]))+0.2
autoplot(data_pca.PC,label=T,data=data_pca,colour="Type", xlim=c(-lim,lim), 
         label.label = pdata[,1])

dev.off()
