setwd("C:/Users/Dell/Desktop/CNR/Black Sea/anoxic/")



trans<-read.csv2("trans2.csv")
trans2<-trans
vari<-read.csv("variables.csv", sep=";")
library("RColorBrewer")
library(ggplot2)
count<-trans2
rownames(count)<-count[,1]
countdata<-count[,-1]

rownames(vari)<-vari[,1]
vari<-vari[,-1]
library(DESeq2)
cts<-count[,-1]
coldata<-vari

rownames(coldata) <- sub("fb", "", rownames(coldata))



colnames(cts)<-rownames(coldata)

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment)
dds
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$treatment <- factor(dds$treatment, levels = c("DA","DO", "LO"))
dds <- DESeq(dds)
res <- results(dds)
res


resultsNames(dds)
resLFC <- lfcShrink(dds, coef="treatment_LO_vs_DA", type="apeglm")
resLFC
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)
res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.01, na.rm=TRUE)
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

#idx <- identify(res$baseMean, res$log2FoldChange)
#rownames(res)[idx]
resultsNames(dds)
resNorm <- lfcShrink(dds, coef=2, type="normal")
library(ashr)
resAsh <- lfcShrink(dds, coef=2, type="ashr")
par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")
plotCounts(dds, gene=which.min(res$padj), intgroup="light")
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
head(assay(vsd), 3)
library("vsn")
library("hexbin")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
library("RColorBrewer")
library("vegan")
sampleDistsMatrix<- as.matrix(dist(t(assay(ntd))))

adonis2(sampleDistsMatrix ~ vari$light*vari$Oxigen)
rownames(sampleDistsMatrix) <- paste(vsd$treatment, vsd$type)
colnames(sampleDistsMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Purples")) )(255)
library("pheatmap")
pM1<-pheatmap(sampleDistsMatrix,
              col=colors)
plotPCA(vsd, intgroup="treatment")


pcaData <- plotPCA(vsd, intgroup="treatment", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=treatment)) +
  geom_point(size=5, alpha=0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c("aquamarine4", "gray50", "darkblue"))+
  coord_fixed()


resultsNames(dds)
reults_file<-as.data.frame(res@listData)
results_file<-cbind(res@rownames,reults_file)
colo<-colorspace::diverge_hsv
library("colorspace")

key<-read.csv2("key_gene.csv")

res_key<-subset(results_file, results_file$`res@rownames`%in%key$gene)

res_key$type<-rep("key")
results_file$type<-rep("all")




tog<-rbind(results_file,res_key)
names(res_key)[names(res_key) == 'res@rownames'] <- 'gene'
dkey<-as.data.frame(merge(res_key, key, by="gene"))
dkey$rep<-rep("1")

library("cowplot")

pAD<-ggplot(data = tog,aes(x = type, y = log2FoldChange))+
  geom_jitter(aes(color=pvalue < 0.01, fill=pvalue < 0.01), shape = 21,size=2, alpha=0.4)+ 
  scale_color_manual(values = c("TRUE" = "grey10", "FALSE" = "grey78"))+ 
  scale_fill_manual(values = c("FALSE" = "grey10", "TRUE" = "darkred"))+ 
  geom_violin(alpha=0.4,size=1,color="grey", position = position_dodge(width = .75)) +
  guides(color="none")+
  ylab(  c("log2 Fold Change LO vs DA")  )  +
  xlab(  c("")  )  +
  theme_minimal_hgrid(12) +
  theme(axis.ticks = element_line(size=0.1,color="grey"),
        axis.ticks.length=unit(0.2,"cm"),
        legend.position = "none")
legend3<-get_legend((pAD+theme(legend.position = "bottom")))
names(res_key)[names(res_key) == 'res@rownames'] <- 'gene'
res_key$`gene` = factor(res_key$`gene`,levels=res_key$`gene`[order(res_key$log2FoldChange)])

pD<-ggplot(res_key,aes(x=log2FoldChange,y=gene,fill=pvalue<0.01))+
  geom_col(alpha=0.4) + 
  scale_fill_manual(values=c("grey30", "darkred"))+
  labs(fill="pvalue<0.01")+
  theme_minimal_hgrid(12) +
  theme(legend.position = "none",
        axis.text.y=element_blank(),
        axis.title.y = element_blank(),
        axis.title.x=element_blank()
  )+
  geom_errorbar(aes(xmin=log2FoldChange-lfcSE, xmax=log2FoldChange+lfcSE), width=0,
                position=position_dodge(0.9)) 


legend <- get_legend(pD + theme(legend.position = 'bottom'))



pG<-ggplot(dkey,aes(x=rep,y=gene,color=Genetype))+
  geom_point(size=3)+
  scale_color_manual(values=c("gold2","sienna2", "orchid3", "aquamarine3","olivedrab","grey34"))+
  guides(color=guide_legend(ncol=2))+
  theme_minimal_hgrid(12) +
  theme(legend.position="none",
        legend.justification = "left",
        axis.text.x=element_blank(),
        axis.title.x = element_blank())


legend2<-get_legend((pG+theme(legend.position = "bottom")))
library("cowplot")

grid1<-plot_grid(pAD,pG, 
                 pD, ncol=3,
                 align="h", rel_widths=c(1,2,2))

grid2<- plot_grid(legend3, legend2,
                  legend,
                  ncol=3,
                  rel_widths = c(1,2,1))

plot_grid(grid1,grid2,ncol=1, rel_heights = c(15,1))




res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.01, na.rm=TRUE)
plotMA(res05, ylim=c(-2,2))
res05<-subset(results_file, results_file$padj<0.001)
res05m<-subset(res05, res05$log2FoldChange<0)

res05m$`res@rownames` = factor(res05m$`res@rownames`,levels=res05m$`res@rownames`[order(res05m$log2FoldChange)])
ggplot(res05m,aes(x=`res@rownames`,y=log2FoldChange,fill=log2FoldChange>0))+
  geom_col() + coord_flip()+
  scale_fill_manual(values=c("grey10","grey59"),
                    labels=c("negative","positive"))+
  labs(fill="LO vs DA", title="LO vs DA <0 (p<0.001)")+
  theme_minimal_hgrid(12) +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x=element_blank())





