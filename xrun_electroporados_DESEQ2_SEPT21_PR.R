## E. multilocularis, electroporation analysis
## Transcriptional effects of electroporation on Echinococcus multilocularis primary cell culture
## Authors: Matias Gaston Perez, Natalia Rego, Markus Spiliotis, Nancy Holroyd, Klaus Brehm, Mara Cecilia Rosenzvit

## September 16, 2021


##########
## LOAD ##
##########

library(gridExtra)
library(RColorBrewer)
library(ggplot2)

library(DESeq2)
library(pheatmap)
library(vsn)

library(FactoMineR)
library(factoextra)

library(seqinr)

library(biomaRt)

library(data.table)

library(clusterProfiler)
library(enrichplot)



##########
## DATA ##
##########


load("electroporacion_SEPT21.RData")
ls()
#[1] "annots"        "counts"        "pheno"         "polycis.genes"
#[5] "sl.genes"      "trs" 
## counts: counts matrix (expression data obtained using HISAT2, StringTie and prepDE.py)
## pheno: phenotype data frame
## annots: gene annotation from Biomart
## trs: fasta file with E. multilocularis transcripts (Emultilocularis_transcripts.fa)
## sl.genes: list of genes target of splice-leader trans-splicing (from Tsai et al 2003)
## polycis.genes: list of genes transcribed in polycistrons (from Tsai et al 2003)



############
## DESeq2 ##
############


# build dds
dds <- DESeqDataSetFromMatrix(countData=counts, colData=pheno, design= ~ treatment)
dds
#class: DESeqDataSet 
#dim: 10663 6 
#metadata(1): version
#assays(1): counts
#rownames(10663): EmuJ_001059300 EmuJ_000212400 ... EmuJ_002178900
#  EmuJ_000425400
#rowData names(0):
#colnames(6): sEPC1 sEPC2 ... sPC2 sPC3
#colData names(2): sample treatment


# Pre-Filtering
keep <- rowSums(counts(dds)) >=10
sum(keep)
#[1] 9197
dds <- dds[keep,]

dds$treatment <- relevel(dds$treatment, ref="PC")  ## establish PC condition as the reference level

# Differential expression analysis
dds <- DESeq(dds)
dds
#class: DESeqDataSet 
#dim: 9197 6 
#metadata(1): version
#assays(3): counts mu cooks
#rownames(9197): EmuJ_001059300 EmuJ_000784100 ... EmuJ_000676300
#  EmuJ_000690600
#rowData names(21): baseMean baseVar ... deviance maxCooks
#colnames(6): sEPC1 sEPC2 ... sPC2 sPC3
#colData names(3): sample treatment sizeFactor

summary(rowSums(counts(dds)))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#    10     935    5998   12274   14536 1355280 

dds$sizeFactor
#     sCE1      sCE2      sCE3     sCSE1     sCSE2     sCSE3 
#0.8872169 1.0770481 1.0962678 1.0530548 0.9838385 0.9449213

resultsNames(dds)
#[1] "Intercept"               "treatment_EPC_vs_PC"


## PCA (FactoMineR, factoextra) ##  

## transforming data (blind dispersion = TRUE / default)
rld <- rlog(dds)
#p.rld <- meanSdPlot(assay(rld))
## variant stabilizing transformation of data (blind dispersion = TRUE / default)
vsd <- vst(dds)
#p.vsd <- meanSdPlot(assay(vsd))
## this gives log2(n + 1)
ntd <- normTransform(dds)
#p.ntd <- meanSdPlot(assay(ntd))

##get the 500 leading genes (with highest row variance)
myvar <- sort(apply(assay(rld),1,var), decreasing=TRUE)
names(myvar[1:500])
res.pca <- PCA(t(assay(rld)[names(myvar[1:1000]), ]), scale.unit=FALSE, ncp=5, graph=FALSE)

get_eigenvalue(res.pca)
#      eigenvalue variance.percent cumulative.variance.percent
#Dim.1  55.321951        53.997411                    53.99741
#Dim.2  16.256731        15.867506                    69.86492
#Dim.3  13.698552        13.370576                    83.23549
#Dim.4  10.892113        10.631330                    93.86682
#Dim.5   6.283622         6.133177                   100.00000


pca.ell <- fviz_pca_ind(res.pca, geom.ind="point", # show points only (nbut not "text")
pointsize=6, 
col.ind = pheno$treatment, mean.point = FALSE, # color by groups
palette= c("coral1","lightsteelblue3"),
addEllipses = TRUE, ellipse.type="confidence", ellipse.level=0.9, ellipse.alpha=0.1,# Concentration ellipses
legend.title = "treatment", title='') + scale_shape_manual(values=c(19,19)) + xlab("PC1: 54.00% variance") + ylab("PC2: 15.87% variance") + theme_bw()
print(pca.ell)


pdf("PCA_con_ellipses_SEPT21.pdf", width=16, height=16)
print(pca.ell)
dev.off()


## heatmap of among-sample distances

df.col <- as.data.frame(pheno[,2])
colnames(df.col) <- "treatment"
rownames(df.col) <- rownames(pheno)

my.colors <- c("coral1","lightsteelblue3")
names(my.colors) <- c("EPC","PC")
mycols.de <- list(treatment=my.colors)

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$sample
colnames(sampleDistMatrix) <- rld$treatment
paletteLength <- 255
myColors <- colorRampPalette( rev(brewer.pal(9, "Greys")) )(paletteLength)
pheatmap(t(sampleDistMatrix),
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=myColors,
	 annotation_col=df.col, annotation_colors= mycols.de, annotation_legend=TRUE,
	 cutree_cols=2, cutree_rows=2,
	 cellwidth=30, cellheight=8,
	 fontsize=9, fontsize_row = 6, fontsize_number=0.6,
	 filename="pheatmap_rld_sample-distance.pdf", width=6, height=6)



#################
## DE: EPCvsPC ##
#################


res <- results(dds, name="treatment_EPC_vs_PC")
res
#log2 fold change (MLE): treatment EPC vs PC 
#Wald test p-value: treatment EPC vs PC 
#DataFrame with 9197 rows and 6 columns
#                 baseMean log2FoldChange      lfcSE       stat     pvalue
#                <numeric>      <numeric>  <numeric>  <numeric>  <numeric>
#EmuJ_001059300 11822.6363     0.18613169 0.10507789  1.7713687 0.07649941
#EmuJ_000784100  1676.9094     0.11773894 0.13104568  0.8984573 0.36894181
#...

summary(res$padj)
#   Min.  1st Qu.  Median    Mean  3rd Qu.    Max.    NA's 
#  0.0000  0.1155  0.4401  0.4499  0.7606  0.9998     205



#diagnostic plot: histogram of the p values (this plot is best formed by excluding genes with very small counts, which otherwise generate spikes in the histogram)

pdf("DE_Pvalues.pdf", height=8, width=10)
hist(res$pvalue[res$baseMean > 1], breaks = 0:50/50,
     col = "grey50", border = "white", main=NULL, xlab="p-value", cex=1.5, cex.lab=1.3)				
dev.off()


summary(res$log2FoldChange)
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -25.26234  -0.19617   0.00442   0.03600   0.20336  25.75525


res05 <- results(dds, name="treatment_EPC_vs_PC", alpha=0.05)
summary(res05)
# out of 9197 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 744, 8.1%
# LFC < 0 (down)     : 901, 9.8%
# outliers [1]       : 26, 0.28%
# low counts [2]     : 179, 1.9%
# (mean count < 3)

## have a look on a combination of padj and log2FoldChange
sum(res05$padj <=0.05 & abs(res$log2FoldChange) >=1, na.rm=TRUE)
#[1] 254
sum(res05$padj <=0.05, na.rm=TRUE)
##[1] 1645


#provide threshold for constructing Wald tests of significance (is abs(log2(FC))>=1?)
res05LFC1 <- results(dds, name="treatment_EPC_vs_PC", alpha=0.05, altHypothesis="greaterAbs", lfcThreshold=1)

summary(res05LFC1)
#out of 9197 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 1.00 (up)    : 29, 0.32%
#LFC < -1.00 (down) : 4, 0.043%
#outliers [1]       : 26, 0.28%
#low counts [2]     : 0, 0%
#(mean count < 2)

dim(subset(res05LFC1, abs(log2FoldChange)>=1 & padj<=0.05))
#[1] 33  6
genes.LFC1 <- rownames(subset(res05LFC1, abs(log2FoldChange)>=1 & padj<=0.05))

## shrunken log fold changes using apeglm 
res.apeglm <- lfcShrink(dds, coef="treatment_EPC_vs_PC", type="apeglm")
res.apeglm.1 <- lfcShrink(dds, coef="treatment_EPC_vs_PC", type="apeglm", lfcThreshold=1) 

dim(subset(res.apeglm, abs(log2FoldChange)>=1))
#[1] 174   5
dim(subset(res.apeglm, abs(log2FoldChange)>=1 & padj<=0.05))
#[1] 174   5

summary(res.apeglm.1$svalue)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0000  0.9123  0.9562  0.8972  0.9708  0.9781
genes.apeglm1 <- rownames(subset(res.apeglm.1, abs(log2FoldChange)>=1 & svalue<=0.05))


genesDE <- rownames(res[res$padj<=0.05 & !is.na(res$padj),])
length(genesDE)
#[1] 1645


## log2FC summary
summary(res[genesDE,]$log2FoldChange)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -11.24679  -0.53356  -0.28546   0.05299   0.43755  12.91472

## shrunken log2FC summary (apeglm)
summary(res.apeglm[genesDE,]$log2FoldChange)
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
# -10.025932  -0.459604  -0.248467   0.005539   0.363996  12.090637


## plotMA with both FoldChange treatments (raw y shrinkage "apeglm") ##REPTIR POR PLOT
pdf("MAplots.pdf", width=15, height=10)
par(mfrow=c(2,2), mar=c(4,4,2,1))
xlim <- c(1,1e6) ; ylim <-c(-13,13)
drawLines <- function() abline(h=c(-1,1),col="dodgerblue",lwd=2)
plotMA(res, alpha=0.05, xlim=xlim, ylim=ylim, main="raw"); drawLines()
plotMA(res05LFC1, alpha=0.05, xlim=xlim, ylim=ylim, main="raw abs(log2(FC))>=1"); drawLines()
plotMA(res.apeglm, alpha=0.05, xlim=xlim, ylim=ylim, main="apeglm"); drawLines()
plotMA(res.apeglm.1, alpha=0.05, xlim=xlim, ylim=ylim, main="apeglm abs(log2(FC))>=1"); drawLines()
dev.off()


## write table DE 
res.out <- as.data.frame(res)
res.out$apeglm <- as.vector(res.apeglm$log2FoldChange)
res.out$apeglm.svalue <- as.vector(res.apeglm.1$svalue)
res.out$LFC1 <- as.vector(res05LFC1$padj)
write.table(res.out,"genes_log2FC_shrinkage-LFC1.csv", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE, dec=",")


## VOLCANO PLOT APEGLM  

topT <- as.data.frame(res.apeglm)
pdf("volcano_shrinkage-apeglm.pdf", width=10, height=10)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, cex=1.0, xlab=bquote(~log[2]~(fold~change)), ylab=bquote(~-log[10]~(adjusted~p~value))))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(padj), pch=20, col="coral1", cex=0.9))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-1, col="black", lty=4, lwd=1.5)
abline(v=1, col="black", lty=4, lwd=1.5)
abline(h=-log10(0.05), col="black", lty=4, lwd=1.5)
dev.off()


## Heatmap 174 genes (apeglm, top DE)  

## genes DE
genes.de.174 <- rownames(subset(topT, padj<0.05 & abs(log2FoldChange)>1))

##assay(rld)
rld.n <- rlog(dds, blind=FALSE)

df <- as.data.frame(colData(dds)[,c("treatment")])
colnames(df) <- c("treatment")
rownames(df) <- colnames(dds)

my.colors <- c("coral1","lightsteelblue3")
names(my.colors) <- c("EPC","PC")
mycols.de <- list(treatment=my.colors)

paletteLength <- 150
myPal <- colorRampPalette(c("blue","white","red"))(paletteLength)

pheatmap(assay(rld.n)[genes.de.174,], color=myPal, scale="row", cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, cutree_cols=2, cutree_rows=2, annotation_col=df, annotation_colors= mycols.de, annotation_legend=TRUE, legend=TRUE, cellwidth=30, cellheight=8, fontsize=9, fontsize_row = 6, fontsize_number=0.6, filename="pheatmap_rld_174-DEgenes.pdf", width=15, height=25)

## genes UP y DOWN
genes.de.UP <- rownames(subset(topT, padj<0.05 & log2FoldChange>=1))
length(genes.de.UP)
#[1] 105
genes.de.DOWN <- rownames(subset(topT, padj<0.05 & log2FoldChange<=-1))
length(genes.de.DOWN)
#[1] 69



## GET TPMs VALUES

## read fasta file
#trs <- read.fasta(file="Emultilocularis_transcripts.fa")
### how many seqs
length(trs)
#[1] 10669
### getLength
trs.length <- getLength(trs)
head(trs.length)
#[1] 1260 5271 1866  498  378 1983
names(trs.length) <- names(trs)

## easy tpm function (https://gist.github.com/slowkow/c6ab0348747f86e2748b)
tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

dim(counts(dds))  ## normalized=FALSE
#[1] 9197    6
length(trs.length[rownames(counts(dds))])
#[1] 9197
trs.length.filt <- trs.length[rownames(counts(dds))]

tpms <- tpm(counts(dds),trs.length.filt)
dim(tpms)
#[1] 9197    6
head(tpms,2)
#                   sEPC1     sEPC2     sEPC3    sPC1    sPC2    sPC3
#EmuJ_001059300 48.78847 52.71058 55.86551 48.19769 49.39976 37.94112
#EmuJ_000784100 20.88623 28.13075 26.61971 27.62236 18.72559 21.62057

## per sample transcripts with TPM>=1
colSums(tpms >=1)
#sEPC1 sEPC2 sEPC3  sPC1  sPC2  sPC3 
# 7202  7327  7360  7339  7303  7264 


tpm.mean.EPC <- rowSums(tpms[,1:3])/3
tpm.mean.PC <- rowSums(tpms[,4:6])/3

tpms.df <- as.data.frame(tpms)
tpms.df$tpm.mean.EPC <- tpm.mean.EPC
tpms.df$tpm.mean.PC <- tpm.mean.PC
head(tpms.df,2)
#                  sEPC1    sEPC2    sEPC3     sPC1     sPC2     sPC3
#EmuJ_001059300 48.78847 52.71058 55.86551 48.19769 49.39976 37.94112
#EmuJ_000784100 20.88623 28.13075 26.61971 27.62236 18.72559 21.62057
#               tpm.mean.EPC tpm.mean.PC
#EmuJ_001059300     52.45485    45.17952
#EmuJ_000784100     25.21223    22.65617


dim(subset(tpms.df, tpm.mean.EPC >=1 | tpm.mean.PC >=1))
#[1] 7448    8
genes.tpms.up1 <- rownames(subset(tpms.df, tpm.mean.EPC >=1 | tpm.mean.PC >=1))

tpms.df.up1 <- tpms.df[genes.tpms.up1,]

## boxplots with TPMs values per sample and condition 
df.tpm <- data.frame(TPM=c(tpms.df.up1$sEPC1,tpms.df.up1$sEPC2,tpms.df.up1$sEPC3,tpms.df.up1$sPC1,tpms.df.up1$sPC2,tpms.df.up1$sPC3), samples=rep(c("sEPC1","sEPC2","sEPC3","sPC1","sPC2","sPC3"),each=7448), treatment=rep(c("EPC","PC"),each=7448*3))

p1 <- ggplot(df.tpm, aes(x=samples, y=log10(TPM+1), fill=treatment)) + geom_boxplot() + scale_fill_manual(values=c("coral1","lightsteelblue3")) + scale_x_discrete(limits=c("sPC1","sPC2","sPC3","sEPC1","sEPC2","sEPC3")) + ylab("log10 (TPM + 1)") + theme_bw()


## sccaterplot mean TPMs between conditions
p2 <- ggplot(tpms.df, aes(x=log10(tpm.mean.PC+1), y=log10(tpm.mean.EPC+1))) + geom_point(alpha=0.4) + xlab("log10 (mean TPM + 1), PC") + ylab("log10 (mean TPM + 1), EPC") + theme_bw()


pdf("TPMs_all_samples_description_7448genesTPM1.pdf", width=15, height=8)
grid.arrange(p1,p2, nrow=1)
dev.off()


##################################################
## GENE EXPRESSION IN PC (primary cell culture) ##
##################################################


tpm.mean.PC <- rowSums(tpms[,4:6])/3
head(tpm.mean.PC,3)
#EmuJ_001059300 EmuJ_000784100 EmuJ_000056900 
#     45.179524      22.656173       2.635898

tpm.PC.order <- names(sort(tpm.mean.PC, decreasing=TRUE))
head(tpm.mean.PC[tpm.PC.order])
#EmuJ_000742900 EmuJ_000381200 EmuJ_000292700 EmuJ_000381500 EmuJ_000982200 
#     3506.1084      2184.4122      1440.6496      1249.8287       908.0712 
#EmuJ_000036300 
#      906.1910

sum(tpm.mean.PC >=1)
#[1] 7324
summary(tpm.mean.PC[tpm.mean.PC >=1])
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#   1.001    5.337   11.129   22.793   22.149 3506.108

quantile(tpm.mean.PC[tpm.mean.PC >=1], probs = seq(0, 1, 0.050))
dim(subset(tpms.df, tpm.mean.PC>=50))
#[1] 577   8
tpms.PC.UP50 <- rownames(subset(tpms.df, tpm.mean.PC>=50))


## top expressed genes in PC that also show DE with electroporation treatment
sum(tpms.PC.UP50 %in% genesDE)
#[1] 202
sum(tpms.PC.UP50 %in% genes.LFC1)
#[1] 2
tpms.PC.UP50[tpms.PC.UP50 %in% genes.LFC1]
#[1] "EmuJ_000127000" "EmuJ_001030300"

annots[tpms.PC.UP50[tpms.PC.UP50 %in% genes.LFC1],]
#                       GeneID               Chromosome GeneStart  GeneEnd
#EmuJ_000127000 EmuJ_000127000 pathogen_EmW_scaffold_04   9607547  9608401
#EmuJ_001030300 EmuJ_001030300 pathogen_EmW_scaffold_02  10108728 10112104
#               Strand GeneName
#EmuJ_000127000      1         
#EmuJ_001030300      1         
#                                                                                            description
#EmuJ_000127000                                                                                         
#EmuJ_001030300 Basic leucine zipper bZIP transcription factor  [Source:UniProtKB/TrEMBL;Acc:A0A068YCT9]
#                      biotype
#EmuJ_000127000 protein_coding
#EmuJ_001030300 protein_coding

res.apeglm[tpms.PC.UP50[tpms.PC.UP50 %in% genes.LFC1],]
#log2 fold change (MAP): treatment CEle vs CSele 
#Wald test p-value: treatment CEle vs CSele 
#DataFrame with 2 rows and 5 columns
#                baseMean log2FoldChange     lfcSE      pvalue        padj
#               <numeric>      <numeric> <numeric>   <numeric>   <numeric>
#EmuJ_000127000  16869.13        1.65755  0.100987 7.57532e-62 6.81173e-58
#EmuJ_001030300   8693.44        1.67751  0.129722 1.69465e-39 5.07942e-36

summary(res05[tpms.PC.UP50[tpms.PC.UP50 %in% genesDE],]$log2FoldChange)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1.2882 -0.4729 -0.3265 -0.2533 -0.2641  1.6966

dim(res.out[tpms.PC.UP50[tpms.PC.UP50 %in% genesDE],])
#[1] 202  9
res.out[tpms.PC.UP50[tpms.PC.UP50 %in% genesDE],]
write.table(res.out[tpms.PC.UP50[tpms.PC.UP50 %in% genesDE],],"genes_TPM_above50_in_DE.csv", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE, dec=",")


## FIG 2
## TPM density plot (highlight TPM50) and geom_segment for specific gene lists

tpm.PC.increase <- sort(tpm.mean.PC)
length(tpm.PC.increase)
#[1] 9197
tpm.PC.increase <- tpm.PC.increase[tpm.PC.increase>=1]
length(tpm.PC.increase)
#[1] 7324

agB <- c("EmuJ_000381200","EmuJ_000381400","EmuJ_000381100","EmuJ_000381500","EmuJ_000381700","EmuJ_000381600","EmuJ_000381800")
length(agB)
#7
length(agB[agB %in% names(tpm.PC.increase)])
#7
#agB <- agB[agB %in% names(tpm.PC.increase)]

fabp <- c("EmuJ_000550000","EmuJ_002165500","EmuJ_000549800","EmuJ_000417200")
length(fabp)
#4
length(fabp[fabp %in% names(tpm.PC.increase)])
#4
#fabp <- fabp[fabp %in% names(tpm.PC.increase)]

histone <- c("EmuJ_000906000","EmuJ_000579800","EmuJ_000472800","EmuJ_000128000","EmuJ_000927700","EmuJ_000439100","EmuJ_001075100","EmuJ_000472900","EmuJ_000472600","EmuJ_000209600",
"EmuJ_000338200","EmuJ_001031500","EmuJ_001142500","EmuJ_001020600","EmuJ_001102300","EmuJ_001016300","EmuJ_000323000","EmuJ_000323100","EmuJ_000688200","EmuJ_001123800",
"EmuJ_000660300","EmuJ_000611800","EmuJ_001181900","EmuJ_000421400","EmuJ_000103800","EmuJ_002153000","EmuJ_000184700","EmuJ_000223500","EmuJ_001114800","EmuJ_001005500",
"EmuJ_000223900","EmuJ_001170200","EmuJ_000970900","EmuJ_000099900","EmuJ_000452300","EmuJ_000674500","EmuJ_000101600","EmuJ_000489600","EmuJ_000750400","EmuJ_000439200",
"EmuJ_000633200","EmuJ_000897900","EmuJ_000922800","EmuJ_000970800","EmuJ_000898000","EmuJ_000105000","EmuJ_000899100","EmuJ_000362900","EmuJ_000993000","EmuJ_000323300",
"EmuJ_000322800","EmuJ_001170100","EmuJ_000210300","EmuJ_000171900","EmuJ_000656800","EmuJ_000971000","EmuJ_000750900","EmuJ_000527350","EmuJ_000918000","EmuJ_000751500",
"EmuJ_000095900","EmuJ_001167700","EmuJ_000566500","EmuJ_002157700")
length(histone)
#64
length(histone[histone %in% names(tpm.PC.increase)])
#63
histone <- histone[histone %in% names(tpm.PC.increase)]

tubulin <- c("EmuJ_000886400","EmuJ_000672200","EmuJ_000413200","EmuJ_000202600","EmuJ_001013400","EmuJ_001099700","EmuJ_000202500","EmuJ_000975200","EmuJ_000718000",
"EmuJ_000096900","EmuJ_000413000","EmuJ_000041800","EmuJ_000325100","EmuJ_000603100","EmuJ_000576100","EmuJ_000616050","EmuJ_000877100","EmuJ_000040300",
"EmuJ_001051800","EmuJ_000075600","EmuJ_001053400","EmuJ_000621500","EmuJ_000922600","EmuJ_001126150","EmuJ_000653400","EmuJ_000581400","EmuJ_000040900",
"EmuJ_000346100","EmuJ_000476400","EmuJ_001125100","EmuJ_000569000","EmuJ_000339900","EmuJ_000041100","EmuJ_000042500","EmuJ_000042200")
length(tubulin)
#35
length(tubulin[tubulin %in% names(tpm.PC.increase)])
#32
tubulin <- tubulin[tubulin %in% names(tpm.PC.increase)]

tetraspanin <- c("EmuJ_001077100","EmuJ_001021300","EmuJ_000355800","EmuJ_000957800","EmuJ_000355500","EmuJ_000354700","EmuJ_000834300","EmuJ_000355400","EmuJ_000355200",
"EmuJ_000355900","EmuJ_000355100","EmuJ_001021500","EmuJ_000833400","EmuJ_000989900","EmuJ_000977700","EmuJ_000363200","EmuJ_000355700","EmuJ_001019200","EmuJ_000355300",
"EmuJ_001065000","EmuJ_001174600","EmuJ_000858500","EmuJ_001021700","EmuJ_000328400","EmuJ_000117700","EmuJ_001077500","EmuJ_000356000","EmuJ_000354850","EmuJ_000354900",
"EmuJ_001077400","EmuJ_001077200","EmuJ_001077300","EmuJ_000736300","EmuJ_000747400")
length(tetraspanin)
#34
length(tetraspanin[tetraspanin %in% names(tpm.PC.increase)])
#34
#tetraspanin <- tetraspanin[tetraspanin %in% names(tpm.PC.increase)]

hsp70 <- unique(c("EmuJ_000008700","EmuJ_001085400","EmuJ_000320800","EmuJ_001190900","EmuJ_000249600","EmuJ_001085100","EmuJ_001065400","EmuJ_001202400","EmuJ_001154300","EmuJ_001154200","EmuJ_001154000","EmuJ_001153700","EmuJ_001153600","EmuJ_000736500","EmuJ_000733500",
"EmuJ_000733300","EmuJ_000733200","EmuJ_000730400","EmuJ_000730000","EmuJ_000710100","EmuJ_000651400","EmuJ_000521700","EmuJ_000521000",
"EmuJ_000366900","EmuJ_000357500","EmuJ_000351000","EmuJ_000297950","EmuJ_000052900","EmuJ_000020200","EmuJ_001085400",
"EmuJ_001085100","EmuJ_000249600","EmuJ_000981800","EmuJ_000331700","EmuJ_000331600","EmuJ_000938600","EmuJ_000917000",
"EmuJ_000910900","EmuJ_0001065400"))
length(hsp70)
#36
length(hsp70[hsp70 %in% names(tpm.PC.increase)])
#11
hsp70 <- hsp70[hsp70 %in% names(tpm.PC.increase)]


inicio <- rep(0,7324)
fin <- rep(0,7324)
names(inicio) <- names(tpm.PC.increase)
names(fin) <- names(tpm.PC.increase)

inicio[tetraspanin] <- 0
fin[tetraspanin] <- 1
inicio[tubulin] <- 1.1
fin[tubulin] <- 2.1
inicio[hsp70] <- 2.3
fin[hsp70] <- 3.3
inicio[histone] <- 3.4
fin[histone] <- 4.4
inicio[fabp] <- 4.5
fin[fabp] <- 5.5
inicio[agB] <- 5.6
fin[agB] <- 6.6

genes <- rep("WW",7324)
names(genes) <- names(tpm.PC.increase)
genes[agB] <- "agB"
genes[fabp] <- "FABPs"
genes[histone] <- "Histones"
genes[hsp70] <- "HSPs"
genes[tubulin] <- "Tubulins"
genes[tetraspanin] <- "Tetraspanins"

fig.df <- data.frame(gen=names(tpm.PC.increase), valor=tpm.PC.increase, tpm=1:7324, inicio=inicio, fin=fin, genes=genes)

dim(fig.df)
#[1] 7324    6


fig.p1 <- ggplot(fig.df, aes(x=log10(valor))) + geom_density(outline.type="full", fill="coral1", alpha=0.3) + geom_vline(xintercept = log10(50), size=0.9, linetype="dotted") + xlab("log10(TPM)") + theme_bw() 

fig.p2 <- ggplot(fig.df, aes(x=tpm, ymin=inicio, ymax=fin, color=genes)) + geom_linerange(size=1.52) + geom_vline(xintercept = (7324-577), size=1, linetype="dotted") + scale_y_discrete(labels = NULL, breaks = NULL) + theme_bw() + theme(legend.position="left") + xlab("gene order by increasing TPM")

library(cowplot)
pdf("fig3_SEPT21.pdf", height=5, width=15)
plot_grid(fig.p1,fig.p2, ncol=1, nrow=2, align="v", axis="lr")

dev.off()


dim(annots)
#[1] 10663     8
head(annots,2)
#           GeneID               Chromosome GeneStart GeneEnd Strand GeneName
#CDKD1.14 CDKD1.14 pathogen_EmW_scaffold_06    664987  666676      1         
#CLK2.7     CLK2.7 pathogen_EmW_scaffold_09   1248470 1274180     -1         
#                                                                 description
#CDKD1.14 Cyclin dependent kinase 1  [Source:UniProtKB/TrEMBL;Acc:A0A087VXA1]
#CLK2.7                                                                      
#                biotype
#CDKD1.14 protein_coding
#CLK2.7   protein_coding


#tpms.df.up1 <- tpms.df[genes.tpms.up1,]
dim(tpms.df.up1)
#[1] 7448    8
dim(annots[genes.tpms.up1,])
#[1] 7448    8
annots.up1 <- annots[genes.tpms.up1,]

tpms.out <- tpms.df.up1
tpms.out$Chromosome <- annots.up1$Chromosome
tpms.out$GeneStart <- annots.up1$GeneStart
tpms.out$GeneEnd <- annots.up1$GeneEnd
tpms.out$Strand <- annots.up1$Strand
tpms.out$GeneName <- annots.up1$GeneName
tpms.out$description <- annots.up1$description
tpms.out$biotype <- annots.up1$biotype

write.table(tpms.out,"genes_TPMs_above1.csv", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE, dec=",")

## update res.out with annotations
head(res.out)
dim(res.out)
#[1] 9197    9
names.res.out <- rownames(res.out)
annots.res <- annots[names.res.out,]
dim(annots.res)
#[1] 9197    8
res.out$Chromosome <- annots.res$Chromosome
res.out$GeneStart <- annots.res$GeneStart
res.out$GeneEnd <- annots.res$GeneEnd
res.out$Strand <- annots.res$Strand
res.out$GeneName <- annots.res$GeneName
res.out$description <- annots.res$description
res.out$biotype <- annots.res$biotype

write.table(res.out,"genes_log2FC_shrinkage-LFC1_INFO.csv", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE, dec=",")



##############################################
#### clusterProfiler - user's annotations ####
##############################################


## SpliceLeader genes Tsai&2013
head(sl.genes,2)
#         gene_id           ID
#1 EmuJ_000040600 SpliceLeader
#2 EmuJ_000041350 SpliceLeader
dim(sl.genes)
#[1] 1468    2


## Polycistron genes Tsai&2013
head(polycis.genes,2)
#  polycistron_id num.cistrons position.cistron        gene_id
#1              1            2                1 EmuJ_000040700
#2              1            2                2 EmuJ_000040600
#                     contig strand ini.position end.position PC_id          ID
#1 pathogen_EMU_contig_60709      -       134382       135263   PC1 polycistron
#2 pathogen_EMU_contig_60709      -       124150       130281   PC1 polycistron
dim(polycis.genes)
#[1] 643  10


## GO from WormBase Parasite biomaRt
head(listAttributes(mart),70)

mart <- useMart("parasite_mart", dataset="wbps_gene", host="https://parasite.wormbase.org", port=443)
genes.go <- getBM(mart = mart, 
               filters = c("species_id_1010", "wbps_gene_id"),
               value = list("ecmultprjeb122", names(tpm.PC.increase)),
               attributes = c("wbps_gene_id", "go_accession", "go_name_1006", "go_definition_1006", "go_linkage_type", "go_namespace_1003"))
head(genes.go)
go.BP <- subset(genes.go, go_namespace_1003=="biological_process")
go.MF <- subset(genes.go, go_namespace_1003=="molecular_function")
go.CC <- subset(genes.go, go_namespace_1003=="cellular_component")
dim(go.BP)
#[1] 5505    6
dim(go.MF)
#[1] 9135    6
dim(go.CC)
#[1] 5441    6


## 1# 577 genes with TPM>50 without treatment (PC)
tpms.PC.UP50

se.up50.GOBP <- enricher(gene=tpms.PC.UP50, universe=names(tpm.PC.increase), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=go.BP$go_accession, gene=go.BP$wbps_gene_id), TERM2NAME=data.frame(term=go.BP$go_accession, name=go.BP$go_name_1006))
as.data.frame(se.up50.GOBP)
as.data.frame(se.up50.GOBP)[,c(1:7,9)]
#                   ID
#GO:0006412 GO:0006412
#GO:0006457 GO:0006457
#GO:0015986 GO:0015986
#GO:0051603 GO:0051603
#GO:0006099 GO:0006099
#GO:0055114 GO:0055114
#GO:0022900 GO:0022900
#GO:0007017 GO:0007017
#GO:0006096 GO:0006096
#GO:0006414 GO:0006414
#GO:0099132 GO:0099132
#GO:0006418 GO:0006418
#GO:0007010 GO:0007010
#                                                          Description GeneRatio
#GO:0006412                                                translation    71/320
#GO:0006457                                            protein folding    19/320
#GO:0015986                     ATP synthesis coupled proton transport    10/320
#GO:0051603 proteolysis involved in cellular protein catabolic process    11/320
#GO:0006099                                   tricarboxylic acid cycle     9/320
#GO:0055114                                oxidation-reduction process    30/320
#GO:0022900                                   electron transport chain    10/320
#GO:0007017                                  microtubule-based process    14/320
#GO:0006096                                         glycolytic process     7/320
#GO:0006414                                   translational elongation     8/320
#GO:0099132      ATP hydrolysis coupled cation transmembrane transport     5/320
#GO:0006418                tRNA aminoacylation for protein translation     7/320
#GO:0007010                                  cytoskeleton organization     5/320
#            BgRatio       pvalue     p.adjust       qvalue Count
#GO:0006412 146/3277 1.168873e-35 7.247016e-34 5.659808e-34    71
#GO:0006457  37/3277 1.314541e-10 4.075076e-09 3.182572e-09    19
#GO:0015986  12/3277 3.821438e-09 7.897639e-08 6.167936e-08    10
#GO:0051603  15/3277 6.249757e-09 9.687124e-08 7.565496e-08    11
#GO:0006099  16/3277 4.460174e-06 5.530615e-05 4.319326e-05     9
#GO:0055114 153/3277 1.129049e-04 1.126164e-03 8.795168e-04    30
#GO:0022900  27/3277 1.271475e-04 1.126164e-03 8.795168e-04    10
#GO:0007017  53/3277 3.787484e-04 2.872514e-03 2.243389e-03    14
#GO:0006096  16/3277 4.169778e-04 2.872514e-03 2.243389e-03     7
#GO:0006414  26/3277 2.452643e-03 1.520639e-02 1.187595e-02     8
#GO:0099132  13/3277 5.722064e-03 3.225163e-02 2.518803e-02     5
#GO:0006418  29/3277 1.870108e-02 9.466850e-02 7.393465e-02     7
#GO:0007010  17/3277 1.984985e-02 9.466850e-02 7.393465e-02     5

write.table(as.data.frame(se.up50.GOBP), file="PC-TPMs-UP50_577genes_GO-BP.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")


pdf("PC-TPMs-UP50_577genes_GO-BP_dotplot.pdf", width=15, height=15)
dotplot(se.up50.GOBP, showCategory=15)
dev.off()
pdf("PC-TPMs-UP50_577genes_GO-BP_upsetplot.pdf", width=15, height=15)
upsetplot(se.up50.GOBP, showCategory=15)
dev.off()


se.up50.GOMF <- enricher(gene=tpms.PC.UP50, universe=names(tpm.PC.increase), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=go.MF$go_accession, gene=go.MF$wbps_gene_id), TERM2NAME=data.frame(term=go.MF$go_accession, name=go.MF$go_name_1006))
as.data.frame(se.up50.GOMF)
as.data.frame(se.up50.GOMF)[,c(1:7,9)]
#                   ID
#GO:0003735 GO:0003735
#GO:0051082 GO:0051082
#GO:0004298 GO:0004298
#GO:0004175 GO:0004175
#GO:0015078 GO:0015078
#GO:0008233 GO:0008233
#GO:0003723 GO:0003723
#GO:0016491 GO:0016491
#GO:0016616 GO:0016616
#GO:0000166 GO:0000166
#GO:0003746 GO:0003746
#GO:0003924 GO:0003924
#GO:0019843 GO:0019843
#GO:0003779 GO:0003779
#GO:0009055 GO:0009055
#                                                                                     Description
#GO:0003735                                                    structural constituent of ribosome
#GO:0051082                                                              unfolded protein binding
#GO:0004298                                                 threonine-type endopeptidase activity
#GO:0004175                                                                endopeptidase activity
#GO:0015078                                             proton transmembrane transporter activity
#GO:0008233                                                                    peptidase activity
#GO:0003723                                                                           RNA binding
#GO:0016491                                                               oxidoreductase activity
#GO:0016616 oxidoreductase activity, acting on the CH-OH group of donors, NAD or NADP as acceptor
#GO:0000166                                                                    nucleotide binding
#GO:0003746                                                translation elongation factor activity
#GO:0003924                                                                       GTPase activity
#GO:0019843                                                                          rRNA binding
#GO:0003779                                                                         actin binding
#GO:0009055                                                            electron transfer activity
#           GeneRatio  BgRatio       pvalue     p.adjust       qvalue Count
#GO:0003735    70/378 119/4150 4.802032e-44 3.361423e-42 2.527386e-42    70
#GO:0051082    15/378  25/4150 2.642844e-10 9.249955e-09 6.954853e-09    15
#GO:0004298    11/378  16/4150 8.958016e-09 2.090204e-07 1.571582e-07    11
#GO:0004175    10/378  14/4150 2.518763e-08 4.407836e-07 3.314162e-07    10
#GO:0015078     7/378  14/4150 9.603726e-05 1.344522e-03 1.010918e-03     7
#GO:0008233    19/378  86/4150 1.913049e-04 2.043608e-03 1.536547e-03    19
#GO:0003723    42/378 266/4150 2.043608e-04 2.043608e-03 1.536547e-03    42
#GO:0016491    20/378  97/4150 3.495993e-04 3.058994e-03 2.299995e-03    20
#GO:0016616     6/378  16/4150 1.992198e-03 1.549487e-02 1.165028e-02     6
#GO:0000166    56/378 426/4150 2.237341e-03 1.566139e-02 1.177548e-02    56
#GO:0003746     7/378  25/4150 5.560920e-03 3.538767e-02 2.660727e-02     7
#GO:0003924    15/378  85/4150 8.830171e-03 4.912055e-02 3.693274e-02    15
#GO:0019843     4/378  10/4150 9.122387e-03 4.912055e-02 3.693274e-02     4
#GO:0003779    10/378  49/4150 1.129786e-02 5.648928e-02 4.247314e-02    10
#GO:0009055     5/378  18/4150 1.930719e-02 9.010020e-02 6.774451e-02     5

write.table(as.data.frame(se.up50.GOMF), file="PCE-TPMs-UP50_577genes_GO-MF.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")


pdf("PC-TPMs-UP50_577genes_GO-MF_dotplot.pdf", width=15, height=15)
dotplot(se.up50.GOMF, showCategory=15)
dev.off()
pdf("PC-TPMs-UP50_577genes_GO-MF_upsetplot.pdf", width=15, height=15)
upsetplot(se.up50.GOMF, showCategory=15)
dev.off()


se.up50.GOCC <- enricher(gene=tpms.PC.UP50, universe=names(tpm.PC.increase), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=go.CC$go_accession, gene=go.CC$wbps_gene_id), TERM2NAME=data.frame(term=go.CC$go_accession, name=go.CC$go_name_1006))
as.data.frame(se.up50.GOCC)
as.data.frame(se.up50.GOCC)[,c(1:7,9)]
#                   ID                  Description GeneRatio  BgRatio
#GO:0005840 GO:0005840                     ribosome    77/314 153/2866
#GO:0005622 GO:0005622                intracellular    56/314 123/2866
#GO:0005737 GO:0005737                    cytoplasm    57/314 238/2866
#GO:0005839 GO:0005839      proteasome core complex    11/314  15/2866
#GO:0015935 GO:0015935      small ribosomal subunit     8/314  11/2866
#GO:0005743 GO:0005743 mitochondrial inner membrane    10/314  26/2866
#GO:0015934 GO:0015934      large ribosomal subunit     6/314  10/2866
#GO:0000502 GO:0000502           proteasome complex    13/314  43/2866
#GO:0005856 GO:0005856                 cytoskeleton    17/314  79/2866
#GO:0030286 GO:0030286               dynein complex    10/314  45/2866
#GO:0005783 GO:0005783        endoplasmic reticulum    11/314  53/2866
#                 pvalue     p.adjust       qvalue Count
#GO:0005840 1.028752e-36 2.880506e-35 1.840925e-35    77
#GO:0005622 1.193880e-23 1.671431e-22 1.068208e-22    56
#GO:0005737 1.588963e-09 1.483032e-08 9.478023e-09    57
#GO:0005839 2.113771e-08 1.479639e-07 9.456342e-08    11
#GO:0015935 2.342917e-06 1.312033e-05 8.385176e-06     8
#GO:0005743 2.302177e-04 9.459093e-04 6.045285e-04    10
#GO:0015934 2.364773e-04 9.459093e-04 6.045285e-04     6
#GO:0000502 4.408792e-04 1.543077e-03 9.861771e-04    13
#GO:0005856 4.215344e-03 1.311441e-02 8.381387e-03    17
#GO:0030286 2.081547e-02 5.828332e-02 3.724874e-02    10
#GO:0005783 2.557400e-02 6.509747e-02 4.160364e-02    11

write.table(as.data.frame(se.up50.GOCC), file="PC-TPMs-UP50_577genes_GO-CC.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")


pdf("PC-TPMs-UP50_577genes_GO-CC_dotplot.pdf", width=15, height=15)
dotplot(se.up50.GOCC, showCategory=15)
dev.off()
pdf("PC-TPMs-UP50_577genes_GO-CC_upsetplot.pdf", width=15, height=15)
upsetplot(se.up50.GOCC, showCategory=15)
dev.off()


se.up50.SLPC <- enricher(gene=tpms.PC.UP50, universe=names(tpm.PC.increase), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=1500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=c(sl.genes$ID,polycis.genes$ID), gene=c(sl.genes$gene_id,polycis.genes$gene_id)))
as.data.frame(se.up50.SLPC)
## NA

sum(tpms.PC.UP50 %in% sl.genes$gene_id)
#[1] 64
sum(names(tpm.PC.increase) %in% sl.genes$gene_id)
#[1] 1378
dim(sl.genes)
#[1] 1468    2

sum(tpms.PC.UP50 %in% polycis.genes$gene_id)
#[1] 26
sum(names(tpm.PC.increase) %in% polycis.genes$gene_id)
#[1] 616
dim(polycis.genes)
#[1] 643    2
subset(polycis.genes, gene_id %in% tpms.PC.UP50)
## only polycistron P242 has 2 in 3 genes among the 577 genes with TPM > 50



## 2# 7324 expressed genes (TPM>1) sorted
tpm.mean.PC.decreasing <- sort(tpm.mean.PC[tpm.mean.PC >=1], decreasing=TRUE)
length(tpm.mean.PC.decreasing)
#[1] 7324

se.up1.GOBP.gsea <- GSEA(tpm.mean.PC.decreasing, minGSSize=10, maxGSSize=500, pvalueCutoff=0.1, pAdjustMethod="BH", TERM2GENE=data.frame(term=go.BP$go_accession, gene=go.BP$wbps_gene_id), TERM2NAME=data.frame(term=go.BP$go_accession, name=go.BP$go_name_1006))
as.data.frame(se.up1.GOBP.gsea)[,1:10]
#                   ID
#GO:0006412 GO:0006412
#GO:0006457 GO:0006457
#GO:0006096 GO:0006096
#GO:0055114 GO:0055114
#GO:0006414 GO:0006414
#GO:0006099 GO:0006099
#GO:0051603 GO:0051603
#GO:0015986 GO:0015986
#                                                          Description setSize
#GO:0006412                                                translation     146
#GO:0006457                                            protein folding      37
#GO:0006096                                         glycolytic process      16
#GO:0055114                                oxidation-reduction process     153
#GO:0006414                                   translational elongation      26
#GO:0006099                                   tricarboxylic acid cycle      16
#GO:0051603 proteolysis involved in cellular protein catabolic process      15
#GO:0015986                     ATP synthesis coupled proton transport      12
#           enrichmentScore      NES       pvalue     p.adjust      qvalues rank
#GO:0006412       0.8325872 1.776688 1.727632e-09 1.744908e-07 1.582147e-07  350
#GO:0006457       0.8280383 1.705895 3.152565e-04 1.592045e-02 1.443543e-02  780
#GO:0006096       0.8892012 1.744349 1.387714e-03 3.685040e-02 3.341308e-02  331
#GO:0055114       0.6923741 1.477855 1.459422e-03 3.685040e-02 3.341308e-02 1308
#GO:0006414       0.8439613 1.720322 1.849836e-03 3.736668e-02 3.388120e-02  422
#GO:0006099       0.8565215 1.680241 5.296530e-03 7.834919e-02 7.104096e-02  820
#GO:0051603       0.8638436 1.686928 6.181089e-03 7.834919e-02 7.104096e-02  591
#GO:0015986       0.8689992 1.668150 6.205877e-03 7.834919e-02 7.104096e-02  524
#                             leading_edge
#GO:0006412  tags=47%, list=5%, signal=46%
#GO:0006457 tags=62%, list=11%, signal=56%
#GO:0006096  tags=44%, list=5%, signal=42%
#GO:0055114 tags=44%, list=18%, signal=37%
#GO:0006414  tags=31%, list=6%, signal=29%
#GO:0006099 tags=81%, list=11%, signal=72%
#GO:0051603  tags=87%, list=8%, signal=80%
#GO:0015986  tags=83%, list=7%, signal=77%


write.table(as.data.frame(se.up1.GOBP.gsea), file="PC-TPMs-UP1_7324genes_GO-BP-GSEA.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")

pdf("PC-TPMs-UP1_7324genes_GO-BP-GSEA_dotplot.pdf", width=15, height=15)
dotplot(se.up1.GOBP.gsea, showCategory=15)
dev.off()
pdf("PC-TPMs-UP1_7324genes_GO-BP-GSEA_upsetplot.pdf", width=15, height=15)
upsetplot(se.up1.GOBP.gsea, showCategory=15)
dev.off()
pdf("PC-TPMs-UP1_7324genes_GO-BP-GSEA_gseaplot2.pdf", width=18, height=7)
gseaplot2(se.up1.GOBP.gsea, geneSetID=c(1,3,4), title=se.up1.GOBP.gsea$Description[c(1,3,4)])
dev.off()


se.up1.GOMF.gsea <- GSEA(tpm.mean.PC.decreasing, minGSSize=10, maxGSSize=500, pvalueCutoff=0.1, pAdjustMethod="BH", TERM2GENE=data.frame(term=go.MF$go_accession, gene=go.MF$wbps_gene_id), TERM2NAME=data.frame(term=go.MF$go_accession, name=go.MF$go_name_1006))
as.data.frame(se.up1.GOMF.gsea)[,1:10]
#             ID                        Description setSize
#GO:0003735 GO:0003735 structural constituent of ribosome     119
#GO:0051082 GO:0051082           unfolded protein binding      25
#           enrichmentScore      NES       pvalue     p.adjust      qvalues rank
#GO:0003735       0.8853635 1.877288 0.0000000001 0.0000000103 9.157895e-09  350
#GO:0051082       0.8478228 1.734220 0.0019056090 0.0981388658 8.725684e-02  780
#                             leading_edge
#GO:0003735  tags=58%, list=5%, signal=56%
#GO:0051082 tags=72%, list=11%, signal=65%

write.table(as.data.frame(se.up1.GOMF.gsea), file="PC-TPMs-UP1_7324genes_GO-MF-GSEA.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")

pdf("PC-TPMs-UP1_7324genes_GO-MF-GSEA_dotplot.pdf", width=15, height=15)
dotplot(se.up1.GOMF.gsea, showCategory=15)
dev.off()
pdf("PC-TPMs-UP1_7324genes_GO-MF-GSEA_upsetplot.pdf", width=15, height=15)
upsetplot(se.up1.GOBP.gsea, showCategory=15)
dev.off()
pdf("PC-TPMs-UP1_7324genes_GO-MF-GSEA_gseaplot2.pdf", width=18, height=7)
gseaplot2(se.up1.GOMF.gsea, geneSetID=1:2, title=se.up1.GOMF.gsea$Description[1:2])
dev.off()


se.up1.GOCC.gsea <- GSEA(tpm.mean.PC.decreasing, minGSSize=10, maxGSSize=500, pvalueCutoff=0.1, pAdjustMethod="BH", TERM2GENE=data.frame(term=go.CC$go_accession, gene=go.CC$wbps_gene_id), TERM2NAME=data.frame(term=go.CC$go_accession, name=go.CC$go_name_1006))
as.data.frame(se.up1.GOCC.gsea)[,1:10]
#                   ID                  Description setSize enrichmentScore
#GO:0005840 GO:0005840                     ribosome     153       0.8623334
#GO:0005622 GO:0005622                intracellular     123       0.8336452
#GO:0015935 GO:0015935      small ribosomal subunit      11       0.9430583
#GO:0015934 GO:0015934      large ribosomal subunit      10       0.8992480
#GO:0005856 GO:0005856                 cytoskeleton      79       0.7426654
#GO:0005839 GO:0005839      proteasome core complex      15       0.8638436
#GO:0005737 GO:0005737                    cytoplasm     238       0.6408266
#GO:0005874 GO:0005874                  microtubule      54       0.6959676
#GO:0005743 GO:0005743 mitochondrial inner membrane      26       0.7523384
#GO:0000502 GO:0000502           proteasome complex      43       0.7012921
#GO:0030286 GO:0030286               dynein complex      45       0.6905099
#GO:0000786 GO:0000786                   nucleosome      17       0.7670175
#                NES       pvalue     p.adjust      qvalues rank
#GO:0005840 1.825024 1.000000e-10 3.400000e-09 2.315789e-09  350
#GO:0005622 1.761919 4.926994e-08 8.375889e-07 5.704940e-07  391
#GO:0015935 1.809365 1.418743e-04 1.607908e-03 1.095170e-03  203
#GO:0015934 1.706192 3.502397e-03 2.451973e-02 1.670075e-02  153
#GO:0005856 1.552120 3.605843e-03 2.451973e-02 1.670075e-02  899
#GO:0005839 1.705469 6.361997e-03 3.593071e-02 2.447293e-02  591
#GO:0005737 1.365278 7.397500e-03 3.593071e-02 2.447293e-02  915
#GO:0005874 1.449603 1.898102e-02 8.066933e-02 5.494505e-02  845
#GO:0005743 1.528917 2.500000e-02 8.645900e-02 5.888848e-02 1183
#GO:0000502 1.451265 2.597403e-02 8.645900e-02 5.888848e-02 1450
#GO:0030286 1.432533 2.797203e-02 8.645900e-02 5.888848e-02  942
#GO:0000786 1.532225 3.434343e-02 9.730640e-02 6.627680e-02 1072
#                             leading_edge
#GO:0005840  tags=50%, list=5%, signal=48%
#GO:0005622  tags=44%, list=5%, signal=42%
#GO:0015935  tags=73%, list=3%, signal=71%
#GO:0015934  tags=60%, list=2%, signal=59%
#GO:0005856 tags=32%, list=12%, signal=28%
#GO:0005839  tags=87%, list=8%, signal=80%
#GO:0005737 tags=36%, list=12%, signal=32%
#GO:0005874 tags=26%, list=12%, signal=23%
#GO:0005743 tags=65%, list=16%, signal=55%
#GO:0000502 tags=72%, list=20%, signal=58%
#GO:0030286 tags=31%, list=13%, signal=27%
#GO:0000786 tags=53%, list=15%, signal=45%

write.table(as.data.frame(se.up1.GOCC.gsea), file="PC-TPMs-UP1_7324genes_GO-CC-GSEA.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")

pdf("PC-TPMs-UP1_7324genes_GO-CC-GSEA_dotplot.pdf", width=15, height=15)
dotplot(se.up1.GOCC.gsea, showCategory=15)
dev.off()
pdf("PC-TPMs-UP1_7324genes_GO-CC-GSEA_upsetplot.pdf", width=15, height=15)
upsetplot(se.up1.GOCC.gsea, showCategory=15)
dev.off()
pdf("PC-TPMs-UP1_7324genes_GO-CC-GSEA_emapplot.pdf", width=15, height=15)
emapplot(se.up1.GOCC.gsea, showCategory=15)
dev.off()
pdf("PC-TPMs-UP1_7324genes_GO-CC-GSEA_gseaplot2.pdf", width=18, height=7)
gseaplot2(se.up1.GOCC.gsea, geneSetID=c(1,2,5), title=se.up1.GOCC.gsea$Description[c(1,2,5)])
dev.off()



se.up1.SLPC.gsea <- GSEA(tpm.mean.PC.decreasing, minGSSize=10, maxGSSize=1500, pvalueCutoff=0.1, pAdjustMethod="BH", TERM2GENE=data.frame(term=c(sl.genes$ID,polycis.genes$ID), gene=c(sl.genes$gene_id,polycis.genes$gene_id)))
as.data.frame(se.up1.SLPC.gsea)
#NA


## 3 # all 1645 genes DE con 0.05 alpha y abs(log2FC)>0

## all genes: genesDE (given log2FoldChange vector)
genesDE.LFC <- res05[genesDE,]$log2FoldChange
names(genesDE.LFC) <- genesDE
head(genesDE.LFC,2)
#EmuJ_000354850 EmuJ_000009000 
#      0.971126       1.015039

de.1645.GOBP <- enricher(gene=genesDE, universe=rownames(res05), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=go.BP$go_accession, gene=go.BP$wbps_gene_id), TERM2NAME=data.frame(term=go.BP$go_accession, name=go.BP$go_name_1006))
as.data.frame(de.1645.GOBP)
as.data.frame(de.1645.GOBP)[,c(1:7,9)]
#                   ID                    Description GeneRatio  BgRatio
#GO:0005975 GO:0005975 carbohydrate metabolic process    23/682  50/3277
#GO:0006096 GO:0006096             glycolytic process    10/682  16/3277
#GO:0055085 GO:0055085        transmembrane transport    49/682 150/3277
#GO:0007017 GO:0007017      microtubule-based process    20/682  53/3277
#                 pvalue    p.adjust      qvalue Count
#GO:0005975 5.043825e-05 0.004741195 0.004353617    23
#GO:0006096 3.364601e-04 0.010864395 0.009976264    10
#GO:0055085 3.467360e-04 0.010864395 0.009976264    49
#GO:0007017 3.267619e-03 0.076789053 0.070511785    20

#carbohydrate metabolic process
summary(res.out[strsplit(de.1645.GOBP$geneID,"/")[[1]],]$log2FoldChange)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1.0954 -0.5856 -0.3747 -0.3648 -0.3104  0.5556 
summary(res.out[strsplit(de.1645.GOBP$geneID,"/")[[1]],]$apeglm)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1.0639 -0.5588 -0.3363 -0.3438 -0.2833  0.5128 
#glycolytic process
summary(res.out[strsplit(de.1645.GOBP$geneID,"/")[[2]],]$apeglm)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-0.6231 -0.5176 -0.4530 -0.3639 -0.3918  0.5128
#microtubule-based process
summary(res.out[strsplit(de.1645.GOBP$geneID,"/")[[4]],]$apeglm)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1.4049 -0.9269 -0.5202 -0.4827 -0.2271  0.9409
## todos bajos excepto EmuJ_000967000 con 0.9409 de apeglm
## EmuJ_000967000: Dynein light chain type 1 2  [Source:UniProtKB/TrEMBL;Acc:A0A068YI72]
#transmembrane transport
summary(res.out[strsplit(de.1645.GOBP$geneID,"/")[[3]],]$apeglm)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1.3986 -0.5424 -0.3431 -0.2131  0.2927  1.9158
#EmuJ_000971400 0.9867939; Phosphate transporter  [Source:UniProtKB/TrEMBL;Acc:A0A068YFD6]
#EmuJ_000931700 1.9158447; ATP binding cassette subfamily B MDR:TAP  [Source:UniProtKB/TrEMBL;Acc:A0A068YEA9]
#EmuJ_000901300 0.5416329; ATP binding cassette subfamily B MDR TAP  [Source:UniProtKB/TrEMBL;Acc:A0A068Y998]
#EmuJ_000714000 0.6277917; Solute carrier family 5  [Source:UniProtKB/TrEMBL;Acc:A0A068Y8M1]
#EmuJ_001165800 0.5793365; Sodium:potassium:calcium exchanger  [Source:UniProtKB/TrEMBL;Acc:A0A068XTP1]



write.table(as.data.frame(de.1645.GOBP), file="DE_1645genes_GO-BP.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")


pdf("DE_1645genes_GO-BP_dotplot.pdf", width=15, height=15)
dotplot(de.1645.GOBP, showCategory=15)
dev.off()
pdf("DE_1645genes_GO-BP_upsetplot.pdf", width=15, height=15)
upsetplot(de.1645.GOBP, showCategory=15)
dev.off()
pdf("DE_1645genes_GO-BP_heatplot.pdf", width=15, height=15)
heatplot(de.1645.GOBP, foldChange=genesDE.LFC)
dev.off()
pdf("DE_1645genes_GO-BP_cnetplot.pdf", width=15, height=15)
cnetplot(de.1645.GOBP, categorySize="pvalue")
dev.off()


de.1645.GOMF <- enricher(gene=genesDE, universe=rownames(res05), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=go.MF$go_accession, gene=go.MF$wbps_gene_id), TERM2NAME=data.frame(term=go.MF$go_accession, name=go.MF$go_name_1006))
as.data.frame(de.1645.GOMF)
as.data.frame(de.1645.GOMF)[,c(1:7,9)]
#                   ID                                          Description
#GO:0004553 GO:0004553 hydrolase activity, hydrolyzing O-glycosyl compounds
#GO:0005509 GO:0005509                                  calcium ion binding
#           GeneRatio  BgRatio       pvalue   p.adjust     qvalue Count
#GO:0004553     8/853  12/4150 0.0006924267 0.05781162 0.05422664     8
#GO:0005509    43/853 135/4150 0.0011447845 0.05781162 0.05422664    43


write.table(as.data.frame(de.1645.GOMF), file="DE_1645genes_GO-MF.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")

pdf("DE_1645genes_GO-MF_heatplot.pdf", width=15, height=15)
heatplot(de.1645.GOMF, foldChange=genesDE.LFC)
dev.off()


de.1645.GOCC <- enricher(gene=genesDE, universe=rownames(res05), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=go.CC$go_accession, gene=go.CC$wbps_gene_id), TERM2NAME=data.frame(term=go.CC$go_accession, name=go.CC$go_name_1006))
as.data.frame(de.1645.GOCC)
as.data.frame(de.1645.GOCC)[,c(1:7,9)]
#                   ID                    Description GeneRatio BgRatio
#GO:0000786 GO:0000786                     nucleosome    11/611 17/2866
#GO:0005694 GO:0005694                     chromosome    12/611 26/2866
#GO:0005875 GO:0005875 microtubule associated complex     6/611 10/2866
#                 pvalue    p.adjust      qvalue Count
#GO:0000786 0.0001326246 0.004376613 0.004048542    11
#GO:0005694 0.0039438127 0.065072910 0.060195036    12
#GO:0005875 0.0087380013 0.096118014 0.088912995     6


write.table(as.data.frame(de.1645.GOCC), file="DE_1645genes_GO-CC.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")

pdf("DE_1645genes_GO-CC_heatplot.pdf", width=15, height=15)
heatplot(de.1645.GOCC, foldChange=genesDE.LFC)
dev.off()
pdf("DE_1645genes_GO-CC_cnetplot.pdf", width=15, height=15)
cnetplot(de.1645.GOCC, categorySize="pvalue")
dev.off()
pdf("DE_1645genes_GO-CC_upsetplot.pdf", width=15, height=15)
upsetplot(de.1645.GOCC)
dev.off()


de.1645.SLPC <- enricher(gene=genesDE, universe=rownames(res05), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=1500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=c(sl.genes$ID,polycis.genes$ID), gene=c(sl.genes$gene_id,polycis.genes$gene_id)))
as.data.frame(de.1645.SLPC)
## NA



## 4 # 174 genes with DE (0.05 alpha y abs(log2FC)>1 considering lfcshrinkage apeglm)

## genes.de.174
genes.174.LFC <- res.apeglm[genes.de.174,]$log2FoldChange
names(genes.174.LFC) <- genes.de.174
head(genes.174.LFC,2)
#EmuJ_002133600 EmuJ_000601800 
#      1.849363      -1.082530

de.174.GOBP <- enricher(gene=genes.de.174, universe=rownames(res.apeglm), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=go.BP$go_accession, gene=go.BP$wbps_gene_id), TERM2NAME=data.frame(term=go.BP$go_accession, name=go.BP$go_name_1006))
as.data.frame(de.174.GOBP)
as.data.frame(de.174.GOBP)[,c(1:7,9)]
#                   ID Description GeneRatio  BgRatio      pvalue   p.adjust
#GO:0006508 GO:0006508 proteolysis      6/34 143/3277 0.003037895 0.08202317
#               qvalue Count
#GO:0006508 0.07674683     6
summary(res.out[strsplit(de.174.GOBP$geneID,"/")[[1]],]$log2FoldChange)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1.7289 -1.6950 -1.6351 -0.7012  0.5205  1.2676
summary(res.out[strsplit(de.174.GOBP$geneID,"/")[[1]],]$apeglm)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-1.5610 -1.5401 -1.5240 -0.6191  0.5192  1.2398
res.out[strsplit(de.174.GOBP$geneID,"/")[[1]],]
#               baseMean log2FoldChange     lfcSE      stat       pvalue
#EmuJ_000654100 1839.131      -1.728947 0.3644630 -4.743822 2.097236e-06
#EmuJ_000824000  906.345      -1.592477 0.1997663 -7.971697 1.565096e-15
#EmuJ_001193200 4315.651       1.224853 0.1339986  9.140792 6.199638e-20
#EmuJ_001193100 4684.751       1.267634 0.1338799  9.468445 2.840412e-21
#EmuJ_000654500 1617.361      -1.677700 0.3566281 -4.704340 2.546879e-06
#EmuJ_000654600 2142.429      -1.700789 0.3616353 -4.703051 2.563027e-06
#                       padj    apeglm apeglm.svalue      LFC1
#EmuJ_000654100 5.267694e-05 -1.561009   0.013776899 1.0000000
#EmuJ_000824000 2.558790e-13 -1.542019   0.000602437 0.5223123
#EmuJ_001193200 1.990970e-17  1.196896   0.015266336 1.0000000
#EmuJ_001193100 1.110477e-18  1.239782   0.007015271 1.0000000
#EmuJ_000654500 6.274394e-05 -1.513771   0.021783595 1.0000000
#EmuJ_000654600 6.296923e-05 -1.534317   0.021043639 1.0000000
#                             Chromosome GeneStart  GeneEnd Strand    GeneName
#EmuJ_000654100 pathogen_EmW_scaffold_05   9182745  9185045      1            
#EmuJ_000824000 pathogen_EmW_scaffold_01   8015047  8017040      1 serpin2Emub
#EmuJ_001193200 pathogen_EmW_scaffold_01  19495210 19496442      1            
#EmuJ_001193100 pathogen_EmW_scaffold_01  19488520 19489752      1            
#EmuJ_000654500 pathogen_EmW_scaffold_05   9210378  9212678     -1      EmCLP1
#EmuJ_000654600 pathogen_EmW_scaffold_05   9215950  9218250      1      EmCLP1
#                                                                                                           description
#EmuJ_000654100                                Cathepsin l cysteine peptidase  [Source:UniProtKB/TrEMBL;Acc:A0A0S4MN34]
#EmuJ_000824000 Estrogen regulated protein EP45; Serine protease inhibitor 2b  [Source:UniProtKB/TrEMBL;Acc:A0A068Y801]
#EmuJ_001193200                                     Serine protease inhibitor  [Source:UniProtKB/TrEMBL;Acc:A0A068XUR3]
#EmuJ_001193100                                     Serine protease inhibitor  [Source:UniProtKB/TrEMBL;Acc:A0A068Y1K4]
#EmuJ_000654500                    Cathepsin L-like proteinase; Cysteine protease  [Source:UniProtKB/TrEMBL;Acc:Q0WYD8]
#EmuJ_000654600                    Cathepsin L-like proteinase; Cysteine protease  [Source:UniProtKB/TrEMBL;Acc:Q0WYD8]
#                      biotype
#EmuJ_000654100 protein_coding
#EmuJ_000824000 protein_coding
#EmuJ_001193200 protein_coding
#EmuJ_001193100 protein_coding
#EmuJ_000654500 protein_coding
#EmuJ_000654600 protein_coding


write.table(as.data.frame(de.174.GOBP), file="DE_174genes_GO-BP.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")

pdf("DE_174genes_GO-BP_heatplot.pdf", width=15, height=15)
heatplot(de.174.GOBP, foldChange=genes.174.LFC)
dev.off()


de.174.GOMF <- enricher(gene=genes.de.174, universe=rownames(res.apeglm), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=go.MF$go_accession, gene=go.MF$wbps_gene_id), TERM2NAME=data.frame(term=go.MF$go_accession, name=go.MF$go_name_1006))
as.data.frame(de.174.GOMF)
as.data.frame(de.174.GOMF)[,c(1:7,9)]
#                   ID                      Description GeneRatio BgRatio
#GO:0008233 GO:0008233               peptidase activity      6/45 86/4150
#GO:0008234 GO:0008234 cysteine-type peptidase activity      3/45 26/4150
#                 pvalue    p.adjust      qvalue Count
#GO:0008233 0.0002831342 0.008777159 0.008643043     6
#GO:0008234 0.0026024878 0.040338561 0.039722182     3

res.out[strsplit(de.174.GOMF$geneID,"/")[[1]],]
#               baseMean log2FoldChange     lfcSE      stat       pvalue
#EmuJ_000654100 1839.131      -1.728947 0.3644630 -4.743822 2.097236e-06
#EmuJ_000824000  906.345      -1.592477 0.1997663 -7.971697 1.565096e-15
#EmuJ_001193200 4315.651       1.224853 0.1339986  9.140792 6.199638e-20
#EmuJ_001193100 4684.751       1.267634 0.1338799  9.468445 2.840412e-21
#EmuJ_000654500 1617.361      -1.677700 0.3566281 -4.704340 2.546879e-06
#EmuJ_000654600 2142.429      -1.700789 0.3616353 -4.703051 2.563027e-06
#                       padj    apeglm apeglm.svalue      LFC1
#EmuJ_000654100 5.267694e-05 -1.561009   0.013776899 1.0000000
#EmuJ_000824000 2.558790e-13 -1.542019   0.000602437 0.5223123
#EmuJ_001193200 1.990970e-17  1.196896   0.015266336 1.0000000
#EmuJ_001193100 1.110477e-18  1.239782   0.007015271 1.0000000
#EmuJ_000654500 6.274394e-05 -1.513771   0.021783595 1.0000000
#EmuJ_000654600 6.296923e-05 -1.534317   0.021043639 1.0000000
#                             Chromosome GeneStart  GeneEnd Strand    GeneName
#EmuJ_000654100 pathogen_EmW_scaffold_05   9182745  9185045      1            
#EmuJ_000824000 pathogen_EmW_scaffold_01   8015047  8017040      1 serpin2Emub
#EmuJ_001193200 pathogen_EmW_scaffold_01  19495210 19496442      1            
#EmuJ_001193100 pathogen_EmW_scaffold_01  19488520 19489752      1            
#EmuJ_000654500 pathogen_EmW_scaffold_05   9210378  9212678     -1      EmCLP1
#EmuJ_000654600 pathogen_EmW_scaffold_05   9215950  9218250      1      EmCLP1
#                                                                                                           description
#EmuJ_000654100                                Cathepsin l cysteine peptidase  [Source:UniProtKB/TrEMBL;Acc:A0A0S4MN34]
#EmuJ_000824000 Estrogen regulated protein EP45; Serine protease inhibitor 2b  [Source:UniProtKB/TrEMBL;Acc:A0A068Y801]
#EmuJ_001193200                                     Serine protease inhibitor  [Source:UniProtKB/TrEMBL;Acc:A0A068XUR3]
#EmuJ_001193100                                     Serine protease inhibitor  [Source:UniProtKB/TrEMBL;Acc:A0A068Y1K4]
#EmuJ_000654500                    Cathepsin L-like proteinase; Cysteine protease  [Source:UniProtKB/TrEMBL;Acc:Q0WYD8]
#EmuJ_000654600                    Cathepsin L-like proteinase; Cysteine protease  [Source:UniProtKB/TrEMBL;Acc:Q0WYD8]
#                      biotype
#EmuJ_000654100 protein_coding
#EmuJ_000824000 protein_coding
#EmuJ_001193200 protein_coding
#EmuJ_001193100 protein_coding
#EmuJ_000654500 protein_coding
#EmuJ_000654600 protein_coding


write.table(as.data.frame(de.174.GOMF), file="DE_174genes_GO-MF.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")

pdf("DE_174genes_GO-MF_heatplot.pdf", width=15, height=15)
heatplot(de.174.GOMF, foldChange=genes.174.LFC)
dev.off()
pdf("DE_174genes_GO-MF_cnetplot.pdf", width=15, height=15)
cnetplot(de.174.GOMF, categorySize="pvalue", circular=TRUE, colorEdge=TRUE)
dev.off()


de.174.GOCC <- enricher(gene=genes.de.174, universe=rownames(res.apeglm), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=go.CC$go_accession, gene=go.CC$wbps_gene_id), TERM2NAME=data.frame(term=go.CC$go_accession, name=go.CC$go_name_1006))
as.data.frame(de.174.GOCC)
as.data.frame(de.174.GOCC)[,c(1:7,9)]
#                   ID                    Description GeneRatio BgRatio
#GO:0005856 GO:0005856                   cytoskeleton      5/46 79/2866
#GO:0005875 GO:0005875 microtubule associated complex      2/46 10/2866
#                pvalue   p.adjust     qvalue Count
#GO:0005856 0.007928098 0.04180276 0.03300218     5
#GO:0005875 0.010450689 0.04180276 0.03300218     2

res.out[strsplit(de.174.GOCC$geneID,"/")[[1]],]
#                baseMean log2FoldChange     lfcSE      stat       pvalue
#EmuJ_000940900 2401.4784      -1.285896 0.1758679 -7.311712 2.637603e-13
#EmuJ_000590100  472.7547      -1.087497 0.2056704 -5.287572 1.239505e-07
#EmuJ_000497500  989.8206      -1.472853 0.2024592 -7.274810 3.469082e-13
#EmuJ_000773400 2252.3106       1.177141 0.1551735  7.585962 3.300281e-14
#EmuJ_000569000  118.4031      -1.502558 0.2668412 -5.630906 1.792651e-08
#                       padj    apeglm apeglm.svalue LFC1
#EmuJ_000940900 2.664868e-11 -1.237873   0.023286840    1
#EmuJ_000590100 4.286781e-06 -1.009468   0.147399277    1
#EmuJ_000497500 3.390650e-11 -1.416795   0.003700508    1
#EmuJ_000773400 4.179736e-12  1.136731   0.047866549    1
#EmuJ_000569000 7.293898e-07 -1.404864   0.013008532    1
#                             Chromosome GeneStart  GeneEnd Strand GeneName
#EmuJ_000940900 pathogen_EmW_scaffold_02   1720332  1720736     -1         
#EmuJ_000590100 pathogen_EmW_scaffold_05   1998970  1999313      1         
#EmuJ_000497500 pathogen_EmW_scaffold_03   1715551  1731718     -1         
#EmuJ_000773400 pathogen_EmW_scaffold_01   4018001  4042513     -1         
#EmuJ_000569000 pathogen_EmW_scaffold_01  15924662 15926066      1         
#                                                                               description
#EmuJ_000940900                Dynein light chain  [Source:UniProtKB/TrEMBL;Acc:A0A068YA92]
#EmuJ_000590100                Dynein light chain  [Source:UniProtKB/TrEMBL;Acc:A0A068Y1W1]
#EmuJ_000497500 E3 ubiquitin protein ligase MYLIP  [Source:UniProtKB/TrEMBL;Acc:A0A068XZ12]
#EmuJ_000773400                      Pleckstrin y  [Source:UniProtKB/TrEMBL;Acc:A0A068Y5U2]
#EmuJ_000569000                Tubulin beta chain  [Source:UniProtKB/TrEMBL;Acc:A0A068Y7X2]
#                      biotype
#EmuJ_000940900 protein_coding
#EmuJ_000590100 protein_coding
#EmuJ_000497500 protein_coding
#EmuJ_000773400 protein_coding
#EmuJ_000569000 protein_coding
res.out[strsplit(de.174.GOCC$geneID,"/")[[2]],]
#                baseMean log2FoldChange     lfcSE      stat       pvalue
#EmuJ_000940900 2401.4784      -1.285896 0.1758679 -7.311712 2.637603e-13
#EmuJ_000590100  472.7547      -1.087497 0.2056704 -5.287572 1.239505e-07
#                       padj    apeglm apeglm.svalue LFC1
#EmuJ_000940900 2.664868e-11 -1.237873    0.02328684    1
#EmuJ_000590100 4.286781e-06 -1.009468    0.14739928    1
#                             Chromosome GeneStart GeneEnd Strand GeneName
#EmuJ_000940900 pathogen_EmW_scaffold_02   1720332 1720736     -1         
#EmuJ_000590100 pathogen_EmW_scaffold_05   1998970 1999313      1         
#                                                                description
#EmuJ_000940900 Dynein light chain  [Source:UniProtKB/TrEMBL;Acc:A0A068YA92]
#EmuJ_000590100 Dynein light chain  [Source:UniProtKB/TrEMBL;Acc:A0A068Y1W1]
#                      biotype
#EmuJ_000940900 protein_coding
#EmuJ_000590100 protein_coding


write.table(as.data.frame(de.174.GOCC), file="DE_174genes_GO-CC.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")

pdf("DE_174genes_GO-CC_heatplot.pdf", width=15, height=15)
heatplot(de.174.GOCC, foldChange=genes.174.LFC)
dev.off()
pdf("DE_174genes_GO-CC_cnetplot.pdf", width=15, height=15)
cnetplot(de.174.GOCC, categorySize="pvalue", circular=TRUE, colorEdge=TRUE)
dev.off()


de.174.SLPC <- enricher(gene=genes.de.174, universe=rownames(res.apeglm), pvalueCutoff=0.1, pAdjustMethod="BH", minGSSize=10, maxGSSize=1500, qvalueCutoff=0.2, TERM2GENE=data.frame(term=c(sl.genes$ID,polycis.genes$ID), gene=c(sl.genes$gene_id,polycis.genes$gene_id)))
as.data.frame(de.174.SLPC)
## NA

sum(genes.de.174 %in% sl.genes$gene_id)
#[1] 7
sum(rownames(res05) %in% sl.genes$gene_id)
#[1] 1423
dim(sl.genes)
#[1] 1468    2

sum(genes.de.174 %in% polycis.genes$gene_id)
#[1] 2
sum(rownames(res05) %in% polycis.genes$gene_id)
#[1] 625
dim(polycis.genes)
#[1] 643    2
subset(polycis.genes, gene_id %in% genes.de.174)
#    polycistron_id num.cistrons position.cistron        gene_id
#193            116            2                2 EmuJ_000448200
#341            206            2                1 EmuJ_000724500
#                          contig strand ini.position end.position PC_id
#193 pathogen_EMU_scaffold_007614      -      2452949      2469957 PC116
#341 pathogen_EMU_scaffold_007765      -        84228        84709 PC206
#             ID
#193 polycistron
#341 polycistron


## 5 # all genes for GSEA but using apeglm lfcshrinkage (alpha 0.05)

## filter NAs
res.apeglm.filt <- subset(res.apeglm, !is.na(padj)) 
dim(res.apeglm.filt)
#[1] 8992    6

head(sort(res.apeglm.filt$log2FoldChange, decreasing=TRUE),2)
#EmuJ_001102600 EmuJ_002194600 
#      12.09064        7.49327 
tail(sort(res.apeglm.filt$log2FoldChange, decreasing=TRUE),2)
#EmuJ_000451500 EmuJ_000274300 
#     -2.173735     -10.025932 
res.apeglm.8992.filt.FC <- sort(res.apeglm.filt$log2FoldChange, decreasing=TRUE)


de.8992.ap.GOBP.gsea <- GSEA(res.apeglm.8992.filt.FC, minGSSize=10, maxGSSize=500, pvalueCutoff=0.1, pAdjustMethod="BH", TERM2GENE=data.frame(term=go.BP$go_accession, gene=go.BP$wbps_gene_id), TERM2NAME=data.frame(term=go.BP$go_accession, name=go.BP$go_name_1006))
as.data.frame(de.8992.ap.GOBP.gsea)
as.data.frame(de.8992.ap.GOBP.gsea)[,1:10]
#                 ID                                  Description setSize
#GO:0006412 GO:0006412                                  translation     146
#GO:0005975 GO:0005975               carbohydrate metabolic process      50
#GO:0008152 GO:0008152                            metabolic process      18
#GO:0055114 GO:0055114                  oxidation-reduction process     153
#GO:0007017 GO:0007017                    microtubule-based process      53
#GO:0007186 GO:0007186 G protein-coupled receptor signaling pathway      24
#           enrichmentScore       NES       pvalue   p.adjust    qvalues rank
#GO:0006412      -0.5295100 -1.862064 0.0002124376 0.02188107 0.01990205 2241
#GO:0005975      -0.6134770 -1.843807 0.0006302509 0.03245792 0.02952228  990
#GO:0008152      -0.7587250 -1.870160 0.0011908577 0.04088611 0.03718819  990
#GO:0055114      -0.4872432 -1.717188 0.0019290275 0.04967246 0.04517985 1626
#GO:0007017      -0.6215259 -1.873226 0.0024545667 0.05036999 0.04581430 1331
#GO:0007186      -0.6871059 -1.795396 0.0029341743 0.05036999 0.04581430  940
#                             leading_edge
#GO:0006412 tags=58%, list=25%, signal=44%
#GO:0005975 tags=38%, list=11%, signal=34%
#GO:0008152 tags=44%, list=11%, signal=40%
#GO:0055114 tags=40%, list=18%, signal=33%
#GO:0007017 tags=42%, list=15%, signal=36%
#GO:0007186 tags=46%, list=10%, signal=41%


summary(res.out[strsplit(de.1645.GOBP$geneID,"/")[[3]],]$apeglm)



write.table(as.data.frame(de.8992.ap.GOBP.gsea), file="DE-log2FC_8992genes_GOBP-GSEA.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")


pdf("DE-log2FC_8992genes_GOBP-GSEA_dotplot.pdf", width=15, height=15)
dotplot(de.8992.ap.GOBP.gsea, showCategory=15)
dev.off()
pdf("DE-log2FC_8992genes_GOBP-GSEA_upsetplot.pdf", width=15, height=15)
upsetplot(de.8992.ap.GOBP.gsea, showCategory=15)
dev.off()
pdf("DE-log2FC_8992genes_GOBP-GSEA_emapplot.pdf", width=15, height=15)
emapplot(de.8992.ap.GOBP.gsea, showCategory=15)
dev.off()
pdf("DE-log2FC_8992genes_GOBP-GSEA_gseaplot2.pdf", width=18, height=7)
gseaplot2(de.8992.ap.GOBP.gsea, geneSetID=c(1,2,5), title=se.up1.GOCC.gsea$Description[c(1,2,5)])
dev.off()
pdf("DE-log2FC_8992genes_GOBP-GSEA_heatplot.pdf", width=18, height=15)
heatplot(de.8992.ap.GOBP.gsea, foldChange=res.apeglm.8992.filt.FC)
dev.off()



de.8992.ap.GOMF.gsea <- GSEA(res.apeglm.8992.filt.FC, minGSSize=10, maxGSSize=500, pvalueCutoff=0.1, pAdjustMethod="BH", TERM2GENE=data.frame(term=go.MF$go_accession, gene=go.MF$wbps_gene_id), TERM2NAME=data.frame(term=go.MF$go_accession, name=go.MF$go_name_1006))
as.data.frame(de.8992.ap.GOMF.gsea)
as.data.frame(de.8992.ap.GOMF.gsea)[,1:10]
#                   ID                                          Description
#GO:0003735 GO:0003735                   structural constituent of ribosome
#GO:0004553 GO:0004553 hydrolase activity, hydrolyzing O-glycosyl compounds
#GO:0016798 GO:0016798         hydrolase activity, acting on glycosyl bonds
#GO:0003824 GO:0003824                                   catalytic activity
#GO:0008234 GO:0008234                     cysteine-type peptidase activity
#           setSize enrichmentScore       NES       pvalue     p.adjust
#GO:0003735     119      -0.6081060 -2.077091 6.840079e-06 0.0007113682
#GO:0004553      12      -0.8283505 -1.857516 7.940247e-04 0.0412892841
#GO:0016798      17      -0.7631618 -1.878072 1.269321e-03 0.0440031316
#GO:0003824     183      -0.4501329 -1.642543 2.804224e-03 0.0721744480
#GO:0008234      26      -0.6787377 -1.824896 3.469925e-03 0.0721744480
#                qvalues rank                   leading_edge
#GO:0003735 0.0006408074 2241 tags=69%, list=25%, signal=52%
#GO:0004553 0.0371937883  990 tags=67%, list=11%, signal=59%
#GO:0016798 0.0396384485  990 tags=47%, list=11%, signal=42%
#GO:0003824 0.0650154441 1623 tags=31%, list=18%, signal=26%
#GO:0008234 0.0650154441  338  tags=23%, list=4%, signal=22%

write.table(as.data.frame(de.8992.ap.GOMF.gsea), file="DE-log2FC_8992genes_GOMF-GSEA.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")


pdf("DE-log2FC_8992genes_GOMF-GSEA_dotplot.pdf", width=15, height=15)
dotplot(de.8992.ap.GOMF.gsea, showCategory=15)
dev.off()
pdf("DE-log2FC_8992genes_GOMF-GSEA_upsetplot.pdf", width=15, height=15)
upsetplot(de.8992.ap.GOMF.gsea, showCategory=15)
dev.off()
pdf("DE-log2FC_8992genes_GOMF-GSEA_emapplot.pdf", width=15, height=15)
emapplot(de.8992.ap.GOMF.gsea, showCategory=15)
dev.off()
pdf("DE-log2FC_8992genes_GOMF-GSEA_gseaplot2.pdf", width=18, height=7)
gseaplot2(de.8992.ap.GOMF.gsea, geneSetID=c(1,2,5), title=se.up1.GOMF.gsea$Description[c(1,2,5)])
dev.off()
pdf("DE-log2FC_8992genes_GOMF-GSEA_heatplot.pdf", width=18, height=15)
heatplot(de.8992.ap.GOMF.gsea, foldChange=res.apeglm.8992.filt.FC)
dev.off()



de.8992.ap.GOCC.gsea <- GSEA(res.apeglm.8992.filt.FC, minGSSize=10, maxGSSize=500, pvalueCutoff=0.1, pAdjustMethod="BH", TERM2GENE=data.frame(term=go.CC$go_accession, gene=go.CC$wbps_gene_id), TERM2NAME=data.frame(term=go.CC$go_accession, name=go.CC$go_name_1006))
as.data.frame(de.8992.ap.GOCC.gsea)
as.data.frame(de.8992.ap.GOCC.gsea)[,1:10]
#                   ID                    Description setSize enrichmentScore
#GO:0005840 GO:0005840                       ribosome     153      -0.5578529
#GO:0005622 GO:0005622                  intracellular     123      -0.4793106
#GO:0030286 GO:0030286                 dynein complex      45      -0.5838638
#GO:0005875 GO:0005875 microtubule associated complex      10      -0.7761922
#GO:0005886 GO:0005886                plasma membrane      37      -0.6219783
#GO:0005874 GO:0005874                    microtubule      54      -0.5659066
#GO:0005856 GO:0005856                   cytoskeleton      79      -0.4770079
#GO:0000139 GO:0000139                 Golgi membrane      25      -0.6143528
#                 NES       pvalue     p.adjust      qvalues rank
#GO:0005840 -1.969748 1.176263e-05 0.0003881668 0.0002847795 2241
#GO:0005622 -1.635998 2.438814e-03 0.0402404242 0.0295224803 2241
#GO:0030286 -1.701588 5.389299e-03 0.0578905375 0.0424715267 1259
#GO:0005875 -1.647927 8.249952e-03 0.0578905375 0.0424715267 1217
#GO:0005886 -1.744646 1.016221e-02 0.0578905375 0.0424715267 1889
#GO:0005874 -1.706533 1.194501e-02 0.0578905375 0.0424715267 1592
#GO:0005856 -1.524364 1.227981e-02 0.0578905375 0.0424715267 1259
#GO:0000139 -1.595943 1.659473e-02 0.0684532431 0.0502208801  474
#                             leading_edge
#GO:0005840 tags=58%, list=25%, signal=44%
#GO:0005622 tags=51%, list=25%, signal=39%
#GO:0030286 tags=40%, list=14%, signal=35%
#GO:0005875 tags=80%, list=14%, signal=69%
#GO:0005886 tags=54%, list=21%, signal=43%
#GO:0005874 tags=39%, list=18%, signal=32%
#GO:0005856 tags=27%, list=14%, signal=23%
#GO:0000139  tags=20%, list=5%, signal=19%


write.table(as.data.frame(de.8992.ap.GOCC.gsea), file="DE-log2FC_8992genes_GOCC-GSEA.csv", quote=FALSE, sep="\t", row.names=FALSE, dec=",")


pdf("DE-log2FC_8992genes_GOCC-GSEA_dotplot.pdf", width=15, height=15)
dotplot(de.8992.ap.GOCC.gsea, showCategory=15)
dev.off()
pdf("DE-log2FC_8992genes_GOCC-GSEA_upsetplot.pdf", width=15, height=15)
upsetplot(de.8992.ap.GOCC.gsea, showCategory=15)
dev.off()
pdf("DE-log2FC_8992genes_GOCC-GSEA_emapplot.pdf", width=15, height=15)
emapplot(de.8992.ap.GOCC.gsea, showCategory=15)
dev.off()
pdf("DE-log2FC_8992genes_GOCC-GSEA_gseaplot2.pdf", width=18, height=7)
gseaplot2(de.8992.ap.GOCC.gsea, geneSetID=c(1,5,7), title=se.up1.GOCC.gsea$Description[c(1,5,7)])
dev.off()
pdf("DE-log2FC_8992genes_GOCC-GSEA_heatplot.pdf", width=18, height=15)
heatplot(de.8992.ap.GOCC.gsea, foldChange=res.apeglm.8992.filt.FC)
dev.off()


de.8992.ap.SLPC.gsea <- GSEA(res.apeglm.8992.filt.FC, minGSSize=10, maxGSSize=1500, pvalueCutoff=0.1, pAdjustMethod="BH", TERM2GENE=data.frame(term=c(sl.genes$ID,polycis.genes$ID), gene=c(sl.genes$gene_id,polycis.genes$gene_id)))
as.data.frame(de.8992.ap.SLPC.gsea)
## NA





sessionInfo()
#R version 4.0.5 (2021-03-31)
#Platform: x86_64-redhat-linux-gnu (64-bit)
#Running under: Fedora 34 (Workstation Edition)

#Matrix products: default
#BLAS/LAPACK: /usr/lib64/libflexiblas.so.3.0

#locale:
# [1] LC_CTYPE=es_ES.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=es_ES.UTF-8        LC_COLLATE=es_ES.UTF-8    
# [5] LC_MONETARY=es_ES.UTF-8    LC_MESSAGES=es_ES.UTF-8   
# [7] LC_PAPER=es_ES.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#[11] LC_MEASUREMENT=es_ES.UTF-8 LC_IDENTIFICATION=C       

#attached base packages:
#[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
#[8] methods   base     

#other attached packages:
# [1] cowplot_1.1.1               enrichplot_1.10.2          
# [3] clusterProfiler_3.18.1      data.table_1.14.0          
# [5] biomaRt_2.46.3              seqinr_4.2-8               
# [7] FactoMineR_2.4              factoextra_1.0.7           
# [9] vsn_3.58.0                  pheatmap_1.0.12            
#[11] DESeq2_1.30.1               SummarizedExperiment_1.20.0
#[13] Biobase_2.50.0              MatrixGenerics_1.2.1       
#[15] matrixStats_0.59.0          GenomicRanges_1.42.0       
#[17] GenomeInfoDb_1.26.7         IRanges_2.24.1             
#[19] S4Vectors_0.28.1            BiocGenerics_0.36.1        
#[21] ggplot2_3.3.3               RColorBrewer_1.1-2         
#[23] gridExtra_2.3              

#loaded via a namespace (and not attached):
#  [1] readxl_1.3.1           shadowtext_0.0.8       backports_1.2.1       
#  [4] fastmatch_1.1-0        BiocFileCache_1.14.0   plyr_1.8.6            
#  [7] igraph_1.2.6           splines_4.0.5          BiocParallel_1.24.1   
# [10] digest_0.6.27          htmltools_0.5.1.1      GOSemSim_2.16.1       
# [13] viridis_0.6.1          GO.db_3.12.1           fansi_0.5.0           
# [16] magrittr_2.0.1         memoise_2.0.0          cluster_2.1.2         
# [19] openxlsx_4.2.4         limma_3.46.0           annotate_1.68.0       
# [22] graphlayouts_0.7.1     bdsmatrix_1.3-4        askpass_1.1           
# [25] prettyunits_1.1.1      colorspace_2.0-1       apeglm_1.12.0         
# [28] blob_1.2.1             rappdirs_0.3.3         ggrepel_0.9.1         
# [31] haven_2.4.1            dplyr_1.0.6            crayon_1.4.1          
# [34] RCurl_1.98-1.3         scatterpie_0.1.6       genefilter_1.72.1     
# [37] survival_3.2-11        glue_1.4.2             polyclip_1.10-0       
# [40] gtable_0.3.0           zlibbioc_1.36.0        XVector_0.30.0        
# [43] DelayedArray_0.16.3    car_3.0-11             abind_1.4-5           
# [46] scales_1.1.1           DOSE_3.16.0            mvtnorm_1.1-2         
# [49] DBI_1.1.1              rstatix_0.7.0          Rcpp_1.0.6            
# [52] emdbook_1.3.12         viridisLite_0.4.0      xtable_1.8-4          
# [55] progress_1.2.2         flashClust_1.01-2      foreign_0.8-81        
# [58] bit_4.0.4              preprocessCore_1.52.1  DT_0.18               
# [61] htmlwidgets_1.5.3      httr_1.4.2             fgsea_1.16.0          
# [64] ellipsis_0.3.2         pkgconfig_2.0.3        XML_3.99-0.6          
# [67] farver_2.1.0           dbplyr_2.1.1           locfit_1.5-9.4        
# [70] utf8_1.2.1             labeling_0.4.2         tidyselect_1.1.1      
# [73] rlang_0.4.11           reshape2_1.4.4         AnnotationDbi_1.52.0  
# [76] cellranger_1.1.0       munsell_0.5.0          tools_4.0.5           
# [79] cachem_1.0.5           downloader_0.4         generics_0.1.0        
# [82] RSQLite_2.2.7          ade4_1.7-17            broom_0.7.7           
# [85] stringr_1.4.0          fastmap_1.1.0          bit64_4.0.5           
# [88] tidygraph_1.2.0        zip_2.2.0              purrr_0.3.4           
# [91] ggraph_2.0.5           DO.db_2.9              leaps_3.1             
# [94] xml2_1.3.2             compiler_4.0.5         curl_4.3.1            
# [97] affyio_1.60.0          ggsignif_0.6.2         tibble_3.1.2          
#[100] tweenr_1.0.2           geneplotter_1.68.0     stringi_1.6.2         
#[103] forcats_0.5.1          lattice_0.20-44        Matrix_1.3-4          
#[106] vctrs_0.3.8            pillar_1.6.1           lifecycle_1.0.0       
#[109] BiocManager_1.30.15    bitops_1.0-7           qvalue_2.22.0         
#[112] R6_2.5.0               affy_1.68.0            rio_0.5.27            
#[115] MASS_7.3-54            assertthat_0.2.1       openssl_1.4.4         
#[118] withr_2.4.2            GenomeInfoDbData_1.2.4 hms_1.1.0             
#[121] grid_4.0.5             coda_0.19-4            tidyr_1.1.3           
#[124] rvcheck_0.1.8          carData_3.0-4          ggpubr_0.4.0          
#[127] bbmle_1.0.24           ggforce_0.3.3          numDeriv_2016.8-1.1   
#[130] scatterplot3d_0.3-41  


