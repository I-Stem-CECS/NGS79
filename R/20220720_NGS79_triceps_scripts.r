# 2022-07-20 NGS79 
# tissus specific : triceps
# RNAseq 


library(RColorBrewer)
library(gplots)
library(DESeq2)
library("pheatmap")
library("gridExtra")
library(reshape)
library("tidyverse")
library("ggforce")
library("cowplot")
library("reshape2")
library(viridis)

rawcts <- read.csv2("../results/generalView/NGS79_readscounts_raw_55471enes.csv", header = T)
rawcts <- rawcts[, c(2:ncol(rawcts))]

names_genes <- read.csv2("../../genomes/GRCm38.99/geneToSymbol_GRCm38.99.csv", sep="\t", header=F , stringsAsFactors=F)
colnames(names_genes) <- c("id_ensembl","name_gene")

coldata <- read.csv2("./coldata.csv", sep="\t", header=T)
rownames(coldata) <- coldata$sampleName

coldata <- coldata[which(coldata$tissus %in% "Triceps"),]

coldata$condition <- as.factor(coldata$condition)
coldata$Genotype <- as.factor(coldata$Genotype)
coldata$tissus <- as.factor(coldata$tissus)
coldata$Breading <- as.factor(coldata$Breading)
coldata$mouse <- as.factor(coldata$mouse)
coldata$Treatment <- as.factor(coldata$Treatment)


rawcts <- rawcts[,c(rownames(coldata))]

cts.1 <- data.frame(rawcts, moy=NA)
for (i in 1:nrow(cts.1)){
  cts.1[i, "moy"] <- mean(as.numeric(cts.1[i , c(1:ncol(rawcts))]))
}
cts.2 <- cts.1[which(cts.1$moy > 5),]
cts.filtre <- cts.2[, c(1:ncol(rawcts))] # 14019
colnames(cts.filtre) <- colnames(rawcts)

##################################
####### dds cts.filtre coldata ###
##################################

cts.filtre <- cts.filtre[,c(order(colnames(cts.filtre)))]
coldata <- coldata[order(coldata$sampleName),]

dds <- DESeqDataSetFromMatrix(countData = cts.filtre, colData = coldata,
                              design = ~ condition )
dds <- DESeq(dds)

###################################
## Comptages bruts et Normalisés ##
###################################

#### Bruts
counts_raw_HP <- counts(dds, norm=F)

counts_raw <- data.frame(symbol=NA, counts_raw_HP)
for (i in 1:nrow(counts_raw) ){
        counts_raw[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(counts_raw)[i]),2]
}

write.table(counts_raw,file="../results/triceps/NGS79_Triceps_readscounts_raw_14019genes.csv", sep=";")

#### Norm
counts_normalise <- counts(dds, norm=T)
counts_normalise_HP <- data.frame(symbol=NA, counts_normalise)
for (i in 1:nrow(counts_normalise) ){
        counts_normalise_HP[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(counts_normalise_HP)[i]),2]
}

write.table(counts_normalise_HP,file="../results/triceps/GS79_Tricpes_readscounts_norm_14019genes.csv", sep=";", dec=",")


#### GENERAL VIEW ###

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)

pdf("../results/triceps/plots_general_view_Tricpes_NGS79.pdf")

condcols=brewer.pal(n = length(unique(coldata$condition)), name = 'Paired')
names(condcols)=unique(coldata$condition)

barplot(colSums(counts(dds, normalized=F)), col=condcols[as.factor(coldata$condition)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[as.factor(coldata$condition)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')

pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )

colnames(sampleDistMatrix) <- c("E1_CTL_WT", "E11_CTL_KO", "E12_CTL_KO", "E13_CTL_KO", 
                           "E14_CTL_KO", "E15_RAPA_KO", "E16_RAPA_KO", "E17_RAPA_KO",
                           "E18_CTL_AAV_KO", "E19_CTL_AAV_KO", "E20_CTL_AAV_KO", 
                           "E22_CTL_AAV_KO", "E23_RAPA_AAV_KO", "E24_RAPA_AAV_KO", "E25_RAPA_AAV_KO",
                           "E26_RAPA_AAV_KO", "E3_CTL_WT", "E5_CTL_WT")

rownames(sampleDistMatrix) <- c("E1_CTL_WT", "E11_CTL_KO", "E12_CTL_KO", "E13_CTL_KO", 
                           "E14_CTL_KO", "E15_RAPA_KO", "E16_RAPA_KO", "E17_RAPA_KO",
                           "E18_CTL_AAV_KO", "E19_CTL_AAV_KO", "E20_CTL_AAV_KO", 
                           "E22_CTL_AAV_KO", "E23_RAPA_AAV_KO", "E24_RAPA_AAV_KO", "E25_RAPA_AAV_KO",
                           "E26_RAPA_AAV_KO", "E3_CTL_WT", "E5_CTL_WT")

colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15), cexRow=0.5, cexCol=0.5)


plotDispEsts(dds)


dev.off()


###########################
### Counts observations ###
###########################

# counts heatmap

# Classic palette BuPu, with 4 colors
coul <- brewer.pal(4, "PuOr") 
# Add more colors to this palette :
coul <- colorRampPalette(coul)(8)

annotation_col <- data.frame( condition = coldata$condition, Breading = coldata$Breading) #, genotype = coldata$genotype )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(condition = brewer.pal(5, "Spectral"), Breading =  coul) #  , genotype = c("mediumpurple", "hotpink4"))
names(mat_colors$condition) <- unique(coldata$condition)
names(mat_colors$Breading) <- unique(coldata$Breading)

#cts.NGS66.2 <- cts.NGS66[,c(2:ncol(cts.NGS66))]
cts.norm.df <- as.data.frame(counts_normalise)
png("../results/triceps/NGS79_Triceps_pheatmap_counts_14019.png" , width = 1000, height = 1000)
pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, 
          cluster_cols=TRUE, annotation_col = annotation_col ,
          annotation_colors = mat_colors , angle_col = "45")
dev.off()


png("../results/triceps/NGS79_Triceps_pheatmap_counts_viridis_14019.png" , width = 1000, height = 1000)
pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, 
          cluster_cols=TRUE, annotation_col = annotation_col ,color = viridis(250),
          annotation_colors = mat_colors , angle_col = "45")
dev.off()

########################
### CTL_WT vs CTL_KO ###
########################

cdition1 <- "CTL_WT"
cdition2Ctrl <- "CTL_KO"
contrastO <- "condition"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
              which((colnames(cts.filtre) %in% coldata[ which(
                coldata$condition %in% cdition1), "sampleName"]) | 
		(colnames(cts.filtre) %in% coldata[ which(
		  coldata$condition %in% cdition2Ctrl ), "sampleName"])
		)], 
		colData = coldata[which((coldata$condition %in% cdition1) | 
		                          (coldata$condition %in% cdition2Ctrl) ),] , 
		design = ~ condition )



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), 
                 lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 205
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 181
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 24

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("../results/triceps/plots_Triceps_",cdition1,"-vs-",cdition2Ctrl,"_NGS79.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$condition)), name = 'Paired')
names(condcols)=unique(coldata$condition)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("condition","tissus"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15))


plotDispEsts(dds)

hits=rownames(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),])

plot(resGA$log2FoldChange,-log(resGA$padj,10),
     ylab='-log10(Adjusted P)',
     xlab="Log2 FoldChange",
     pch=19,cex=0.5, col = "dimgray"
     )      

points(resGA[hits,'log2FoldChange'],
       -log(resGA[hits,'padj'],10),
       pch=19,
       cex=0.5,
       col="firebrick"
       )
abline(h=-log10(0.05),lty=3)
abline(v=-0.4,lty=3)
abline(v=0.4,lty=3)


#plot the -log10(p-val) from all genes over the normalized mean counts 
plot(resGA$baseMean+1, -log10(resGA$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- resGA$baseMean > 10
table(use)
h1 <- hist(resGA$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(resGA$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


counts_normalise <- counts(dds, norm=T)
cts.norm.df <- as.data.frame(counts_normalise)

annotation_col <- data.frame( condition = coldata$condition )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(condition = c("firebrick", "steelblue"))
names(mat_colors$condition) <- c(cdition1, cdition2Ctrl)


pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

pheatmap( log10(cts.norm.df[hits,] + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

dev.off()

## Formatage gene symbol

res_sig <- as.data.frame(resGA)
mat_sig <- data.frame(symbol=NA, res_sig)
for (i in 1:nrow(mat_sig) ){
	mat_sig[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(mat_sig)[i]),2]
}

write.csv2(mat_sig, paste0("../results/triceps/",cdition1,"-vs-", cdition2Ctrl,"_Triceps_log2FC0-4.csv"))

## Formatage I-Stem

coldt <- coldata[which((coldata$condition %in% cdition1) | 
                         (coldata$condition %in% cdition2Ctrl) ),]

counts_norm_HP <- counts(dds, norm=T)
counts_norm_HP <- counts_norm_HP[,c(order(coldt$condition))]

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA)


for (i in 1:nrow(finalDE) ){
	if (is.na(finalDE[i,"log2FoldChange"])) {
		finalDE[i, "differential_expression"] <- NA
	} else if (finalDE[i, "log2FoldChange"] >= 0 ) {
		finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
	} else {
		finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
	}
	
}

finalDE <- merge(finalDE, counts_norm_HP, by = "row.names")
rownames(finalDE) <- finalDE$Row.names
finalDE <- finalDE[,c(2:ncol(finalDE))]

write.csv2(finalDE, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))



############################
### CTL_AAV_KO vs CTL_KO ###
############################

cdition1 <- "CTL_AAV_KO"
cdition2Ctrl <- "CTL_KO"
contrastO <- "condition"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$condition %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$condition %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$condition %in% cdition1) | 
                                                        (coldata$condition %in% cdition2Ctrl) ),] , 
                              design = ~ condition )



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), 
                 lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 0
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 0
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 0

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("../results/triceps/plots_Triceps_",cdition1,"-vs-",cdition2Ctrl,"_NGS79.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$condition)), name = 'Paired')
names(condcols)=unique(coldata$condition)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15))


plotDispEsts(dds)

hits=rownames(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),])

plot(resGA$log2FoldChange,-log(resGA$padj,10),
     ylab='-log10(Adjusted P)',
     xlab="Log2 FoldChange",
     pch=19,cex=0.5, col = "dimgray"
)      

points(resGA[hits,'log2FoldChange'],
       -log(resGA[hits,'padj'],10),
       pch=19,
       cex=0.5,
       col="firebrick"
)
abline(h=-log10(0.05),lty=3)
abline(v=-0.4,lty=3)
abline(v=0.4,lty=3)


#plot the -log10(p-val) from all genes over the normalized mean counts 
plot(resGA$baseMean+1, -log10(resGA$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- resGA$baseMean > 10
table(use)
h1 <- hist(resGA$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(resGA$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


counts_normalise <- counts(dds, norm=T)
cts.norm.df <- as.data.frame(counts_normalise)

annotation_col <- data.frame( condition = coldata$condition )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(condition = c("firebrick", "steelblue"))
names(mat_colors$condition) <- c(cdition1, cdition2Ctrl)


pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

pheatmap( log10(cts.norm.df[hits,] + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

dev.off()

## Formatage gene symbol

res_sig <- as.data.frame(resGA)
mat_sig <- data.frame(symbol=NA, res_sig)
for (i in 1:nrow(mat_sig) ){
  mat_sig[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(mat_sig)[i]),2]
}

write.csv2(mat_sig, paste0("../results/triceps/",cdition1,"-vs-", cdition2Ctrl,"_Triceps_log2FC0-4.csv"))

## Formatage I-Stem

coldt <- coldata[which((coldata$condition %in% cdition1) | 
                         (coldata$condition %in% cdition2Ctrl) ),]

counts_norm_HP <- counts(dds, norm=T)
counts_norm_HP <- counts_norm_HP[,c(order(coldt$condition))]

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA)


for (i in 1:nrow(finalDE) ){
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

finalDE <- merge(finalDE, counts_norm_HP, by = "row.names")
rownames(finalDE) <- finalDE$Row.names
finalDE <- finalDE[,c(2:ncol(finalDE))]

write.csv2(finalDE, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))





############################
### RAPA_AAV_KO vs CTL_KO ###
############################

cdition1 <- "RAPA_AAV_KO"
cdition2Ctrl <- "CTL_KO"
contrastO <- "condition"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$condition %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$condition %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$condition %in% cdition1) | 
                                                        (coldata$condition %in% cdition2Ctrl) ),] , 
                              design = ~ condition )



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), 
                 lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 33
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 30
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 3

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("../results/triceps/plots_Triceps_",cdition1,"-vs-",cdition2Ctrl,"_NGS79.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$condition)), name = 'Paired')
names(condcols)=unique(coldata$condition)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15))


plotDispEsts(dds)

hits=rownames(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),])

plot(resGA$log2FoldChange,-log(resGA$padj,10),
     ylab='-log10(Adjusted P)',
     xlab="Log2 FoldChange",
     pch=19,cex=0.5, col = "dimgray"
)      

points(resGA[hits,'log2FoldChange'],
       -log(resGA[hits,'padj'],10),
       pch=19,
       cex=0.5,
       col="firebrick"
)
abline(h=-log10(0.05),lty=3)
abline(v=-0.4,lty=3)
abline(v=0.4,lty=3)


#plot the -log10(p-val) from all genes over the normalized mean counts 
plot(resGA$baseMean+1, -log10(resGA$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- resGA$baseMean > 10
table(use)
h1 <- hist(resGA$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(resGA$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


counts_normalise <- counts(dds, norm=T)
cts.norm.df <- as.data.frame(counts_normalise)

annotation_col <- data.frame( condition = coldata$condition )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(condition = c("firebrick", "steelblue"))
names(mat_colors$condition) <- c(cdition1, cdition2Ctrl)


pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

pheatmap( log10(cts.norm.df[hits,] + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

dev.off()

## Formatage gene symbol

res_sig <- as.data.frame(resGA)
mat_sig <- data.frame(symbol=NA, res_sig)
for (i in 1:nrow(mat_sig) ){
  mat_sig[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(mat_sig)[i]),2]
}

write.csv2(mat_sig, paste0("../results/triceps/",cdition1,"-vs-", cdition2Ctrl,"_Triceps_log2FC0-4.csv"))

## Formatage I-Stem

coldt <- coldata[which((coldata$condition %in% cdition1) | 
                         (coldata$condition %in% cdition2Ctrl) ),]

counts_norm_HP <- counts(dds, norm=T)
counts_norm_HP <- counts_norm_HP[,c(order(coldt$condition))]

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA)


for (i in 1:nrow(finalDE) ){
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

finalDE <- merge(finalDE, counts_norm_HP, by = "row.names")
rownames(finalDE) <- finalDE$Row.names
finalDE <- finalDE[,c(2:ncol(finalDE))]

write.csv2(finalDE, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))



#################################
### RAPA_AAV_KO vs CTL_AAV_KO ###
#################################

cdition1 <- "RAPA_AAV_KO"
cdition2Ctrl <- "CTL_AAV_KO"
contrastO <- "condition"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$condition %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$condition %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$condition %in% cdition1) | 
                                                        (coldata$condition %in% cdition2Ctrl) ),] , 
                              design = ~ condition )



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), 
                 lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 1
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 0
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 1

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("../results/triceps/plots_Triceps_",cdition1,"-vs-",cdition2Ctrl,"_NGS79.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$condition)), name = 'Paired')
names(condcols)=unique(coldata$condition)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15))


plotDispEsts(dds)

hits=rownames(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),])

plot(resGA$log2FoldChange,-log(resGA$padj,10),
     ylab='-log10(Adjusted P)',
     xlab="Log2 FoldChange",
     pch=19,cex=0.5, col = "dimgray"
)      

points(resGA[hits,'log2FoldChange'],
       -log(resGA[hits,'padj'],10),
       pch=19,
       cex=0.5,
       col="firebrick"
)
abline(h=-log10(0.05),lty=3)
abline(v=-0.4,lty=3)
abline(v=0.4,lty=3)


#plot the -log10(p-val) from all genes over the normalized mean counts 
plot(resGA$baseMean+1, -log10(resGA$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- resGA$baseMean > 10
table(use)
h1 <- hist(resGA$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(resGA$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


counts_normalise <- counts(dds, norm=T)
cts.norm.df <- as.data.frame(counts_normalise)

annotation_col <- data.frame( condition = coldata$condition )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(condition = c("firebrick", "steelblue"))
names(mat_colors$condition) <- c(cdition1, cdition2Ctrl)


pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

pheatmap( log10(cts.norm.df[hits,] + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

dev.off()

## Formatage gene symbol

res_sig <- as.data.frame(resGA)
mat_sig <- data.frame(symbol=NA, res_sig)
for (i in 1:nrow(mat_sig) ){
  mat_sig[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(mat_sig)[i]),2]
}

write.csv2(mat_sig, paste0("../results/triceps/",cdition1,"-vs-", cdition2Ctrl,"_Triceps_log2FC0-4.csv"))

## Formatage I-Stem

coldt <- coldata[which((coldata$condition %in% cdition1) | 
                         (coldata$condition %in% cdition2Ctrl) ),]

counts_norm_HP <- counts(dds, norm=T)
counts_norm_HP <- counts_norm_HP[,c(order(coldt$condition))]

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA)


for (i in 1:nrow(finalDE) ){
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

finalDE <- merge(finalDE, counts_norm_HP, by = "row.names")
rownames(finalDE) <- finalDE$Row.names
finalDE <- finalDE[,c(2:ncol(finalDE))]

write.csv2(finalDE, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))


##############################
### RAPA_AAV_KO vs RAPA_KO ###
##############################

cdition1 <- "RAPA_AAV_KO"
cdition2Ctrl <- "RAPA_KO"
contrastO <- "condition"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$condition %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$condition %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$condition %in% cdition1) | 
                                                        (coldata$condition %in% cdition2Ctrl) ),] , 
                              design = ~ condition )



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), 
                 lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 25
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 23
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 2

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("../results/triceps/plots_Triceps_",cdition1,"-vs-",cdition2Ctrl,"_NGS79.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$condition)), name = 'Paired')
names(condcols)=unique(coldata$condition)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15))


plotDispEsts(dds)

hits=rownames(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),])

plot(resGA$log2FoldChange,-log(resGA$padj,10),
     ylab='-log10(Adjusted P)',
     xlab="Log2 FoldChange",
     pch=19,cex=0.5, col = "dimgray"
)      

points(resGA[hits,'log2FoldChange'],
       -log(resGA[hits,'padj'],10),
       pch=19,
       cex=0.5,
       col="firebrick"
)
abline(h=-log10(0.05),lty=3)
abline(v=-0.4,lty=3)
abline(v=0.4,lty=3)


#plot the -log10(p-val) from all genes over the normalized mean counts 
plot(resGA$baseMean+1, -log10(resGA$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- resGA$baseMean > 10
table(use)
h1 <- hist(resGA$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(resGA$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


counts_normalise <- counts(dds, norm=T)
cts.norm.df <- as.data.frame(counts_normalise)

annotation_col <- data.frame( condition = coldata$condition )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(condition = c("firebrick", "steelblue"))
names(mat_colors$condition) <- c(cdition1, cdition2Ctrl)


pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

pheatmap( log10(cts.norm.df[hits,] + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

dev.off()

## Formatage gene symbol

res_sig <- as.data.frame(resGA)
mat_sig <- data.frame(symbol=NA, res_sig)
for (i in 1:nrow(mat_sig) ){
  mat_sig[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(mat_sig)[i]),2]
}

write.csv2(mat_sig, paste0("../results/triceps/",cdition1,"-vs-", cdition2Ctrl,"_Triceps_log2FC0-4.csv"))

## Formatage I-Stem

coldt <- coldata[which((coldata$condition %in% cdition1) | 
                         (coldata$condition %in% cdition2Ctrl) ),]

counts_norm_HP <- counts(dds, norm=T)
counts_norm_HP <- counts_norm_HP[,c(order(coldt$condition))]

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA)


for (i in 1:nrow(finalDE) ){
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

finalDE <- merge(finalDE, counts_norm_HP, by = "row.names")
rownames(finalDE) <- finalDE$Row.names
finalDE <- finalDE[,c(2:ncol(finalDE))]

write.csv2(finalDE, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))



##############################
### CTL_WT vs RAPA_KO ###
##############################

cdition1 <- "CTL_WT"
cdition2Ctrl <- "RAPA_KO"
contrastO <- "condition"

message(cdition1)
message(cdition2Ctrl)
message(contrastO)

dds <- DESeqDataSetFromMatrix(countData = cts.filtre[ ,
                                                      which((colnames(cts.filtre) %in% coldata[ which(
                                                        coldata$condition %in% cdition1), "sampleName"]) | 
                                                          (colnames(cts.filtre) %in% coldata[ which(
                                                            coldata$condition %in% cdition2Ctrl ), "sampleName"])
                                                      )], 
                              colData = coldata[which((coldata$condition %in% cdition1) | 
                                                        (coldata$condition %in% cdition2Ctrl) ),] , 
                              design = ~ condition )



dds <- DESeq(dds)

resGA <- results(dds, contrast=c(contrastO,cdition1,cdition2Ctrl), 
                 lfcThreshold=0.4, altHypothesis="greaterAbs")

message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),]))
# 229
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange < 0 ),])) # Down
# 175
message(nrow(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 & resGA$log2FoldChange > 0 ),])) # Up
# 54

rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vstMat = assay(vsd)


pdf(paste0("../results/triceps/plots_Triceps_",cdition1,"-vs-",cdition2Ctrl,"_NGS79.pdf"))

hmcol = colorRampPalette(brewer.pal(9, 'GnBu'))(100)
condcols=brewer.pal(n = length(unique(coldata$condition)), name = 'Paired')
names(condcols)=unique(coldata$condition)

barplot(colSums(counts(dds, normalized=F)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Pre Normalised Counts')

barplot(colSums(counts(dds, normalized=T)), col=condcols[c(cdition1, cdition2Ctrl)], 
        las=2,cex.names=0.4,
        main='Post Normalised Counts')



ylim <- c(-10,10)
drawLines <- function() abline(h=c(-0.4,0.4),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()


pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()


sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours, margins=c(15,15))


plotDispEsts(dds)

hits=rownames(resGA[which((resGA$padj < 0.05) & resGA$baseMean > 20 ),])

plot(resGA$log2FoldChange,-log(resGA$padj,10),
     ylab='-log10(Adjusted P)',
     xlab="Log2 FoldChange",
     pch=19,cex=0.5, col = "dimgray"
)      

points(resGA[hits,'log2FoldChange'],
       -log(resGA[hits,'padj'],10),
       pch=19,
       cex=0.5,
       col="firebrick"
)
abline(h=-log10(0.05),lty=3)
abline(v=-0.4,lty=3)
abline(v=0.4,lty=3)


#plot the -log10(p-val) from all genes over the normalized mean counts 
plot(resGA$baseMean+1, -log10(resGA$pvalue),
     log="x", xlab="mean of normalized counts",
     ylab=expression(-log[10](pvalue)),
     ylim=c(0,30),
     cex=.4, col=rgb(0,0,0,.3))

use <- resGA$baseMean > 10
table(use)
h1 <- hist(resGA$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(resGA$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))


counts_normalise <- counts(dds, norm=T)
cts.norm.df <- as.data.frame(counts_normalise)

annotation_col <- data.frame( condition = coldata$condition )
rownames(annotation_col) <- rownames(coldata)
mat_colors <- list(condition = c("firebrick", "steelblue"))
names(mat_colors$condition) <- c(cdition1, cdition2Ctrl)


pheatmap( log10(cts.norm.df + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

pheatmap( log10(cts.norm.df[hits,] + 1), 
          cluster_rows=TRUE, show_rownames=FALSE, color = viridis(250),
          cluster_cols=TRUE, annotation_col = annotation_col , #fontsize = 12,
          annotation_colors = mat_colors , angle_col = "45"
)

dev.off()

## Formatage gene symbol

res_sig <- as.data.frame(resGA)
mat_sig <- data.frame(symbol=NA, res_sig)
for (i in 1:nrow(mat_sig) ){
  mat_sig[i,1] <- names_genes[which(names_genes$id_ensembl %in% rownames(mat_sig)[i]),2]
}

write.csv2(mat_sig, paste0("../results/triceps/",cdition1,"-vs-", cdition2Ctrl,"_Triceps_log2FC0-4.csv"))

## Formatage I-Stem

coldt <- coldata[which((coldata$condition %in% cdition1) | 
                         (coldata$condition %in% cdition2Ctrl) ),]

counts_norm_HP <- counts(dds, norm=T)
counts_norm_HP <- counts_norm_HP[,c(order(coldt$condition))]

finalDE <- data.frame(mat_sig[,c(1:3, 6,7)], differential_expression=NA)


for (i in 1:nrow(finalDE) ){
  if (is.na(finalDE[i,"log2FoldChange"])) {
    finalDE[i, "differential_expression"] <- NA
  } else if (finalDE[i, "log2FoldChange"] >= 0 ) {
    finalDE[i, "differential_expression"] <- 2^(finalDE[i, "log2FoldChange"])
  } else {
    finalDE[i, "differential_expression"] <- (-1)*(1/(2^(finalDE[i,"log2FoldChange"])))
  }
  
}

finalDE <- merge(finalDE, counts_norm_HP, by = "row.names")
rownames(finalDE) <- finalDE$Row.names
finalDE <- finalDE[,c(2:ncol(finalDE))]

write.csv2(finalDE, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4.csv"))
finalDE_padj05 <- finalDE[which(finalDE$padj < 0.05 ),]
finalDE_padj05_BM20 <- finalDE_padj05[which(finalDE_padj05$baseMean >= 20 ),]
write.csv2(finalDE_padj05_BM20, paste0("../results/triceps/NGS79_Triceps_", cdition1, "-vs-", cdition2Ctrl,"_log2FC0-4_padj0-5_BM20.csv"))











