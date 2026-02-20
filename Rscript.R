file.exists("C:/Users/ROSET/OneDrive/Desktop/RNAseqdata/GSE164073.tsv")

gene_counts <- read.delim("C:/Users/roset/OneDrive/Desktop/RNAseqdata/GSE164073.TSV", header=TRUE, stringsAsFactors=FALSE)
head(gene_counts)
dim(gene_counts)
colnames(gene_counts)

rownames(gene_counts) <- gene_counts$GeneID
top_genes <- gene_counts[order(rowSums(gene_counts[,-1]), decreasing=TRUE), ][1:10, ]
gene_names <- as.character(top_genes$GeneID)
counts <- as.numeric(rowSums(top_genes[,-1]))
df_plot <- data.frame(Gene = gene_names, Counts = counts)
library(ggplot2)
ggplot(df_plot, aes(x = Gene, y = Counts)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=10)) +
  labs(title="Top 10 Genes by Total Counts", y="Sum of Counts", x="Gene")


library(DESeq2)
sample_info <- data.frame(
  sample = colnames(gene_counts)[-1]
)
dds <- DESeqDataSetFromMatrix(countData = gene_counts[,-1],
                              colData = sample_info,
                              design = ~ 1)

res <- DESeq(dds)
res <- results(res)
res <- res[order(res$pvalue), ]
sig_genes <- res[which(res$padj < 0.05), ]
res_clean <- res[!is.na(res$padj), ]
res_clean$threshold <- as.factor(res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1)
ggplot(res_clean, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
  geom_point(alpha=0.5, size=1.5) +
  scale_color_manual(values=c("grey","red")) +
  theme_minimal() +
  xlab("log2 Fold Change")







library(tidyr)
library(ggplot2)
sample_cols <- colnames(gene_counts)[-1]
top_genes <- gene_counts[order(rowSums(gene_counts[, sample_cols]), decreasing=TRUE), ][1:20, ]
df_melted <- pivot_longer(top_genes, cols = all_of(sample_cols), names_to = "Sample", values_to = "Counts")
df_melted$Sample <- factor(df_melted$Sample, levels = sample_cols)
df_melted$GeneID <- factor(df_melted$GeneID, levels = rev(top_genes$GeneID))

ggplot(df_melted, aes(x = Sample, y = GeneID, fill = Counts)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Top 20 Genes Heatmap", x = "Sample", y = "Gene", fill = "Counts")





library(ggplot2)
res_clean <- res[!is.na(res$padj), ]
res_clean$threshold <- as.factor(res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1)

ggplot(res_clean, aes(x = baseMean, y = log2FoldChange, color = threshold)) +
  geom_point(alpha = 0.5, size = 1.5) +
  scale_x_log10() +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  labs(title = "MA Plot", x = "Mean Expression", y = "log2 Fold Change")

  

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
mat <- mat[rowSums(mat) > 0, ]
pca_res <- prcomp(t(mat))

pca_data <- data.frame(PC1 = pca_res$x[,1],
                       PC2 = pca_res$x[,2],
                       Sample = colnames(mat))

ggplot(pca_data, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(title = "PCA Plot")



library(pheatmap)
top_genes <- head(res[order(res$padj), ], 20)
mat <- assay(vsd)[rownames(top_genes), ]
pheatmap(mat, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE,
         color = colorRampPalette(c("white", "red"))(100))






top_genes <- head(res[order(res$padj), ], 5)
mat <- assay(vsd)[rownames(top_genes), ]
df <- as.data.frame(t(mat))
df$Sample <- rownames(df)
df_long <- tidyr::pivot_longer(df, cols = -Sample, names_to = "Gene", values_to = "Expression")

ggplot(df_long, aes(x = Gene, y = Expression, fill = Gene)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Top 5 Gene Expression")





sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)
rownames(sample_dist_matrix) <- colnames(vsd)
colnames(sample_dist_matrix) <- colnames(vsd)

pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         color = colorRampPalette(c("white", "red"))(100))



sample_cor <- cor(assay(vsd))
pheatmap(sample_cor,
         color = colorRampPalette(c("white", "blue"))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE)






library(clusterProfiler)
library(org.Hs.eg.db)

gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_res <- gseGO(geneList = gene_list,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  nPerm = 1000,
                  minGSSize = 10,
                  maxGSSize = 500,
                  pvalueCutoff = 0.05,
                  verbose = FALSE)

dotplot(gsea_res, showCategory = 20)





library(WGCNA)
library(flashClust)
options(stringsAsFactors = FALSE)

datExpr <- t(assay(vsd))
gsg <- goodSamplesGenes(datExpr, verbose=3)
datExpr <- datExpr[, apply(datExpr, 2, function(x) all(!is.na(x)) & sd(x) > 0)]
gsg <- goodSamplesGenes(datExpr, verbose=3)

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2",
     type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.9, col="red")

topVarGenes <- head(order(rowVars(datExpr), decreasing=TRUE), 5000)
datExpr_small <- datExpr[, topVarGenes]

softPower <- 6
bwnet <- blockwiseModules(datExpr_small, power=softPower,
                          TOMType="unsigned", minModuleSize=30,
                          reassignThreshold=0, mergeCutHeight=0.25,
                          numericLabels=TRUE, saveTOMs=TRUE,
                          saveTOMFileBase="TOM", verbose=3)

moduleColors <- labels2colors(bwnet$colors)
plotDendroAndColors(bwnet$dendrograms[[1]],
                    moduleColors[bwnet$blockGenes[[1]]],
                    "Dynamic Tree Cut", dendroLabels=FALSE,
                    hang=0.03, addGuide=TRUE, guideHang=0.05)






MEs <- moduleEigengenes(datExpr_small, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs)

moduleTraitCor <- cor(MEs, traitData, use = "pairwise.complete.obs")



moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples = nrow(datExpr_small))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)



png("ModuleTraitHeatmap.png", width=2000, height=1600, res=150)

par(mar=c(5,4,4,2)+0.1)

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traitData),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = "Module-trait relationships")

dev.off()

getwd()












