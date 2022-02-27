#####
#Loading Libraries
library(affy)
library(simpleaffy)
library(affyPLM)
library(affyQCReport)
library(sva)
library(ggplot2)
library(pheatmap)
library(hgu133plus2.db)
library(AnnotationDbi)
library(dplyr)
library(WGCNA)
library(genefilter)
library(limma)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)

#####
#Loading in Data

raw <- ReadAffy(celfile.path = "FinalSamples")
# meta <- read.csv(file = "metadata.csv")

#####
#Normalizing the data

normdata <- rma(raw)

#####
#Quality Control

QCReport(raw, file = "QCReportAfterRemoval")

plmdata <- fitPLM(raw)

rledata <- RLE(plmdata, type = "stats")
hist(rledata)
boxplot(rledata)


nusedata <- NUSE(plmdata, type = "stats")
hist(nusedata)
boxplot(nusedata)

#####
#Box Plots

boxplot(raw)
boxplot(normdata)

# #####
# #Batch Correction
# normmatrix <- exprs(normdata)
# 
# modeldata <- model.matrix(~CN-1, data = meta)
# 
# combatdata <- ComBat(dat = normmatrix, batch = meta$BATCH, mod = modeldata[,1])
# 
# write.csv(combatdata, file = "CombatData.csv")

#####
#PCA Plot
prcompdata <- prcomp(normdata)
prcompdf <- as.data.frame(prcompdata$rotation)
group <- factor(c(rep("control", 16), rep("cancer", 17))) 
name <- factor(c(rownames(prcompdf)))
label <- c(1:33)
ggplot(prcompdf) + geom_point(aes(x = PC1, y = PC2, color = group)) + ggtitle("PC1 versus PC2 PCA Plot") + geom_text(aes(x = PC1, y = PC2, label = label), hjust = 1.5)

#####
#Heat Maps
intensitydata <- exprs(normdata)

rld_cor <- 1 - cor(intensitydata)
heatmap <- pheatmap(rld_cor)

#####
#Gene Annotation

normdataga <- exprs(normdata)
normdataga <- data.frame(PROBEID = c(rownames(normdataga)), normdataga)

nmdatarn <- rownames(normdataga)
selectdata <- AnnotationDbi::select(hgu133plus2.db, keys = nmdatarn, columns = c("SYMBOL"), keytypes = "PROBEID")

finalgenematrix <- selectdata
finalgenematrix <- finalgenematrix[!duplicated(selectdata$PROBEID), ]
# finalgenematrix <- finalgenematrix[!duplicated(selectdata$SYMBOL), ]

normdatags <- full_join(finalgenematrix, normdataga, by = "PROBEID", copy = F)

rownames(normdatags) <- normdatags[["PROBEID"]]

normdatags <- normdatags[complete.cases(normdatags),]

gene_vector_new <- normdatags[["SYMBOL"]]
gene_id_vector <- normdatags[["PROBEID"]]

drops <- c("PROBEID","SYMBOL")
normdatags <- normdatags[ , !(names(normdatags) %in% drops)]

new_matrix <- collapseRows(normdatags, gene_vector_new, gene_id_vector)

new_matrix2 <- as.matrix(new_matrix$datETcollapsed)

#####
#Gene Filtering

gene_matrix <- na.omit(new_matrix2, na.action = "omit", fill = NULL)
new_eset <- ExpressionSet(assayData=as.matrix(gene_matrix))
first_gene_vector <- as.vector(gene_matrix[1,])
quantile(first_gene_vector)
filt5 <- genefilter::varFilter(new_eset, var.func = IQR, var.cutoff = 0.02, filterByQuantile = TRUE)

#####
#Limma Analysis

new_gene_matrix <- exprs(filt5)

degem <- filt5
type <-  c(rep("control", 16), rep("cancer", 17))

fl <- as.factor(type)

degem$description <- fl
design <- model.matrix(~description + 0, degem)

colnames(design) <- levels(fl)

fit <- lmFit(degem, design)
contrast_matrix <- makeContrasts(cancer-control, levels = design)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

final_matrix <- topTable(fit2,coef=1, n=20542, sort = "logFC")

volcano_plot <- EnhancedVolcano(final_matrix,
                                lab = rownames(final_matrix),
                                pCutoff = 0.01,
                                FCcutoff = 0.5,
                                x = "logFC",
                                y = "adj.P.Val",
                                pointSize = 1,
                                legendLabSize = 10,
                                labSize = 3.0)
volcano_plot

genevector <- as.data.frame(final_matrix[,1])
rownames(genevector) <- rownames(final_matrix)
colnames(genevector) <- c("logFC")
genevector <- arrange(genevector, -logFC)

upreg <- rownames(genevector)[genevector$logFC < -0.5]

topgenes <- topTable(fit2, coef=1, p.value = 0.05, number = 50)

heatmap_matrix <- new_gene_matrix[(rownames(new_gene_matrix) %in% rownames(topgenes)), ]

#heatmap_cor <- 1 - cor(heatmap_matrix)
pheatmap(heatmap_matrix)

#####
#Gene Ontology Analysis

degene <- rownames(genevector)[genevector$logFC > 1.5]
View(degene)
write(degene, file = "genenames.txt", sep = "\n")

# entrezid <- bitr(geneID = degene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop = TRUE)
entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys = degene, columns = c("ENTREZID"), keytype = "SYMBOL")

entrezidvector <- entrezid[["ENTREZID"]]

entrezidvector2 <- filter(genevector, logFC > 1.5)
rownames(entrezidvector2) <- entrezid$ENTREZID

entrezidvector2 <- arrange(entrezidvector2, -logFC)

entrezidvector3 <- entrezidvector2$logFC
names(entrezidvector3) <- as.character(entrezid$ENTREZID)
entrezidvector3 <- sort(entrezidvector3, decreasing = TRUE)


gogroup <- groupGO(entrezidvector, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", level = 3, ont = "MF")

gogroupread <- setReadable(gogroup, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

barplot(gogroupread)

enrichgo <- enrichGO(entrezidvector, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = TRUE)

summary(enrichgo)

write.csv(as.data.frame(summary(enrichgo)), file = "enrichGOResults.csv")

enrichgo_noreadable <- enrichGO(entrezidvector, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID", readable = FALSE, ont = "MF")
enrichgo_readable <- setReadable(enrichgo_noreadable, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

dotplot(enrichgo_readable)
barplot(enrichgo_readable)

plotGOgraph(enrichgo_readable)

keggenrich <- enrichKEGG(entrezidvector)
write.csv(as.data.frame(summary(keggenrich)), file = "KeggEnrichResults.csv")
dotplot(keggenrich)

keggenrichreadable <- setReadable(keggenrich, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(keggenrichreadable, circular = TRUE)
?cnetplot

# gmtdata <- read.gmt("hallmark.gmt")
# gmtmatrix <- as.data.frame(gmtdata)
# 
# gseadata <- GSEA(entrezidvector3, TERM2GENE = gmtdata, pvalueCutoff = 0.5)

msigdata <- msigdbr(species = "Homo sapiens", category = "H")

newmsigdata <- msigdata %>% select(gs_name, entrez_gene)

ent <- topTable(fit2, coef=1, p.value = 0.05, number = 25000)

entrows <- rownames(ent)

entdata <- AnnotationDbi::select(org.Hs.eg.db, keys = entrows, columns = c("ENTREZID"), keytype = "SYMBOL")
entdata <- entdata[complete.cases(entdata),]
entdata <- entdata[!duplicated(entdata$SYMBOL), ]


entvector <- ent[["logFC"]]
names(entvector) <- as.character(entdata[["ENTREZID"]])
entvector <- sort(entvector, decreasing = TRUE)

msiggseadata <- GSEA(entvector, TERM2GENE = newmsigdata)

msigdataread <- setReadable(msiggseadata, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")

gseaplot2(msigdataread, geneSetID = 1:5, pvalue_table = T, ES_geom = "dot", subplots = 1:3)

heatplot(msigdataread, showCategory = 30, foldChange = entvector)




# symbolsforcnet <- rownames(final_matrix)[final_matrix$logFC > 1]
# 
# cnetdata <- AnnotationDbi::select(org.Hs.eg.db, keys = symbolsforcnet, columns = c("ENTREZID"), keytype = "SYMBOL")
# 
# cnetvector <- cnetdata[["ENTREZID"]]
# 
# cnetvector2 <- filter(genevector, logFC > 1)
# rownames(cnetvector2) <- cnetdata$ENTREZID
# 
# cnetvector2 <- arrange(cnetvector2, -logFC)
# 
# cnetvector3 <- cnetvector2$logFC
# names(cnetvector3) <- as.character(cnetdata$ENTREZID)
# cnetvector3 <- sort(cnetvector3, decreasing = TRUE)

gmtdata <- read.gmt("c32.gmt")

enrich_c3 <- enricher(names(entrezidvector3),TERM2GENE = gmtdata)

enricherreadable <- setReadable(enrich_c3, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnetplot(enricherreadable, showCategory = 5, foldChange = entrezidvector3, colorEdge = T)

