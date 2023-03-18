##WGCNA optimization script optimizes TOM_calculation_time
#luoshaofan edited in 20230301
options(stringsAsFactors = FALSE)
library(WGCNA)
library(reshape2)
library(stringr)
allowWGCNAThreads()
dataExpr<-read.csv("Expression_data.csv",sep=",",header=T,row.names=1)
dataTrait<-read.delim("Trait_data.txt",sep="\t",header=T,row.names=1)
dataExpr = as.data.frame(t(dataExpr))
dataTrait = as.data.frame(t(dataTrait))
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)
nGenes = ncol(dataExpr)    
nSamples = nrow(dataExpr)
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK)
{
      if (sum(!gsg$goodGenes)>0)
       printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ", ")));
        if (sum(!gsg$goodSamples)>0)
      printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ", ")));
         dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
gsg = goodSamplesGenes(dataExpr, verbose = 3) 
gsg$allOK
sampleTree = hclust(dist(dataExpr), method = "average")
pdf(file="sampleTree.pdf",width=9,height=5)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, verbose=5)
pdf(file="power.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"));text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
power = sft$powerEstimate
power


adjacency_matrix = adjacency(datExpr, power = power);
source("Replace_the_TOMsimilarity_function_of_WGCNA.R")
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
pdf(file="geneTree.pdf",width=9,height=5)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04);
dev.off()
minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize);
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="DendroAndColors.pdf",width=9,height=5)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()
#model with trait

MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf(file="Clustering_of_module_eigengenes.pdf",width=9,height=5)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")
dev.off()
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
pdf(file="CplotDendroAndColors.pdf",width=9,height=5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
save(MEs, moduleLabels, moduleColors, geneTree, file = "AS-green-FPKM-02-networkConstruction-stepByStep.RData")
unique(moduleColors) -> modules
rownames(dataExpr)
sampleName = rownames(dataExpr)
dataTrait = dataTrait[match(sampleName, rownames(dataTrait)), ]
modTraitCor = cor(MEs_col, dataTrait, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(dataTrait),yLabels = colnames(MEs_col),cex.lab = 0.5,ySymbols = colnames(MEs_col), colorLabels = FALSE,colors = blueWhiteRed(50),textMatrix = textMatrix, setStdMargins = FALSE,cex.text = 0.5, zlim = c(-1,1),main = paste("Module-trait relationships"))
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = ("AS-green-FPKM-Step-by-step-CytoscapeInput-edges-all_model.txt"),
                               nodeFile = ("AS-green-FPKM-Step-by-step-CytoscapeInput-nodes-all_model.txt"),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);