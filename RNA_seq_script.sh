#The flow of script is  Next-generation sequencing for RNA-seq, including comparison, quantification, output fpkm, TPM, and differential gene analysis.
#author:LUOsf 
#date:2023.2.15 change hisat2 scrpit name
#2022.11.07 edit add DEseq2 and PCA 
############################################################################################################
DIR=~/$PREFIX/your_NGS_RNA_seq_cleanData
GENOME = ~/$PREFIX/YOUR_GENOME_ASSEMBLY&ANNOATION_FILE
#Use the scripts extract_splice_sites.py and extract_ xon.py that come from hisat2 to extract alternative splicing and exon information in the reference genome annotation file *.GFF or .GTF file*.
hisat2_extract_splice_sites.py $GENOME/genome.gtf > genome.ss
hisat2_extract_exons.py $GENOME/genome.gtf > genome.exon
#use hisat2 building genome index
hisat2-build -p 12 --ss genome.ss --exon genome.exon $GENOME/assembly.fna  Genome
#hisat2 mapping RNA-seq reads to reference genome's annotation.
hisat2 -x $GENOME/Genome -p 40 --known-splicesite-infile $GENOME/genome.ss --no-unal --new-summary --summary-file  -t -1 $DIR/R1_clean.fq.gz -2  $DIR/R2_clean.fq.gz | samtools view -bS - > r.bam;
#converent bam file to sortec bam file.
samtools sort -@ 20  r.bam  -o r.sorted.bam
#use featureCounts to counts read in each exon .
featureCounts -T 50 -p -t exon -g gene_id -a $GENOME/genome.gtf  -o counts.txt r.sorted.bam 
#Here should be your all sorted bam files.
#Finally get counts result *counts.txt*.
############################################################################################################
#prepare your different expression genes analysis data.
# all by R, use Bioconductor to install follow packages.
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#BiocManager:install(c("eoffice","DESeq2","ggplot2","VennDiagram","tidyverse","dplyr")).
library(eoffice) #For output office format result.
library(DESeq2) #For differenet expression genes analysis.
library(ggplot2) 
library(VennDiagram) # For Venn Plot.
library(tidyverse) # tidyverse
library(dplyr)

data<-read.table("./counts.txt",sep="\t",skip=1,header=T)
countData <- as.matrix(data[7:12])
rownames(countData) <- data$Geneid
group<-read.table("./group.txt",sep="\t",header=T) #group.txt is your samples info file.
#if your dont have group file, you can build it by R.
#For Example:
#sampleNames <- c("1_1","1_2","1_3","2_a","2_b","2_c","3_1","3_2","3_3")
#group <- data.frame(name=sampleNames,condition=c("1","1","1","2","2","2","3","3","3"))
#rownames(group) <- sampleNames
g<-as.vector(group$Group);
names(g)<-group$Name;
condition<-factor(g[colnames(dat)])
condition<-factor(group[,2])
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
dds<-DESeq(dds)
#product versus list, If you have three samples, the value of the list is three according to the pairwise comparison=n?-n.
ll<-combn(unique(as.vector(group[,2])),2,simplify=F)
res<-lapply(ll,function(x)assign(paste(x,sep="",collapse="_vs_"),results(dds, contrast=c("condition",x))))
names(res)<-unlist(lapply(ll,function(x)paste(x,sep="",collapse="_vs_"))) 
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
dir.create("DE")
write.table(resdata,file="DE/All_DE_group.txt",sep="\t",quote=F);
for(i in names(res)){
    write.table(res[[i]],file=paste(i,".txt",sep=""),quote=F)
}
rld <- rlogTransformation(dds)

library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])
sampleDists <- as.matrix(dist(t(assay(rld))))

#valcano_plot
library(ggplot2)
resdata$change <- as.factor(
	ifelse(
		resdata$padj<0.01 & abs(resdata$log2FoldChange)>1,
		ifelse(resdata$log2FoldChange>1, "Up", "Down"),
		"NoDiff"
	)
)
valcano <- ggplot(data=resdata, aes(x=log2FoldChange, y=-log10(padj), color=change)) + 
	geom_point(alpha=0.8, size=1) + 
	theme_bw(base_size=15) + 
	theme(
		panel.grid.minor=element_blank(),
		panel.grid.major=element_blank()
	) + 
	ggtitle("DESeq2 Valcano") + 
	scale_color_manual(name="", values=c("red", "green", "black"), limits=c("Up", "Down", "NoDiff")) + 
	geom_vline(xintercept=c(-1, 1), lty=2, col="gray", lwd=0.5) + 
	geom_hline(yintercept=-log10(0.01), lty=2, col="gray", lwd=0.5)

topptx(valcano,filenames="valcano.pptx",width=6,hight=6)

########sample_pca_plot#######
#https://github.com/sachitsaksena/nfkb-tag-seq/blob/master/rld_pca.R
library(gplots)
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  print(length(fac))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  # with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  # legend(legendpos, legend=levels(fac), col=colors, pch=20, ncol = 4)
  # rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #       pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                   terldt = list(levels(fac)), rep = FALSE)))
  return(pca)
}
pdf(file="group_pca.pdf")
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-300, 300),ylim=c(-300,300),textcx=0.5)
dev.off()

############################################################################################################
#Convert featureCounts result "Count.txt" to "Gene_expression_FPKM_or_TPM".
## First you could use wheel package to convert counts result to FPKM, then Manual convert to the TPM format.
## Install package.
## Stable version, install countToFPKM from CRAN.
install.packages("countToFPKM")
## Lastest version, install countToFPKM from GitHub
if(!require(devtools)) install.packages("devtools")
devtools::install_github("AAlhendi1707/countToFPKM", build_vignettes = TRUE)
library(countToFPKM)
##Example from Github AAlhendi1707/countToFPKM".
##file.readcounts <- system.file("extdata", "RNA-seq.read.counts.csv", package="countToFPKM")
##file.annotations <- system.file("extdata", "Biomart.annotations.hg38.txt", package="countToFPKM")
##file.sample.metrics <- system.file("extdata", "RNA-seq.samples.metrics.txt", package="countToFPKM")
# Import the read count matrix data into R.
counts <- as.matrix(read.csv(counts.txt))
# Import feature annotations. 
# Assign feature lenght into a numeric vector.
gene.annotations <- read.table(gff_file, sep="\t", header=TRUE)
featureLength <- gene.annotations$length
# Import sample metrics. 
# Assign mean fragment length into a numeric vector.
samples.metrics <- read.table(file.sample.metrics, sep="\t", header=TRUE)
meanFragmentLength <- samples.metrics$meanFragmentLength
# Return FPKM into a numeric matrix.
fpkm_matrix <- fpkm (counts, featureLength, meanFragmentLength)
############################################################################################################
countdata<-read.table("counts.txt",skip = 1,sep="\t",header = T,row.names = 1)
metadata <- countdata[,1:5]
countdata <- countdata[,6:ncol(countdata)]
kb <- metadata$Length / 1000
rpk <- countdata / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
write.csv(avg_tpm,paste0(prefix,"_avg_tpm.csv"))
write.csv(tpm,paste0(prefix,"_tpm.csv"))
###################################################################################################
## functions for rpkm and tpm
## from https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44
## from https://www.biostars.org/p/171766/
rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e9
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}

## read table from featureCounts output
args <- commandArgs(T)
tag <- tools::file_path_sans_ext(args[1])
ftr.cnt <- read.table(args[1], sep="\t", stringsAsFactors=FALSE,header=TRUE)
library(dplyr)
library(tidyr)

ftr.rpkm <- ftr.cnt %>%
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%
  group_by(sample) %>%
  mutate(rpkm=rpkm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, rpkm)
write.table(ftr.rpkm, file=paste0(tag, "_rpkm.txt"), sep="\t", row.names=FALSE, quote=FALSE)

ftr.tpm <- ftr.cnt %>%
  gather(sample, cnt, 7:ncol(ftr.cnt)) %>%
  group_by(sample) %>%
  mutate(tpm=tpm(cnt, Length)) %>%
  select(-cnt) %>%
  spread(sample, tpm)
write.table(ftr.tpm, file=paste0(tag, "_tpm.txt"), sep="\t", row.names=FALSE, quote=FALSE)