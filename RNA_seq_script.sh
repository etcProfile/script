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
