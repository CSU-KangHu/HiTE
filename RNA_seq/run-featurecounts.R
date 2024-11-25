#!/usr/bin/env Rscript
# parse parameter ---------------------------------------------------------
# 检查并安装缺失的包
load_or_install <- function(package, from_bioconductor=FALSE) {
  if (!require(package, character.only = TRUE)) {
    if (from_bioconductor) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(package)
    } else {
      install.packages(package, dependencies = TRUE)
    }
    library(package, character.only = TRUE)
  }
}

# 加载所需的包
load_or_install("argparser")                   # 从 CRAN 安装
load_or_install("Rsubread", from_bioconductor=TRUE)  # 从 Bioconductor 安装
load_or_install("limma", from_bioconductor=TRUE)     # 从 Bioconductor 安装
load_or_install("edgeR", from_bioconductor=TRUE)     # 从 Bioconductor 安装


# Create a parser
p <- arg_parser("run featureCounts and calculate FPKM/TPM")

# Add command line arguments
p <- add_argument(p, "--bam", help="input: bam file", type="character")
p <- add_argument(p, "--gtf", help="input: gtf file", type="character")
p <- add_argument(p, "--featureType", help="a character string or a vector of character strings giving the feature type or types used to select rows in the GTF annotation which will be used for read summarization", type="character", default="exon")
p <- add_argument(p, "--attrType", help="a character string giving the attribute type in the GTF annotation which will be used to group features (eg. exons) into meta-features", type="character", default="gene_id")
p <- add_argument(p, "--isPairedEnd", help="indicating whether libraries contain paired-end reads or not", type="logical", default=TRUE)
p <- add_argument(p, "--strandSpecific", help="0 (unstranded), 1 (stranded) and 2 (reversely stranded)", type="numeric", default=0)
p <- add_argument(p, "--output", help="output prefix", type="character")

# Parse the command line arguments
argv <- parse_args(p)

bamFile <- argv$bam
gtfFile <- argv$gtf
nthreads <- 1
outFilePref <- argv$output

outStatsFilePath  <- paste(outFilePref, '.log',  sep = ''); 
outCountsFilePath <- paste(outFilePref, '.count', sep = ''); 

fCountsList = featureCounts(bamFile, annot.ext=gtfFile, isGTFAnnotationFile=TRUE, nthreads=nthreads,
                            GTF.featureType=argv$featureType, GTF.attrType=argv$attrType, isPairedEnd=argv$isPairedEnd,
                            strandSpecific=argv$strandSpecific, countMultiMappingReads=FALSE)
dgeList = DGEList(counts=fCountsList$counts, genes=fCountsList$annotation)
cpm = cpm(dgeList)
fpkm = rpkm(dgeList, dgeList$genes$Length)
tpm = exp(log(fpkm) - log(sum(fpkm)) + log(1e6))

write.table(fCountsList$stat, outStatsFilePath, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)

featureCounts = cbind(fCountsList$annotation[,1], fCountsList$counts, fpkm, tpm, cpm)
colnames(featureCounts) = c('gene_id', 'counts', 'fpkm','tpm', 'cpm')
write.table(featureCounts, outCountsFilePath, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
