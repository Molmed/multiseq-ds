
# Script to preform aditional preprocessing on output from cellranger-atac count and aggr. 
# Calculate QC and gene activity 
# Saves seurat object as RDS file  

#------------------------------------
#Install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

#install.packages("Signac", repos = "http://cran.us.r-project.org", lib = "/domus/h1/alvaa/R")

library(Seurat)
library(Signac, lib.loc = "/domus/h1/alvaa/R")
library(ggplot2)
library(patchwork)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
set.seed(1234) 

#--- args from comandline ----------

#args[1] = absolute/path/to/output folder form count/aggr 
#args[2] = output_file_name

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("provide path to dir output folder form count/aggr", call.=FALSE)
}else if(length(args)==1){
  args[1] = "seurat_out"
}

path_feat_matrix =  paste(args[1], "/outs/filtered_peak_bc_matrix.h5", sep="")
path_singlecell =  paste(args[1], "/outs/singlecell.csv", sep="")
path_frag =  paste(args[1], "/outs/fragments.tsv.gz", sep="")

#--- DATA --------------------------

counts <- Read10X_h5(filename = path_feat_matrix)

metadata <- read.csv(file = path_singlecell, header = TRUE, row.names = 1)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'GRCh38',
  fragments = path_frag,
  min.cells = 10,
  min.features = 200
)

#create seurat object 
data <- CreateSeuratObject(
  counts = chrom_assay, 
  assay = "peaks", 
  meta.data = metadata
)

print("DATA OK")

#--- ANNOTATION --------------------

#gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "GRCh38"

#change to UCSC style 
seqlevelsStyle(annotations) <- 'UCSC'

#add the gene information to object
Annotation(data) <- annotations

print("ANNOTATION OK")

#--- CALCULATE QC ----------------------------
#1. Nucleosome banding pattern
#2. Transcriptional start site (TSS) enrichment score 
#3. Total number of fragments in peaks
#4. Fraction of fragments in peaks
#5. Ratio reads in genomic blacklist regions

#compute nucleosome signal score per cell (Nucleosome banding pattern)
data <- NucleosomeSignal(object = data)

#compute TSS enrichment score per cell
data <- TSSEnrichment(object = data, fast = FALSE)

#fraction of reads in peaks
data$pct_reads_in_peaks <- data$peak_region_fragments / data$passed_filters * 100

#blacklist ratio 
data$blacklist_ratio <- data$blacklist_region_fragments / data$peak_region_fragments

print("QC OK")

#--- CALCULATE GENE ACTIVITY -------

gene.activities <- GeneActivity(data)

#add the gene activity matrix to the Seurat object as a new assay and normalize it
data[['gene.act']] <- CreateAssayObject(counts = gene.activities)

print("GENE ACT OK")

#--- SAVE SEURAT OBJECT -----------

saveRDS(data, file = args[2])

