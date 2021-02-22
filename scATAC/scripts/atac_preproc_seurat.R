
# Script to preform aditional preprocessing and filtering 
# on output from cellranger-atac count and aggr. 
# 1. Calculate QC 
# 2. Filter on QC 
# 3. Normalize 
# 4. Dim redction 

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
  stop("provide path to dir outs", call.=FALSE)
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

#--- QC ----------------------------
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

#--- PLOT QC -----------------------

v_plot = VlnPlot(
  object = data,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5,
  cols = "orange"
)

png(file="QC_plot", width=600, height=350)
v_plot
dev.off()

#--- FILTER QC ---------------------

data <- subset(
  x = data,
  subset = peak_region_fragments > 3000 &
  peak_region_fragments < 20000 &
  pct_reads_in_peaks > 15 &
  blacklist_ratio < 0.05 &
  nucleosome_signal < 4 &
  TSS.enrichment > 2
)

print("FILTER OK")

#--- PLOT FILTERED -----------------

v_plot = VlnPlot(
  object = data,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5,
  cols = "gold"
)

png(file="QC_plot_filtered", width=600, height=350)
v_plot
dev.off()

#--- NORMALIZATION ----------------

data <- RunTFIDF(data)
data <- FindTopFeatures(data, min.cutoff = 'q0')

print("NORMALIZATION OK")

#--- DIM REDUCTION ----------------

data <- RunSVD(data)

print("DIM REDUCTION OK")

#--- SAVE SEURAT OBJECT -----------

saveRDS(data, file = args[2])

