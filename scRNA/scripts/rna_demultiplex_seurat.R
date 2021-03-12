#!/usr/bin/env Rscript

#USE: 
# Rscript --vanilla demultiplex_seurat.R \
# "path/to/filtered_feature_bc_matrix" #rna \
# "path/to/filtered_feature_bc_matrix" #hto\
# "output_file_name"

library(Seurat)

#--- args from comandline ----------

#args[1] = path/to/filtered_feature_bc_matrix for RNA
#args[2] = path/to/filtered_feature_bc_matrix for HTO
#args[3] = output_file_name

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("provide paths to RNA and HTO feature matrix", call.=FALSE)
}else if(length(args)==2){
  args[3] = "out"
}

dir_rna = paste(args[1], "/outs/filtered_feature_bc_matrix", sep = "")
dir_hto = paste(args[2], "/outs/filtered_feature_bc_matrix", sep = "")

print(dir_rna)
print(dir_hto)

#--- DATA ----------------------

#read umi matrix 
rna.data = Read10X(data.dir = dir_rna)

#read hashtag matrix 
hto.data = Read10X(data.dir = dir_hto)

#--- CELL BARCODES -------------

#select cell barcodes found in both gex and ht matrix  
joint = intersect(colnames(rna.data), colnames(hto.data))

#subset on joint cell barcodes
rna.data = rna.data[, joint]
hto.data = hto.data[, joint]

#--- SEURAT OBJECT -----------

#create seurat object with RNA assay 
data = CreateSeuratObject(counts = rna.data)

#add hashtag data as a new assay 
data[["HTO"]] = CreateAssayObject(counts = hto.data)

#normalize HTO data (centered log-ratio (CLR) transformation) 
data = NormalizeData(data, assay = "HTO", normalization.method = "CLR")

#--- DEMULTIPLEX --------------

#samples can be found in data@meta.data

data = HTODemux(data, assay = "HTO", positive.quantile = 0.99)

#--- SAVE ---------------------

#save seurat object as RDS file 
saveRDS(data, file = args[3])



