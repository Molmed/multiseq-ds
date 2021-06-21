#!/usr/bin/env Rscript

#USE: 
# Rscript --vanilla demultiplex_seurat.R \
# "path/to/filtered_feature_bc_matrix" #rna \
# "path/to/filtered_feature_bc_matrix" #hto\
# "output_file_name"

library(Seurat)

#--- args from comandline ----------

#args[1] = path/to/filtered_feature_bc_matrix for RNA and HTO 
#args[2] = output_file_name

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("provide paths to RNA and HTO feature matrix", call.=FALSE)
}else if(length(args)==1){
  args[3] = "out"
}

dir_data = paste(args[1], "/outs/filtered_feature_bc_matrix", sep = "")

print(dir_data)
print(args[2])

#--- DATA ----------------------

#read umi matrix 
data.count = Read10X(data.dir = dir_data)

print("DATA OK")

#--- CELL BARCODES -------------

#select cell barcodes found in both gex and ht matrix  
joint = intersect(colnames(data.count$`Gene Expression`), colnames(data.count$Custom))

#subset on joint cell barcodes
rna.data = data.count$`Gene Expression`[,joint]
hto.data = data.count$Custom[,joint]

print("BARCODES OK")

#--- SEURAT OBJECT -----------

#create seurat object with RNA assay 
data = CreateSeuratObject(counts = rna.data)

#add hashtag data as a new assay 
data[["HTO"]] = CreateAssayObject(counts = hto.data)

#normalize HTO data (centered log-ratio (CLR) transformation) 
data = NormalizeData(data, assay = "HTO", normalization.method = "CLR")

print("SEURAT OBJECT OK")

#--- DEMULTIPLEX --------------

#samples can be found in data@meta.data

data = HTODemux(data, assay = "HTO", positive.quantile = 0.90)

print(data)

print("DEMULTIPLEX OK")

#--- SAVE ---------------------

#save seurat object as RDS file 
saveRDS(data, file = args[2])

print("DONE")
 


