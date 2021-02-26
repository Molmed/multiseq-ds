if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()

library(Seurat)
library(Signac, lib.loc = "/domus/h1/alvaa/R")

#--- args from comandline ----------

#args[1] = output file name 
#args[2:] = path/to/seurat objects 

args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("provide path seurat objects to merge", call.=FALSE)
}else if(length(args)==1){
  args[1] = "seurat_out.rds"
} 

#--- MERGE DATA ---------------------

obj = list()

print("READ DATA")

for(i in 2:length(args)){
  obj[[i-1]] = readRDS(args[i])
  obj[[i-1]]$dataset = args[i]
}

print("MERGE")

combined = merge(x = obj[[1]], y = obj[2:length(obj)])

#--- SAVE -----------------------------

saveRDS(combined, file=args[1])

