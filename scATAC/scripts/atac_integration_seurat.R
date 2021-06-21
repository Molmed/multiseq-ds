
library(Seurat)
library(Signac)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
set.seed(1234) 

#--- args from command line -----

#args[1] = output_file_name
#args[2] = path to output from merge combined peak set 

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("provide path to output from atac_merge_combined_peaks.R and output file name", call.=FALSE)
}

#--- DATA -------------------
 
data = readRDS(args[2])
print("DATA OK") 
print(data)  

DefaultAssay(data) = "ATAC"
Idents(data) = "dataset"
sets = levels(Idents(data))

print(DefaultAssay(data))
print(Idents(data))
print(sets)

obj = list()

#make a list of seurat objects for each dataset 
for(i in 1:length(sets)){
   obj[i] = subset(data, idents = sets[i])
}

print("LIST OBJECTS OK")

#--- FIND ANCHORS -----------

anchors = FindIntegrationAnchors(
  object.list = obj,
  anchor.features = rownames(data),
  k.filter = NA
)

saveRDS(anchors, "anchors.rds")

print("ANCHORS OK") 

#--- INTEGRATE --------------

integrated = IntegrateData(
  anchorset = anchors,
  preserve.order = TRUE
)

print("INTEGRATE OK") 

#--- SAVE DATA --------------

saveRDS(integrated, args[2])


