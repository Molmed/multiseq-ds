#script to create common peak set  

library(Signac, lib.loc = "/domus/h1/alvaa/R")
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
set.seed(1234) 

#--- args from command line -----

#args[1] = output_file_name
#args[2] = /path/to/output folder form count

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("provide path to dir output folder form count and output file name", call.=FALSE)
}

#--- PEAK BED FILES -------------

bed.REH = paste(args[2], "/outs/peaks.bed", sep="")
peaks.REH = makeGRangesFromDataFrame((read.table(file = bed.REH, col.names = c("chr", "start", "end"))))

bed.GM = paste(args[3], "/outs/peaks.bed", sep="")
peaks.GM = makeGRangesFromDataFrame((read.table(file = bed.GM, col.names = c("chr", "start", "end"))))

print("READ PEAKS OK")

#--- COMBINE PEAKS --------------

combined.peaks = reduce(x = c(peaks.REH, peaks.GM))

print("COMBINE PEAKS OK")

#--- FILTER ---------------------

peakwidths = width(combined.peaks)
combined.peaks = combined.peaks[peakwidths  < 10000 & peakwidths > 20]

#--- DATA -----------------------

path.singlecell.REH = path_frag =  paste(args[2], "/outs/singlecell.csv", sep="")
path.frag.REH =  paste(args[2], "/outs/fragments.tsv.gz", sep="")

metadata.REH = read.table(
  file = path.singlecell.REH ,
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1,]

path.singlecell.GM = path_frag =  paste(args[3], "/outs/singlecell.csv", sep="")
path.frag.GM =  paste(args[3], "/outs/fragments.tsv.gz", sep="")

metadata.GM = read.table(
  file = path.singlecell.GM ,
  stringsAsFactors = FALSE,
  sep = ",",
  header = TRUE,
  row.names = 1
)[-1,]

frag.REH = CreateFragmentObject(
  path = path.frag.REH,
  cells = rownames(metadata.REH),
  validate.fragments = FALSE)

frag.GM = CreateFragmentObject(
  path = path.frag.GM,
  cells = rownames(metadata.GM),
  validate.fragments = FALSE)

print("DATA OK")

#--- QUANTIFY PEAKS -------------

counts.REH = FeatureMatrix(
  fragments = frag.REH,
  features = combined.peaks,
  cells = rownames(metadata.REH)
)

counts.GM = FeatureMatrix(
  fragments = frag.GM,
  features = combined.peaks,
  cells = rownames(metadata.GM)
)

#--- SEURAT OBJECT --------------

chr.assay.REH = CreateChromatinAssay(counts.REH, fragments = frag.REH)
seurat.REH = CreateSeuratObject(chr.assay.REH, assay = "ATAC")
seurat.REH$dataset = "REH"
print(seurat.REH)

chr.assay.GM = CreateChromatinAssay(counts.GM, fragments = frag.GM)
seurat.GM = CreateSeuratObject(chr.assay.GM, assay = "ATAC")
seurat.GM$dataset = "GM12078"
print(seurat.GM)

print("SEURAT OBJECT OK")

#--- MERGE ---------------------

combined <- merge(
  x = seurat.REH,
  y = seurat.GM,
  add.cell.ids = c("REH", "GM12078")
)

print("MERGE OK")

#--- SAVE -----------------------

saveRDS(combined, file = args[1])


