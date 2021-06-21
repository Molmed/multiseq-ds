
library(deMULTIplex, lib.loc = "/domus/h1/alvaa/R")
library(ggplot2)
library(readr)
library(Seurat)


#--- DATA ----------------------

#read umi matrix
data.count = Read10X(data.dir = "/crex/proj/uppstore2017165/nobackup/alva/REH_MULTI-seq/out/FU180-MULTI-seq/FU180-MULTI-seq_count/outs/filtered_feature_bc_matrix")

print("DATA OK")

#--- CELL BARCODES -------------

joint = intersect(colnames(data.count$`Gene Expression`), colnames(data.count$Custom))

#select cell barcodes found in both gex and ht matrix
bar.table = as.data.frame(t(as.data.frame(data.count$Custom[,joint])))

#--- BARCODE TSNE -------------

bar.tsne = barTSNE(bar.table)

pdf("bc.check.pdf")
for (i in 3:ncol(bar.tsne)) {
  g = ggplot(bar.tsne, aes(x = TSNE1, y = TSNE2, color = bar.tsne[,i])) +
    geom_point() +
    scale_color_gradient(low = "black", high = "red") +
    ggtitle(colnames(bar.tsne)[i]) +
    theme(legend.position = "none") 
  print(g)
}
dev.off()


#--- THRESHOLD --------------------

bar.table_sweep.list <- list()
n <- 0
for (q in seq(0.01, 0.99, by=0.02)) {
  print(q)
  n <- n + 1
  bar.table_sweep.list[[n]] <- classifyCells(bar.table, q=q)
  names(bar.table_sweep.list)[n] <- paste("q=",q,sep="")
}

#Identify ideal inter-maxima quantile to set barcode-specific thresholds
threshold.results1 <- findThresh(call.list=bar.table_sweep.list)

pdf("threshold.check.pdf")
thresh = ggplot(data=threshold.results1$res, aes(x=q, y=Proportion, color=Subset)) + 
  geom_line() + 
  theme(legend.position = "none") +
  geom_vline(xintercept=threshold.results1$extrema, lty=2) + 
  scale_color_manual(values=c("red","black","blue"))
print(thresh)
dev.off()

#--- CLASSIFICATIONS --------------

#This step might need to be repeated multiple times until no negative cells are present 
#Write loop to do this automatically (if better data available)

round1.calls <- classifyCells(bar.table, q=findQ(threshold.results1$res, threshold.results1$extrema))
neg.cells <- names(round1.calls)[which(round1.calls == "Negative")]
bar.table <- bar.table[-which(rownames(bar.table) %in% neg.cells), ]

#--- FINAL CLASSIFICATION -----------------

final.calls <- c(round1.calls, rep("Negative",length(neg.cells)))
names(final.calls) <- c(names(round1.calls),neg.cells)

#--- SAVE ----------------------------------

saveRDS(final.calls, "cell.classificaion.MULTI-seq.rds")



