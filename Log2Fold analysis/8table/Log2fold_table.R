# Install and load R pakcages that are necessary for the analysis - Packages are collections of R functions, data, 
# and compiled code in a well-defined format. Remove the hash sign to download and install the packages.

#source('http://bioconductor.org/biocLite.R')
#biocLite('phyloseq')
library("phyloseq")
packageVersion("phyloseq")

#biocLite("biomformat")
library("biomformat")
packageVersion("biomformat")

#install.packages("ggplot2")
library("ggplot2")
packageVersion("ggplot2")

#install.packages("vegan")
library("vegan")
packageVersion('vegan')

#install.packages("grid")
library("grid")
packageVersion('grid')

#install.packages("magrittr")
library(magrittr)
packageVersion('magrittr')

library(dplyr)
packageVersion('dplyr')

library('DESeq2')
packageVersion('DESeq2')

library("phyloseq")

#Make sure the data has the sample names on rows
abund_table<-read.csv("NicTable_otu.csv",row.names=1,check.names=FALSE)
meta_table<-read.csv("table_meta.csv",row.names=1,check.names=FALSE)
#Transpose the data to have sample names on rows
abund_table<-t(abund_table)

meta_table<-data.frame(meta_table)
#Filter out samples not present in meta_table
abund_table<-abund_table[rownames(abund_table) %in% rownames(meta_table),]

#Now load the taxonomy
OTU_taxonomy<-read.csv("taxonomy_phyloseq2.csv",row.names=1,check.names=FALSE)

OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)
physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM)

###for Box and whisker plot###
#bnw=tax_glom(physeq, "Genus" )
#OTU1 = as(otu_table(bnw), "matrix")
#OTUdf = as.data.frame(OTU1)
#OTUdf<-t(OTUdf)
#write.csv(OTUdf, file = "GenusOTU.csv")

#ANALYSIS
diagdds = phyloseq_to_deseq2(physeq, ~ category)
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="local")
res = results(diagdds)
res = res[order(res$padj, na.last=NA), ]
alpha = 0.01
sigtab = res[(res$padj < alpha),]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
#sigtab

library("ggplot2")
pdf("table_0T_N_5_7_20.pdf")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
sigtabgen2= sigtabgen[order(sigtab$log2FoldChange < -7 && sigtab$log2FoldChange > 7,decreasing = TRUE)]
# # Phylum order
# x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
# x = sort(x, TRUE)
# sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))
sigtabgen2=sigtabgen[abs(sigtabgen$log2FoldChange) > 6.0,]
sigtabgen3= sigtabgen2[order(sigtabgen2$log2FoldChange,decreasing = TRUE),]

# mapping = aes(x = reorder(Name, Number), Number)

ggplot(sigtabgen3, aes(y=Genus, x=log2FoldChange))+
  ggtitle("Table Microbiome: Smoking/Non-smoking") +
  geom_vline(xintercept = 0.0, color = "gray", size = 0.01) +
  xlab(bquote(~Log[2]~'( fold change )' ))+
  geom_point(size=5) + theme_bw()+
  theme(text= element_text(size=10),axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5),legend.position = "none")

dev.off()



positive= sigtabgen[sigtabgen$log2FoldChange>0,]
negative= sigtabgen[sigtabgen$log2FoldChange<0,]
pGcounts=summary(positive['Genus'],maxsum = 40)
nGcounts=summary(negative['Genus'],maxsum = 40)
write.table(pGcounts,"GCounts2.csv",col.names=F , row.names = FALSE)
write.table( nGcounts,  
             file="./GCounts2.csv", 
             append = T, 
             sep=',', 
             row.names=T, 
             col.names=F )
# 
# ###Histogram###
# top25otus = names(sort(taxa_sums(physeq), TRUE)[1:25])
# taxtab25 = cbind(tax_table(physeq), family25 = NA)
# taxtab25[top25otus, "family25"] <- as(tax_table(physeq)[top25otus, "Family"], 
#                                       "character")
# tax_table(physeq) <- tax_table(taxtab25)
# 
# title = "Table: top OTUs; Family taxa; 01-019=smoker; 101-110=non-smoker"
# #plot_bar(physeq, "Sample", fill = "family19", title = title) + coord_flip()
# p19 = prune_taxa(top25otus, physeq)
# plot_bar(p19, "Sample", fill = "family25", title = title) + coord_flip()
# 
# 
