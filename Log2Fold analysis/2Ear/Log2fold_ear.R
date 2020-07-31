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
#abund_table<-read.csv("earReduced1.csv",row.names=1,check.names=FALSE)
abund_table<-read.csv("earNict5_2_20v2.csv",row.names=1,check.names=FALSE)

meta_table<-read.csv("ear_meta_Nic.csv",row.names=1,check.names=FALSE)

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
#t=otu_table(physeq)

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
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
#sigtab

#put in order, absolute value of 2 in log2Fold.
#t=sigtab[abs(sigtab$log2FoldChange)>2,]


#PLOT GRAPH & SAVE

library("ggplot2")
pdf("Ear0T_N_5_2_20.pdf")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

sigtabgen2=sigtabgen[abs(sigtabgen$log2FoldChange) > 6.0,]
sigtabgen3= sigtabgen2[order(sigtabgen2$log2FoldChange,decreasing = TRUE),]

ggplot(sigtabgen3, aes(y=Genus, x=log2FoldChange))+
  ggtitle("Ear Microbiome: Smoking/NonSmoking") +
  geom_vline(xintercept = 0.0, color = "gray", size = 1) +
  geom_point(size=8, colour="black") + 
  # geom_point(position=position_jitter(h=0.1, w=0.1),
  #            shape = 21, alpha = 0.5, size = 8, fill="black")+
  #scale_color_manual(aesthetics = c("colour", "fill"))+
  theme_bw(base_size = 18)+
  labs(x=expression(Ratio:THS-polluted~vs~THS-free~Homes~'(log)'[2] ))+
  #xlab("Ratio: THS-polluted vs THS-free Homes Log[2]") +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15, angle = 90),
    axis.text.x = element_text( vjust = 1, size=14),
    axis.text.y = element_text(face = "italic",size=14), 
    legend.position = "none")
#x=expression(Production~rate~" "~mu~moles~NO[3]^{"-"}-N~Kg^{-1})
#theme(axis.text.y = blue.bold.italic.16.text)
#xlab(bquote(~Log[2]~'( fold change )' ))
dev.off()


positive= sigtabgen3[sigtabgen3$log2FoldChange>0,]
negative= sigtabgen3[sigtabgen3$log2FoldChange<0,]
pGcounts=summary(positive['Genus'],maxsum = 40)
nGcounts=summary(negative['Genus'],maxsum = 40)
write.table(pGcounts,"GCounts2.csv",col.names=F , row.names = FALSE)
write.table( nGcounts,  
             file="./GCounts2.csv", 
             append = T, 
             sep=',', 
             row.names=T, 
             col.names=F )

##PERMANOVA##
library(vegan)
library("ape")
metadata <- as(sample_data(physeq), "data.frame")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
#plot(random_tree)
physeq2 = phyloseq(OTU, TAX, SAM, random_tree)
unifrac.dist.ear <- UniFrac(physeq2, 
                               weighted = FALSE, 
                               normalized = TRUE,  
                               parallel = FALSE, 
                               fast = TRUE)
metadf <- data.frame(sample_data(physeq))
permanova <- adonis(unifrac.dist.ear ~ category, data = metadf)

permanova

ps.disper <- betadisper(unifrac.dist.ear, metadf$category)
permutest(ps.disper, pairwise = TRUE)


# 
# figure <- ggarrange(sp, bp + font("x.text", size = 10),
#                     ncol = 1, nrow = 2)
# annotate_figure(figure,
#                 top = text_grob("Visualizing mpg", color = "red", face = "bold", size = 14),
#                 bottom = text_grob("Data source: \n mtcars data set", color = "blue",
#                                    hjust = 1, x = 1, face = "italic", size = 10),
#                 left = text_grob("Figure arranged using ggpubr", color = "green", rot = 90),
#                 right = "I'm done, thanks :-)!",
#                 fig.lab = "Figure 1", fig.lab.face = "bold"
# )
# 
# 
# ###HISTOGRAM GENERATETOR###
# top25otus = names(sort(taxa_sums(physeq), TRUE)[1:25])
# taxtab25 = cbind(tax_table(physeq), family25 = NA)
# taxtab25[top25otus, "Family"] <- as(tax_table(physeq)[top25otus, "Family"], 
#                                       "character")
# tax_table(physeq) <- tax_table(taxtab25)
# title = "Ear: top 25 genera with family taxa colored; 01-019=smoker; 101-110=non-smoker"
# #plot_bar(physeq, "Sample", fill = "family25", title = title) + coord_flip()
# p19 = prune_taxa(top25otus, physeq)
# plot_bar(p19, "Sample", fill = "family25", title = title) + coord_flip()
# 
# erie_phylum <- erie %>%
#   tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
#   transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
#   psmelt() %>%                                         # Melt to long format
#   filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
#   arrange(Phylum)                                      # Sort data frame alphabetically by phylum
# library(tidyverse)
# library(lubridate)
# #`df` is `CPOD_Time` saved as `df<-as.data.frame(CPOD_Time)`
# df<-as.data.frame(CPOD_Time)


###Random Forest###
# library(randomForest)
# library(caret)
# abund_table<-read.csv("MyData.csv",row.names=1,check.names=FALSE)
# #abund_table<-t(abund_table)
# meta_table<-data.frame(abund_table)
# 
# set.seed(123)
# ind1 <- sample(2,nrow(meta_table),replace = TRUE, prob = c(0.83,0.17))
# train1 <-meta_table[ind1==1,]
# test1 <- meta_table[ind1==2,]
# 
# rf1 <- randomForest(Target~.,data=train1,ntree= 460,mtry = 150,importance=TRUE,proximity=TRUE)
# prediction1 <- predict(rf1,test1)
# 
# plot(rf1)
# tuning = tuneRF(train1[,-15],train1[,15],stepFactor = 0.8,plot = TRUE,ntreeTry = 460,trace = TRUE,improve = 0.05)
# varImpPlot(rf1,sort = T,n.var = 5,main = "Feature Importance")
# print("No. of Trees 460 and split ration of 0.83:")
# print(confusionMatrix(prediction1,abund_table$Target))
# 
# plot(varImp(rf1), top = 10)
# tm = rpart(Target~., train1, method="class")
# library(rpart.plot)
# 
