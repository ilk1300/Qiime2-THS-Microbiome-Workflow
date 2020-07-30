# # Install and load R pakcages that are necessary for the analysis - Packages are collections of R functions, data, 
# # and compiled code in a well-defined format. Remove the hash sign to download and install the packages.
# 
# #source('http://bioconductor.org/biocLite.R')
# #biocLite('phyloseq')
 library("phyloseq")
# packageVersion("phyloseq")
# 
# #biocLite("biomformat")
# library("biomformat")
# packageVersion("biomformat")
# 
# #install.packages("ggplot2")
 library("ggplot2")
# packageVersion("ggplot2")
# 
# #install.packages("vegan")
# library("vegan")
# packageVersion('vegan')
# 
# #install.packages("grid")
# library("grid")
# packageVersion('grid')
# 
# #install.packages("magrittr")
# library(magrittr)
# packageVersion('magrittr')
# 
# library(dplyr)
# packageVersion('dplyr')
# 
library('DESeq2')
# packageVersion('DESeq2')
# 
library("phyloseq")

#Make sure the data has the sample names on rows

# 
# location= readline(prompt="Please enter microbial community location: " )
# 
# if (location="finger") {
#       #1. Finger
#       otu = "new_finger_otu.csv" 
#       meta = "Finger_metadataM.csv"
#       }
# if (location="ear"){
#       #2. Ear
#       otu="ear_otuTrimmed.csv"
#       meta = "ear_metaM.csv"
# }
# 
# if(location="nose"){
#       #3. Nose
#       otu= "new_nose_otu.csv "
#       meta="nose_meta.csv"
# }     
# if(location="mouth"){
#       #4. Mouth
#       otu= "new_mouth_otu1.csv"
#       meta="mouth_metaM.csv"
# }
# if(location="pillowcase"){
#       #5. Pillowcase
#       otu= "pillowcase_otu.csv"
#       meta= "pillow_metaM.csv"
# }
# if(location="bedframe"){
#       #6. Bed framed
#   otu="bedframe_otu.csv"
#   meta="bedframe_metaM.csv"
# }
# if(location="bottomsheet"){
#   #7. Bottom Sheet
#   otu="bottom_otuT.csv"
#   meta="bottom_sheet_meta.csv"
# }
# if(location="table"){
#   #8. Table
#   otu="table_otuT.csv"
#   meta="table_meta.csv"
# }
# if(location="armrest"){
#   #9. Arm rest
#   otu="arm_otuT.csv"
#   meta="armRest_meta.csv"
# }
# if(location="floor"){
#   #10. floor
#   otu="floor_otuT.csv"
#   meta="floor_metaMild.csv"
# }
#       
#       
 otu = "FingerOTU5_15.csv" 
 meta = "Finger_meta_Nic.csv"
# 

#abund_table<-read.csv("finger_otu.csv",row.names=1,check.names=FALSE)
abund_table<-read.csv(otu,row.names=1,check.names=FALSE)

meta_table<-read.csv(meta,row.names=1,check.names=FALSE)
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
alpha = 0.05
sigtab = res[(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(physeq)[rownames(sigtab), ], "matrix"))
#sigtab

library("ggplot2")
pdf("FingerReducedM.pdf")
theme_set(theme_bw())
sigtabgen = subset(sigtab, !is.na(Genus))
x = tapply(sigtabgen$log2FoldChange, sigtabgen$Genus, function(x) max(x))
x = sort(x, FALSE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels=names(x))

sigtabgen2=sigtabgen[abs(sigtabgen$log2FoldChange) > 6.0,]
sigtabgen3= sigtabgen2[order(sigtabgen2$log2FoldChange,decreasing = TRUE),]

#plot_bar(subsetA,"category","Abundance" )

#newdata <- sigtabgen[order(sigtabgen$log2FoldChange),]
#sigtabgen2= sigtab[order(sigtab$log2FoldChange < -2 & sigtab$log2FoldChange > 2,decreasing = TRUE)]
# Phylum order
#x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
#x = sort(x, TRUE)
#sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
# Genus order

ggplot(sigtabgen3, aes(y=Genus, x=log2FoldChange))+
  ggtitle("Finger Microbiome") +
  geom_vline(xintercept = 0.0, color = "gray", size = 1) +
  #geom_point(size=8, colour="black") + 
  geom_point(position=position_jitter(h=0.1, w=0.1),
             shape = 21, alpha = 0.5, size = 8, fill="black")+
  scale_color_manual(aesthetics = c("colour", "fill"))+
  theme_bw(base_size = 18)+
  labs(x=expression("Ratio:"~"THS-polluted/THS-free"~"Homes"~(log[2])  ))+
  #xlab("Ratio: THS-polluted vs THS-free Homes Log[2]") +
  theme(
    plot.title = element_text(size = 16, face = "bold",hjust=0.5),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15, angle = 90),
    axis.text.x = element_text( vjust = 1, size=14),
    axis.text.y = element_text(face = "italic",size=14), 
    legend.position = "none")
#x=expression(Production~rate~" "~mu~moles~NO[3]^{"-"}-N~Kg^{-1})
#theme(axis.text.y = blue.bold.italic.16.text)
#xlab(bquote(~Log[2]~'( fold change )' ))
dev.off()




##PERMANOVA##
library(vegan)
library("ape")
metadata <- as(sample_data(physeq), "data.frame")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
#plot(random_tree)
physeq2 = phyloseq(OTU, TAX, SAM, random_tree)
unifrac.dist.finger <- UniFrac(physeq2, 
                        weighted = FALSE, 
                        normalized = TRUE,  
                        parallel = FALSE, 
                        fast = TRUE)
metadf <- data.frame(sample_data(physeq))
permanova <- adonis(unifrac.dist.finger ~ category, data = metadf,permutations = 9999)

permanova

ps.disper <- betadisper(unifrac.dist, metadf$category)
permutest(ps.disper, pairwise = TRUE)



# 
# ###option 2###
#   top25otus = names(sort(taxa_sums(physeq), TRUE)[1:25])
#   taxtab25 = cbind(tax_table(physeq), family25 = NA)
#   taxtab25[top25otus, "family25"] <- as(tax_table(physeq)[top25otus, "Family"], 
#                                         "character")
#   tax_table(physeq) <- tax_table(taxtab25)
# 
#   title = "Finger: top OTUs; Family taxa; 01-019=smoker; 101-110=non-smoker"
#   #plot_bar(physeq, "Sample", fill = "family19", title = title) + coord_flip()
#   p19 = prune_taxa(top25otus, physeq)
#   plot_bar(p19, "Sample", fill = "family25", title = title) + coord_flip()
#   
#   restroomRm = merge_samples(restroomR, "SURFACE")
#   
#   ##option3###
#   data[data$phylum %in% remainder,]$Phylum <- "Phyla < 1% abund."
#   #rename phyla with < 1% relative abundance
#   data$phylum[data$Abundance < 0.01] <- "Phyla < 1% abund."
#   
#   #plot with condensed phyla into "< 1% abund" category
#   p <- ggplot(data=data, aes(x=Sample, y=Abundance, fill=phylum))
#   p + geom_bar(aes(), stat="identity", position="stack") +
#     scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue",
#                                  "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
#     theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
#   
#   
#   p <- ggplot(data=top25otus, aes(x=Sample, y=Abundance, fill=family))
#   p + geom_bar(aes(), stat="identity", position="stack") +
#     scale_fill_manual(values = c("darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey", "darkblue", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "lightskyblue", "darkgreen", "deeppink", "khaki2", "firebrick", "brown1", "darkorange1", "cyan1", "royalblue4", "darksalmon", "darkblue", "royalblue4", "dodgerblue3", "steelblue1", "lightskyblue", "darkseagreen", "darkgoldenrod1", "darkseagreen", "darkorchid", "darkolivegreen1", "brown1", "darkorange1", "cyan1", "darkgrey")) +
#     theme(legend.position="bottom") + guides(fill=guide_legend(nrow=5))
#   
#   
# ###Random Forest###
#   library(randomForest)
#   library(caret)
#   abund_table<-read.csv("new_finger_otu.csv",row.names=1,check.names=FALSE)
#   abund_table<-t(abund_table)
#   meta_table<-as.data.frame(abund_table)
#   
#   set.seed(123)
#   ind1 <- sample(2,nrow(meta_table),replace = TRUE, prob = c(0.83,0.17))
#   train1 <-meta_table[ind1==1,]
#   test1 <- meta_table[ind1==2,]
#   
#   rf1 <- randomForest(Target~.,data=train1,ntree= 460,mtry = 15,importance=TRUE,proximity=TRUE)
#   prediction1 <- predict(rf1,test1)
#   
#   plot(rf1)
#   tuning = tuneRF(train1[,-15],train1[,15],stepFactor = 0.8,plot = TRUE,ntreeTry = 460,trace = TRUE,improve = 0.05)
#   varImpPlot(rf1,sort = T,n.var = 5,main = "Feature Importance")
#   print("No. of Trees 460 and split ration of 0.83:")
#   print(confusionMatrix(prediction1,abund_table$Target))
#   
#   plot(varImp(rf1), top = 10)
#   tm = rpart(Target~., train1, method="class")
#   library(rpart.plot)
#   
#   testing <- as.factor(test1$class)
#   
#   boxplot(X~ faith_pd*X.1 , data=test, notch=TRUE,
#           col=(c("gold","darkgreen")),
#           main="faith-pd", xlab="Samples")
#   
#   
#     