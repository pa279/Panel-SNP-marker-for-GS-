####selectcion of marker agre
setwd("/Users/user/Desktop/Panel_marker/")
library(tidyverse)
library(adegenet)
library(viridis)
library(VariantAnnotation)
library(Rsamtools)
library(ggplot2)
library(ape)
library(dplyr)
library(cluster)
library(arules)
library(snpReady)
library(digest)
library(Matrix)
library(gdsfmt)
library(crayon)
library(RUnit)
library(knitr)
library(markdown)
library(BiocGenerics)
library("devtools")
library(SNPRelate)
library(qgraph)
library(gplots)
library(lattice)
library(FactoMineR)
library(factoextra)
library(mice)
library(missMDA)
library(ade4)
library(VIM)
library(caret)
library(lme4)
library(psych)
library(ppcor)
library(missForest)
library(philentropy)
library(proxy)
library(arules)
library(vegan)
library(ggfortify)
library(rrBLUP)
library(ggplot2)
library(cowplot)
library(apcluster)
library(corrplot)
library(lattice)
library(ppcor)
library(viridis)
library(dendextend)
library(ggsignif)
library(qgraph)
library(FactoMineR)
library(VIM)
library(GGally)
library(cluster)
library(graphics)
library(philentropy)
library(proxy)
library(biotools)
library(rpanel)
library(tcltk)
library(NbClust)
library(GenAlgo)
library(stats)
library(StatMatch)
library(ClassDiscovery)
library(xtable)
source("marker_stats.R")
source("getNumGeno.R")
source("evaluate_marker_sets.R")

bg <- "331_new_data_vcf_IBRC.vcf"
rds <- "yam336.m2M2vsnps_missing0.9.recode.rds" # unfiltered markers; optional
refgenome <- FaFile("TDr96_F1_v2_PseudoChromosome.rev07.fasta")
numgen <- getNumGeno(bg)
numgen[1:10,1:10]
??interIndividualIBS
mydist <- interIndividualIBS(numgen)

hist(mydist, xlab = "Euclidian distance", main = "Histogram of inter-individual distances")
par(cex=0.7, font =1)
plot(nj(mydist), type = "unrooted", show.tip.label = F)
numgen <- removeClones(numgen, mydist)
alfreq <- rowMeans(numgen, na.rm = TRUE) / 2
hist(alfreq, xlab = "Allele frequency", main = "Histogram of alternative allele frequency")
He <- Expected_het(alfreq)
Ho <- rowMeans(numgen == 1, na.rm = TRUE)
hist(Ho / He)
hohe_cutoff <- 1.5

ggplot(mapping = aes(x = alfreq, y = Ho/He)) +
  geom_point(alpha = 0.05) +
  geom_hline(yintercept = hohe_cutoff, color = "blue") +
  labs(x = "Allele frequency") +
  ggtitle("Filtering for putative paralogs")
mean(Ho/He > hohe_cutoff)
## [1] 0.110541
We’ll do that filtering now, then start a data frame to track statistics about each SNP.

keep1 <- which(Ho/He < hohe_cutoff)

numgen <- numgen[keep1,]
snpdf <- data.frame(row.names = rownames(numgen),
                    AlFreq = alfreq[keep1],
                    Ho = Ho[keep1],
                    He = He[keep1])

snpdf$MissingRate <- rowMeans(is.na(numgen))
summary(snpdf)

hetByInd <- colMeans(numgen == 1, na.rm = TRUE)
hist(hetByInd)

numgen <- numgen[, hetByInd < 0.5]

myclust <- find.clusters(t(numgen), n.pca = 100, n.clust = 10)
myclust$size
mydapc <- dapc(t(numgen), grp = myclust$grp, scale = TRUE, n.pca = 50, n.da = 2)

#mydapc <- dapc(t(numgen), grp = myclust$grp, scale = TRUE, n.pca = 50, n.da = 2)

scatter(mydapc, posi.da = "topright")
colnames(numgen)[mydapc$assign == 9]
colnames(numgen)[mydapc$assign == 4]

grp_alfreq <- sapply(levels(mydapc$assign),
                     function(x) rowMeans(numgen[, mydapc$assign == x], na.rm = TRUE) / 2)
colnames(grp_alfreq) <- paste0("AlFreq_Grp", colnames(grp_alfreq))
head(grp_alfreq)

myvcf <- readVcf(bg,
                 param = ScanVcfParam(geno = NA))

myvcf <- readRDS(rds)

rowRanges(myvcf)
hist(rowRanges(myvcf)$QUAL, xlab = "Quality score",
     main = "Histogram of quality scores in large VCF")

temp <- paste(seqnames(myvcf), start(myvcf), sep = "_")

myvcf <- myvcf[rowRanges(myvcf)$QUAL > 900 | 
                 temp %in% rownames(snpdf),]

snpdf$GCcontent <- gcContent(myvcf, rownames(snpdf), refgenome)

hist(snpdf$GCcontent, xlab = "GC content", main = "Histogram of GC content")
mean(snpdf$GCcontent >= 0.4 & snpdf$GCcontent <= 0.6)

snpdf2 <- filter(snpdf, GCcontent >= 0.4, GCcontent <= 0.6)
snpdf2$Nflanking <- nFlankingSNPs(myvcf, rownames(snpdf2))

hist(snpdf2$Nflanking, xlab = "Number of flanking SNPs",
     main = "Histogram of number of flanking SNPs")

table(snpdf2$Nflanking)
snpdf3 <- filter(snpdf2, Nflanking <= 2)

### choose set of SNP markers 
rr3 <- rowRanges(myvcf)[match(rownames(snpdf3),
                              paste(seqnames(myvcf), start(myvcf), sep = "_"))]
rr3a <- GRanges(sub("chrom", "OM", seqnames(rr3)), ranges(rr3))
subvcf <- readVcf(bg, genome = seqinfo(rr3a),
                  param = ScanVcfParam(which = rr3a))

writeVcf(subvcf, filename = "results/marker_subset.vcf")

markers_purity <- read.delim("results/Galaxy9-[Purity_(beta)_on_data_8].tabular")$Name[1:50]

sort(markers_purity)

markers_simanneal <- findMarkerSet(grp_alfreq[rownames(snpdf3),], nSNP = 3000)$Set
DiffScore(grp_alfreq[markers_simanneal,])
markerseq2 <- formatKasp(myvcf, markers_simanneal, refgenome)
head(markerseq2)
##write the list selected SNP markers 
write.csv(markerseq2, file = "3000_panel_marker_simanneal.csv",
          row.names = FALSE)
###lod marker ID, chrom and position to plot the marker distribution across the genome
selected_marker <- read.csv("position_all_panel_SNP_markers.csv", header = T)
CMplot(selected_marker, plot.type="d", bin.size=1e6, chr.den.col=c("darkgreen", "yellow", "red"),  file="jpg",memo="",dpi=300,
       file.output=TRUE, verbose=TRUE)
####list marker per chromosome
table(selected_marker$Chrom)



































