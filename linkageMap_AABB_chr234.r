setwd("/Users/longhuaguo/Dropbox/0Aging/QTLmapping/QTL_F2_genotyping_2020/novaseq/vcf2")
library(qtl)
library(qtlTools)
library(ASMap)
library(DataCombine)
library(LinkageMapView)
library(tidyverse)

mapthis<-read.cross("csv", file="f2.2parents.aabb.asmpa.rqtl.transformed_noparents.txt", BC.gen=0, F.gen=2, genotypes=c("A","X","B"),na.strings=c("U"))
mapthis

#remove badly genotyped individuals and markers (missing data)
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")
nt.bymar <- ntyped(mapthis, "mar")
todrop <- names(nt.bymar[nt.bymar < 220])
mapthis <- drop.markers(mapthis, todrop)
mapthis <- subsetCross(mapthis, ind=(ntyped(mapthis)>2200))
mapthis
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

#identify and remove clones
cg <- comparegeno(mapthis)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])
wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh
#remove duplicated individuals
g <- pull.geno(mapthis)
table(g[1,], g[4,])
mapthis <- subset(mapthis, ind=-wh[,2])
mapthis

#alternatively, use the function from ASMap#
#identify clones and merge genotypes
gc <- genClones(mapthis, tol = 0.9,id="CHROM_POS_REF_ALT")
gc$cgd
write.table(gc$cgd, file = "f2_2parents_aabb_asmpa_rqtl_transformed_noparents_clones.txt",sep = "\t")
mapthis <- fixClones(mapthis, gc$cgd, id = "CHROM_POS_REF_ALT", consensus = TRUE)
mapthis
plot(ntyped(mapthis), ylab="No. typed markers", main="No. genotypes by individual") 
plot(ntyped(mapthis, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

profileMark(mapthis,stat.type = c("seg.dist","prop","miss"),crit.val = "bonf",layout=c(1,8),type="l",cex=0,display.markers = FALSE)

#remove markers with segregation distortion
sg <- statMark(mapthis, stat.type = "marker")
mm <- statMark(mapthis, stat.type = "marker")$marker$AA
mapthis_seg <- drop.markers(mapthis, c(markernames(mapthis)[mm > 0.35],
                                            markernames(mapthis)[mm < 0.15]))
mapthis_seg
profileMark(mapthis_seg,stat.type = c("seg.dist","prop","miss"),crit.val = "bonf",layout=c(1,8),type="l",cex=0,display.markers = FALSE)

#alternatively, use MSTmap functions for segregation distortion
#####MSTmap variant filtering########
mapBC <- pullCross(mapthis, type = "missing", pars = list(miss.thresh =0.3))
mapBC2 <- pullCross(mapBC, type = "seg.distortion")
mapBC3 <- pullCross(mapBC2, type = "co.located")
sum(ncol(mapBC3$missing$data), ncol(mapBC3$seg.dist$data), ncol(mapBC3$co.located$data))
mapBC
mapBC2
mapBC3

profileMark(mapBC3,stat.type = c("seg.dist","prop","miss"),crit.val = "bonf",layout=c(1,8),type="l",cex=0,display.markers = FALSE)

#examine extent of allele switching
rf <- pull.rf(mapBC3)
lod <- pull.rf(mapBC3, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#construct map
mapBC4 <- mstmap(mapBC3, bychr = F, trace = TRUE, dist.fun =
                   "kosambi", p.value = 1e-10,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)
chrlen(mapBC4)
plotRF(mapBC4,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue",cex=0.2)

####switch alleles and improve the linkages
#heatMap(mapBC4, lmax = 6,color = rev(rainbow(256, start = 0, end = 2/3)))
toswitch <- markernames(mapBC4, chr=c("L.10","L.11","L.17","L.14","L.17","L.7"))
mapBC5 <- switchAlleles(mapBC4, toswitch)
mapBC6 <- mstmap(mapBC5, bychr = F, trace = TRUE, dist.fun =
                   "kosambi", p.value = 1e-10,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)
chrlen(mapBC6)
plotRF(mapBC6,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue")
plotRF(mapBC6,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue",chr=c("L.1","L.10","L.11","L.12","L.13","L.14","L.15","L.16","L.17","L.18","L.19","L2","L.21","L.22","L.3"))

#iterate
toswitch <- markernames(mapBC6, chr=c("L.6","L.10","L.11","L.3","L.15"))
mapBC7 <- switchAlleles(mapBC6, toswitch)
mapBC8 <- mstmap(mapBC7, bychr = F, trace = TRUE, dist.fun =
                   "kosambi", p.value = 1e-10,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)
chrlen(mapBC8)
plotRF(mapBC8,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue")

toswitch <- markernames(mapBC8, chr=c("L.19","l.4","l.5","l.7"))
mapBC9 <- switchAlleles(mapBC8, toswitch)
mapBC10 <- mstmap(mapBC9, bychr = F, trace = TRUE, dist.fun =
                   "kosambi", p.value = 1e-10,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)
chrlen(mapBC10)
plotRF(mapBC10,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue")

toswitch <- markernames(mapBC10, chr=c("L.14"))
mapBC11 <- switchAlleles(mapBC10, toswitch)
mapBC12 <- mstmap(mapBC11, bychr = F, trace = TRUE, dist.fun =
                    "kosambi", p.value = 1e-10,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)
chrlen(mapBC12)
plotRF(mapBC12,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue")

#pull map to analyze scaffolds within clusters
marker_chr <- pull.map(mapBC12,as.table = T)
write.table(marker_chr, file = "f2_qtl_mapBC12_linkagegroups.txt",sep = "\t")

#merge broken linkages into chromosomes
mapDHm <- mergeCross(mapBC12, merge = list(`chr1` = c("L.1", "L.10","L.15","L.9"), `chr2`= c("L.13", "L.14","L.16","L.3","L.4","L.7","L.8"),`chr3`=c("L.12","L.6","L.5"),`chr4`=c("L.11","L.17","L.2")))
chrlen(mapDHm)
plotRF(mapDHm,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue")

#reorder markers
mapBC13 <- mstmap(mapDHm, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 2,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)
chrlen(mapBC13)
plotRF(mapBC13,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue", cex=0.1)
marker_chr <- pull.map(mapBC13,as.table = T)
write.table(marker_chr, file = "f2_qtl_mapBC12_linkagegroups_merged.txt",sep = "\t")
