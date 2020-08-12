setwd("~/Dropbox/7spermseq_linkagemap/s2_feb22")
library(qtl)
library(qtlTools)
library(ASMap)
library(DataCombine)
library(LinkageMapView)
library(tidyverse)

#convert the MSTmap input file to CSV file for read.cross. 
gts <- read.delim("s2_rad_sperm_rqtlInput", sep = "\t", header = TRUE)
as.matrix(gts)
gts2 <- t(gts)
write.csv(gts2, file = "s2_rad_sperm_rqtlInput.transformed")

mapthis <- read.cross("csv",file="s2_rad_sperm_rqtlInput.transformed", sep=",",genotypes = c("0","1"))
class(mapthis)[1] <- "dh" 
mapthis

#missing data
#plotMissing(mapthis)
plot(ntyped(mapthis), ylab= "No. typed markers", main="Typed genotypes by individual") 
nt.bymar <- ntyped(mapthis, "mar")
plot(nt.bymar, ylab= "No. typed individuals", main="No. individuals by marker") 
todrop <- names(nt.bymar[nt.bymar < 45])
todrop
sg <- statGen(mapthis, bychr = FALSE, stat.type = "miss",id="CHROM_POS_REF_ALT")
mapthis_miss <- subset(mapthis, ind = sg$miss < 50000)
mapthis_miss <- drop.markers(mapthis_miss, todrop)
mapthis_miss
#genotyping clones; remove clones if needed
gc <- genClones(mapthis_miss, tol = 0.95,id="locus_name")
#segregation distortion
profileMark(mapthis,stat.type = c("seg.dist","prop","miss"),crit.val = "bonf",layout=c(1,4),type="l",cex=0)
profileMark(mapthis_miss,stat.type = c("seg.dist","prop","miss"),crit.val = "bonf",layout=c(1,3),type="l",cex=0)
sg <- statMark(mapthis_miss, stat.type = "marker")
mm <- statMark(mapthis_miss, stat.type = "marker")$marker$AA
mapthis_seg <- drop.markers(mapthis_miss, c(markernames(mapthis_miss)[mm > 0.7],
                                            markernames(mapthis_miss)[mm < 0.3]))
mapthis_seg
profileMark(mapthis_seg,stat.type = c("seg.dist","prop","miss"),crit.val = "bonf",layout=c(1,4),type="l",cex=0.4)
rf <- pull.rf(mapthis_seg)
lod <- pull.rf(mapthis_seg, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")


#####MSTmap variant filtering########
mapBC <- pullCross(mapthis_miss, type = "missing", pars = list(miss.thresh =0.3))
mapBC2 <- pullCross(mapBC, type = "seg.distortion", pars =list(seg.ratio = "80:20"))
mapBC3 <- pullCross(mapBC2, type = "co.located")
sum(ncol(mapBC3$missing$data), ncol(mapBC3$seg.dist$data),
    ncol(mapBC3$co.located$data))
mapBC
mapBC2
mapBC3
profileMark(mapBC3,stat.type = c("seg.dist","prop","miss"),crit.val = "bonf",layout=c(1,4),type="l",cex=0.4)
rf <- pull.rf(mapBC3)
lod <- pull.rf(mapBC3, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

#construct map
mapBC4 <- mstmap(mapBC3, bychr = F, trace = TRUE, dist.fun =
                   "kosambi", p.value = 1e-5,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)
chrlen(mapBC4)
plotRF(mapBC4,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue")

####switch alleles and improve the linkages
heatMap(mapBC4, lmax = 6,color = rev(rainbow(256, start = 0, end = 2/3)))
toswitch <- markernames(mapBC4, chr=c("L.14","L.27","L.3","L.30","L.32","L.33","L.35","L.40","L.43","L.48","L.49","L.57","L.59","L.6","L.9"))
mapBC5 <- switchAlleles(mapBC4, toswitch)
mapBC6 <- mstmap(mapBC5, bychr = F, trace = TRUE, dist.fun =
                   "kosambi", p.value = 1e-5,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)

toswitch <- markernames(mapBC6, chr=c("L.9","L.19","L.26","L.27","L.33","L.34","L.35","L.40","L.38","L.39"))
mapBC7 <- switchAlleles(mapBC6, toswitch)
mapBC8 <- mstmap(mapBC7, bychr = F, trace = TRUE, dist.fun =
                   "kosambi", p.value = 1e-5,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)
mapBC8 <- mstmap(mapBC8, bychr = T, trace = TRUE, dist.fun =
                   "kosambi", p.value = 2,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)

#iterate
heatMap(mapBC6, lmax = 6,color = rev(rainbow(256, start = 0, end = 2/3)))
plotRF(mapBC6,alternate.chrid = T,zmax = 6,col.scheme = "redblue")
plotRF(mapBC8,alternate.chrid = T,zmax = 6,col.scheme = "redblue")

chrlen(mapBC6)

rf <- pull.rf(mapBC8)
lod <- pull.rf(mapBC8, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

marker_chr <- pull.map(mapBC8,as.table = T)
write.table(marker_chr, file = "s2_linkage_rad_august_48sperms.txt",sep = "\t")

#merge broken linkages into chromosomes
mapDHm <- mergeCross(mapBC8, merge = list(`chr1` = c("L.2", "L.21","L.24","L.33","L.34","L.36","L.4","L.38"), `chr2`= c("L.1", "L.15","L.12","L.13","L.14","L.16","L.17","L.18","L.20","L.25","L.26","L.28","L.27","L.5","L.6","L.7","L.8"),`chr3`=c("L.11","L.19","L.22","L.23","L.31","L.32","L.30","L.37","L.9"),`chr4`=c("L.10","L.35")))

chrlen(mapDHm)
plotRF(mapDHm,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue")

#reorder markers
mapBC14 <- mstmap(mapDHm, bychr = T, trace = TRUE, dist.fun = "kosambi", p.value = 2,id="CHROM_POS_REF_ALT", objective.fun = "COUNT", mvest.bc = T, detectBadData = T, return.imputed = T)
chrlen(mapBC14)
plotRF(mapBC14,alternate.chrid = T,zmax = 5.5,col.scheme = "redblue", cex=0.1,chr = c("chr1","chr2","chr3","chr4"))

marker_chr <- pull.map(mapBC14,as.table = T)
write.table(marker_chr, file = "s2_linkage_rad_august_48sperms_merged.txt",sep = "\t")

saveRDS(mapBC14, "s2_linkage_rad_august_48sperms_merged.rds")
