##################PCA######################################
library(SNPRelate)


########Dataset A (no R filter, min-maf =0)
vcf.fn <- "DatasetA.vcf"

snpgdsVCF2GDS(vcf.fn, "DatasetA.gds",  method="biallelic.only")

genofile <- openfn.gds("DatasetA.gds")

ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

pdf("DatasetA.pdf")  
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2], pch=19, col="yellow")
text(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],labels=ccm_pca$sample.id, cex=0.2)
dev.off()

tq1<-ccm_pca$eigenvect[,1]
tq2<-ccm_pca$eigenvect[,2]
pop<-ccm_pca$sample.id
tq<-cbind(tq1,tq2)
tq<-cbind(pop,tq)
tq<-as.DatasetA.frame(tq)

write.csv(tq,"pca_DatasetA.csv")




###Dataset B (no R filter, min-maf =0.05)
vcf.fn <- "DatasetB.vcf"

snpgdsVCF2GDS(vcf.fn, "DatasetB.gds",  method="biallelic.only")

genofile <- openfn.gds("DatasetB.gds")

ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

pdf("DatasetB.pdf")  
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2], pch=19, col="yellow")
text(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],labels=ccm_pca$sample.id, cex=0.2)
dev.off()

tq1<-ccm_pca$eigenvect[,1]
tq2<-ccm_pca$eigenvect[,2]
pop<-ccm_pca$sample.id
tq<-cbind(tq1,tq2)
tq<-cbind(pop,tq)
tq<-as.DatasetB.frame(tq)

write.csv(tq,"pca_DatasetB.csv")




###Dataset C (R=0.9, no min-maf filter)vcf.fn <- "DatasetC.vcf"

snpgdsVCF2GDS(vcf.fn, "DatasetC.gds",  method="biallelic.only")

genofile <- openfn.gds("DatasetC.gds")

ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

pdf("DatasetC.pdf")  
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2], pch=19, col="yellow")
text(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],labels=ccm_pca$sample.id, cex=0.2)
dev.off()

tq1<-ccm_pca$eigenvect[,1]
tq2<-ccm_pca$eigenvect[,2]
pop<-ccm_pca$sample.id
tq<-cbind(tq1,tq2)
tq<-cbind(pop,tq)
tq<-as.DatasetC.frame(tq)

write.csv(tq,"pca_DatasetC.csv")



###Dataset D (R=0.9, min-maf =0.05)

vcf.fn <- "DatasetD.vcf"

snpgdsVCF2GDS(vcf.fn, "DatasetD.gds",  method="biallelic.only")

genofile <- openfn.gds("DatasetD.gds")

ccm_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

pdf("DatasetD.pdf")  
plot(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2], pch=19, col="yellow")
text(ccm_pca$eigenvect[,1],ccm_pca$eigenvect[,2],labels=ccm_pca$sample.id, cex=0.2)
dev.off()

tq1<-ccm_pca$eigenvect[,1]
tq2<-ccm_pca$eigenvect[,2]
pop<-ccm_pca$sample.id
tq<-cbind(tq1,tq2)
tq<-cbind(pop,tq))
tq<-as.DatasetD.frame(tq)

write.csv(tq,"pca_DatasetD.csv")

##I exported the PCA four Dataset to plot with excel for easier color manipulation


