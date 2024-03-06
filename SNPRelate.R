library(SNPRelate)

##################Estimate kinship#######################################

###Dataset C (R=0.9, no min-maf filter)

vcf.fn <- "DatasetC.vcf"
snpgdsVCF2GDS(vcf.fn, "DatasetC.gds",  method="biallelic.only")
genofile <- openfn.gds("DatasetC.gds")
tqML <- snpgdsIBDMLE(genofile, autosome.only = FALSE, kinship = TRUE)
DatasetC_MLE.ibd <- snpgdsIBDSelection(tqML)
write.csv(DatasetC_MLE.ibd, "DatasetC_MLE.csv")

##NP_R0.9##
NPvcf.fn <- "NP_R0.9_p1.vcf"
snpgdsVCF2GDS(NPvcf.fn, "NP.gds",  method="biallelic.only")
genofile <- openfn.gds("NP.gds")
tqML <- snpgdsIBDMLE(genofile, autosome.only = FALSE, kinship = TRUE)
NP_MLE.ibd <- snpgdsIBDSelection(tqML)
write.csv(NP_MLE.ibd, "NP_MLE.csv")

##NS_R0.9##
NSvcf.fn <- "NS_R0.9_p1.vcf"
sNSgdsVCF2GDS(NSvcf.fn, "NS.gds",  method="biallelic.only")
genofile <- openfn.gds("NS.gds")
tqML <- sNSgdsIBDMLE(genofile, autosome.only = FALSE, kinship = TRUE)
NS_MLE.ibd <- sNSgdsIBDSelection(tqML)
write.csv(NS_MLE.ibd, "NS_MLE.csv")

#MC_R0.9##
MCvcf.fn <- "MC_R0.9_p1.vcf"  
snpgdsVCF2GDS(MCvcf.fn, "MC.gds",  method="biallelic.only")
genofile <- openfn.gds("MC.gds")
tqML <- snpgdsIBDMLE(genofile, autosome.only = FALSE, kinship = TRUE)
MC_MLE.ibd <- snpgdsIBDSelection(tqML)
write.csv(MC_MLE.ibd, "MC_MLE.csv")

