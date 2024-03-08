###########################################################Population-based statistics######################################

#######################################################Nucleotide diversity and Tajima's D##############################
##Load the library
library(PopGenome)

#################################Nucleotide diversity and Tajima's D of two genetic populations:GX and YN province
DatasetC <- readData("DIRECTORY/DatasetC", format="VCF", SDatasetC.DATA = TRUE)

##Asign individuals to population (Using T. tridentatus as example)###
pop.GX<-as.character(read.table("GX_popmap.txt")[[1]]) 
pop.YN<-as.character(read.table("YN_popmap.txt")[[1]]) 

DatasetC <- set.populations(DatasetC,list(pop.GX, pop.YN))
DatasetC_2pop <- neutrality.stats(DatasetC)

##Get nucleotide diversity
DatasetC_2pop@nuc.diversity.within

##Get Tajima's D
DatasetC_2pop@Tajima.D



#################################Nucleotide diversity and Tajima's D of three habitat populations:NP,NS,MC
DatasetC <- readData("DIRECTORY/DatasetC", format="VCF", SDatasetC.DATA = TRUE)

##Asign individuals to population (Using T. tridentatus as example)###
pop.NP<-as.character(read.table("NP_popmap.txt")[[1]]) 
pop.NS<-as.character(read.table("NS_popmap.txt")[[1]]) 
pop.MC<-as.character(read.table("MC_popmap.txt")[[1]])
DatasetC <- set.populations(DatasetC,list(pop.NP, pop.NS, pop.MC))
DatasetC_3pop <- neutrality.stats(DatasetC)

##Get nucleotide diversity
DatasetC_3pop@nuc.diversity.within

##Get Tajima's D
DatasetC_3pop@Tajima.D



##################################Nucleotide diversity and Tajima's D of NP
NP <- readData("DIRECTORY/NP_R0.9_p1", format="VCF", SNP.DATA = TRUE)

##Asign individuals to population (Using T. tridentatus as example)###ad means adult, sa means sapling, se means seedlings
pop.NP_ad <- c("Adults_samplename")
pop.NP_sa <- c("Saplings_samplename")
pop.NP_se <- c("Seedlings_samplename")

NP <- set.populations(NP,list(pop.NP_ad, pop.NP_sa, pop.NP_se))
NP <- neutrality.stats(NP)

##Get nucleotide diversity
NP@nuc.diversity.within

##Get Tajima's D
NP@Tajima.D


####################################Nucleotide diversity and Tajima's D of NS
NS <- readData("DIRECTORY/NS_R0.9_p1", format="VCF", SNS.DATA = TRUE)

##Asign individuals to population (Using T. tridentatus as example)###ad means adult, sa means sapling, se means seedlings
pop.NS_ad <- c("Adults_samplename")
pop.NS_se <- c("Seedlings_samplename")
NS <- set.populations(NSH,list(pop.NS_ad, pop.NS_se))
NS <- set.populations(NS,list(pop.NS_ad, pop.NS_sa, pop.NS_se))
NS <- neutrality.stats(NS)

##Get nucleotide diversity
NS@nuc.diversity.within

##Get Tajima's D
NS@Tajima.D


######################################Nucleotide diversity and Tajima's D of MC
MC <- readData("DIRECTORY/MC_R0.9_p1", format="VCF", SNP.DATA = TRUE)
##Asign individuals to population (Using T. tridentatus as example)###mad means adult, msa means sapling, mse means seedlings
pop.mad <- c("Adults_samplename")
pop.msa <- c("Saplings_samplenamen")
pop.mse <- c("Seedlings_samplename")

MC <- set.populations(MC,list(pop.mad, pop.msa, pop.mse))
MC <- neutrality.stats(MC)

##Get nucleotide diversity
MC@nuc.diversity.within

##Get Tajima's D
MC@Tajima.D


####################################################################Calculate Ho, Fis, Ar####################
library(hierfstat)

##########################################Two genetic populations: GX and YN province
####Dataset C (R=0.9, no min-maf filter)
DatasetC <- import2genind("DatasetC.gen")
DatasetC_population <- read.csv("two_province_popmap.csv", header=T)
DatasetC_dat <- as.factor(DatasetC_population$pop)
DatasetC_genind@pop <- DatasetC_dat
DatasetC_genind@pop
DatasetC_basic_stat <- basic.stats(DatasetC_genind, diploid=TRUE)
write.csv(basic_stat$perloc,"DatasetC_two_pop_basicstats.csv")

####Ho
Ho <- basic_stat$Ho
mean(Ho[,1], na.rm=TRUE)
mean(Ho[,2], na.rm=TRUE)

####Fis
Fis <- basic_stat$Fis
mean(Fis[,1], na.rm=TRUE)
mean(Fis[,2], na.rm=TRUE)

####Ar
DatasetC_Ar <- DatasetCelic.richness(DatasetC_genind,min.n=NULL,diploid=TRUE)
DatasetCelic.richness <- DatasetC_Ar$Ar
mean(DatasetCelic.richness[,1], na.rm=TRUE)
mean(DatasetCelic.richness[,2], na.rm=TRUE)




##########################################three habitat populations: NP, NS and MC
####Dataset C (R=0.9, no min-maf filter)
DatasetC <- import2genind("DatasetC.gen")
DatasetC_population <- read.csv("three_popmap.csv", header=T)
DatasetC_dat <- as.factor(DatasetC_population$pop)
DatasetC_genind@pop <- DatasetC_dat
DatasetC_genind@pop
DatasetC_basic_stat <- basic.stats(DatasetC_genind, diploid=TRUE)
write.csv(basic_stat$perloc,"DatasetC_three_pop_basicstats.csv")

####Ho
Ho <- basic_stat$Ho
mean(Ho[,1], na.rm=TRUE)
mean(Ho[,2], na.rm=TRUE)
mean(Ho[,3])

####Fis
Fis <- basic_stat$Fis
mean(Fis[,1], na.rm=TRUE)
mean(Fis[,2], na.rm=TRUE)
mean(Fis[,3], na.rm=TRUE)

####Ar
DatasetC_Ar <- DatasetCelic.richness(DatasetC_genind,min.n=NULL,diploid=TRUE)
DatasetCelic.richness <- DatasetC_Ar$Ar
mean(DatasetCelic.richness[,1], na.rm=TRUE)
mean(DatasetCelic.richness[,2], na.rm=TRUE)
mean(DatasetCelic.richness[,3], na.rm=TRUE)


##########################################NP population
NP <- import2genind("NP.gen")
NP_population <- read.csv("NP_popmap.csv", header=T)
NP_dat <- as.factor(NP_population$pop)
NP_genind@pop <- NP_dat
NP_genind@pop
NP_basic_stat <- basic.stats(NP_genind, diploid=TRUE)
write.csv(basic_stat$perloc,"NP_pop_basicstats.csv")

####Ho
Ho <- basic_stat$Ho
mean(Ho[,1], na.rm=TRUE)
mean(Ho[,2], na.rm=TRUE)
mean(Ho[,3])

####Fis
Fis <- basic_stat$Fis
mean(Fis[,1], na.rm=TRUE)
mean(Fis[,2], na.rm=TRUE)
mean(Fis[,3], na.rm=TRUE)

####Ar
NP_Ar <- NPelic.richness(NP_genind,min.n=NULL,diploid=TRUE)
NPelic.richness <- NP_Ar$Ar
mean(NPelic.richness[,1], na.rm=TRUE)
mean(NPelic.richness[,2], na.rm=TRUE)
mean(NPelic.richness[,3], na.rm=TRUE)


##########################################NS population
NS <- import2genind("NS.gen")
NS_population <- read.csv("NS_popmap.csv", header=T)
NS_dat <- as.factor(NS_population$pop)
NS_genind@pop <- NS_dat
NS_genind@pop
NS_basic_stat <- basic.stats(NS_genind, diploid=TRUE)
write.csv(basic_stat$perloc,"NS_pop_basicstats.csv")

####Ho
Ho <- basic_stat$Ho
mean(Ho[,1], na.rm=TRUE)
mean(Ho[,2], na.rm=TRUE)
mean(Ho[,3])

####Fis
Fis <- basic_stat$Fis
mean(Fis[,1], na.rm=TRUE)
mean(Fis[,2], na.rm=TRUE)
mean(Fis[,3], na.rm=TRUE)

####Ar
NS_Ar <- NSelic.richness(NS_genind,min.n=NULL,diploid=TRUE)
NSelic.richness <- NS_Ar$Ar
mean(NSelic.richness[,1], na.rm=TRUE)
mean(NSelic.richness[,2], na.rm=TRUE)
mean(NSelic.richness[,3], na.rm=TRUE)

##########################################MC population
MC <- import2genind("MC.gen")
MC_population <- read.csv("MC_popmap.csv", header=T)
MC_dat <- as.factor(MC_population$pop)
MC_genind@pop <- MC_dat
MC_genind@pop
MC_basic_stat <- basic.stats(MC_genind, diploid=TRUE)
write.csv(basic_stat$perloc,"MC_pop_basicstats.csv")

####Ho
Ho <- basic_stat$Ho
mean(Ho[,1], na.rm=TRUE)
mean(Ho[,2], na.rm=TRUE)
mean(Ho[,3])

####Fis
Fis <- basic_stat$Fis
mean(Fis[,1], na.rm=TRUE)
mean(Fis[,2], na.rm=TRUE)
mean(Fis[,3], na.rm=TRUE)

####Ar
MC_Ar <- MCelic.richness(MC_genind,min.n=NULL,diploid=TRUE)
MCelic.richness <- MC_Ar$Ar
mean(MCelic.richness[,1], na.rm=TRUE)
mean(MCelic.richness[,2], na.rm=TRUE)
mean(MCelic.richness[,3], na.rm=TRUE)

