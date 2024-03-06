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
pop.NP_ad <- c("VGNP34.clean")
pop.NP_sa <- c("VGNP35.clean","VGNP105.clean")
pop.NP_se <- c("VGNP100.clean","VGNP102.clean","VGNP103.clean","VGNP104.clean","VGNP106.clean","VGNP107.clean","VGNP108.clean","VGNP10_all.clean","VGNP110.clean","VGNP13.clean","VGNP14.clean","VGNP16.clean","VGNP17.clean","VGNP1.clean","VGNP20.clean","VGNP21_all.clean","VGNP22.clean","VGNP23.clean","VGNP24.clean","VGNP25.clean","VGNP26_all.clean","VGNP27.clean","VGNP28.clean","VGNP2.clean","VGNP31.clean","VGNP32.clean","VGNP36.clean","VGNP37.clean","VGNP38.clean","VGNP39.clean","VGNP3.clean","VGNP40.clean","VGNP41.clean","VGNP42.clean","VGNP44.clean","VGNP46.clean","VGNP48.clean","VGNP4.clean","VGNP50.clean","VGNP56.clean","VGNP59.clean","VGNP5.clean","VGNP63.clean","VGNP64.clean","VGNP65.clean","VGNP66.clean","VGNP69.clean","VGNP6.clean","VGNP70.clean","VGNP72.clean","VGNP73.clean","VGNP74.clean","VGNP75.clean","VGNP76.clean","VGNP77.clean","VGNP78.clean","VGNP79.clean","VGNP7.clean","VGNP80.clean","VGNP82.clean","VGNP83.clean","VGNP84.clean","VGNP85.clean","VGNP86.clean","VGNP87.clean","VGNP88.clean","VGNP89.clean","VGNP90.clean","VGNP91.clean","VGNP92.clean","VGNP93.clean","VGNP94.clean","VGNP96.clean","VGNP97.clean","VGNP98.clean","VGNP99.clean","VGNP9.clean")

NP <- set.populations(NP,list(pop.NP_ad, pop.NP_sa, pop.NP_se))
NP <- neutrality.stats(NP)

##Get nucleotide diversity
NP@nuc.diversity.within

##Get Tajima's D
NP@Tajima.D


####################################Nucleotide diversity and Tajima's D of NS
NS <- readData("DIRECTORY/NS_R0.9_p1", format="VCF", SNS.DATA = TRUE)

##Asign individuals to population (Using T. tridentatus as example)###ad means adult, sa means sapling, se means seedlings
pop.NS_ad <- c("VXRN135.clean","VXRN160.clean","VXRN161.clean","VXRN33.clean")
pop.NS_se <- c("VXRN02.clean","VXRN03.clean","VXRN04.clean","VXRN100.clean","VXRN104.clean","VXRN106.clean","VXRN107.clean","VXRN113.clean","VXRN114.clean","VXRN116.clean","VXRN117.clean","VXRN119.clean","VXRN120.clean","VXRN121.clean","VXRN122.clean","VXRN126.clean","VXRN128.clean","VXRN12.clean","VXRN130.clean","VXRN131.clean","VXRN134.clean","VXRN136.clean","VXRN137.clean","VXRN138.clean","VXRN139.clean","VXRN141.clean","VXRN142.clean","VXRN143.clean","VXRN144.clean","VXRN145.clean","VXRN146.clean","VXRN147.clean","VXRN148.clean","VXRN149.clean","VXRN150.clean","VXRN152.clean","VXRN153.clean","VXRN154.clean","VXRN155.clean","VXRN159.clean","VXRN162.clean","VXRN163.clean","VXRN16.clean","VXRN24.clean","VXRN25.clean","VXRN28.clean","VXRN30.clean","VXRN34.clean","VXRN52.clean","VXRN57.clean","VXRN65.clean")
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
pop.mad <- c("VXRM178.clean","VXRM179.clean","VXRM183.clean","VXRM199.clean","VXRM219.clean","VXRM230.clean","VXRM231.clean","VXRM232.clean","VXRM233.clean","VXRM235.clean","VXRM236.clean","VXRM237.clean","VXRM238.clean","VXRM239.clean","VXRM241.clean","VXRM242.clean","VXRM243.clean","VXRM244.clean","VXRM254.clean","VXRM256.clean","VXRM257.clean","VXRM259.clean")
pop.msa <- c("VXRM165.clean","VXRM167.clean","VXRM168.clean","VXRM169.clean","VXRM170.clean","VXRM171.clean","VXRM175.clean","VXRM180.clean","VXRM181.clean","VXRM185.clean","VXRM186.clean","VXRM193.clean","VXRM194.clean","VXRM198.clean","VXRM208.clean","VXRM210.clean","VXRM211.clean","VXRM212.clean","VXRM213.clean","VXRM214.clean","VXRM222.clean","VXRM229.clean","VXRM234.clean","VXRM240.clean","VXRM253.clean","VXRM258.clean","VXRM262.clean","VXRM32.clean")
pop.mse <- c("VXRM166.clean","VXRM172.clean","VXRM173.clean","VXRM174.clean","VXRM176.clean","VXRM177.clean","VXRM182.clean","VXRM184.clean","VXRM188.clean","VXRM189.clean","VXRM190.clean","VXRM191.clean","VXRM192.clean","VXRM195.clean","VXRM197.clean","VXRM200.clean","VXRM201.clean","VXRM202.clean","VXRM203.clean","VXRM204.clean","VXRM205.clean","VXRM206.clean","VXRM209.clean","VXRM215.clean","VXRM216.clean","VXRM218.clean","VXRM220.clean","VXRM223.clean","VXRM224.clean","VXRM225.clean","VXRM227.clean","VXRM245.clean","VXRM247.clean","VXRM248.clean","VXRM249.clean","VXRM251.clean","VXRM252.clean","VXRM255.clean","VXRM264.clean","VXRM266.clean","VXRM267.clean","VXRM57.clean","VXRM62.clean","VXRM97.clean")

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

