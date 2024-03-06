# Scripts used for the RADSeq data of the Vatica guangxiensis

Manuscript title: The last stand: demographic and population genomic analysis reveals terminal endangerment in tropical timber species Vatica guangxiensis.

Additional information:


## Step 1. Bioinformatic pipeline and variant calling (SNPs.txt)

Programs used: 

FASTP: https://github.com/OpenGene/fastp

STACKS 2.55 (denovo_map.pl): https://catchenlab.life.illinois.edu/stacks/comp/denovo_map.php  

BCFTOOLS: https://www.htslib.org/doc/1.1/bcftools.html 


## Step 2. Genetic divergence with fastStructure (fastStructure.txt) and principal component analysis (PCA.R)

Package used: 

FASTSTRUCTURE: https://rajanil.github.io/fastStructure/

SNPRelate: https://bioconductor.org/packages/release/bioc/html/SNPRelate.html

We also ran an analysis of molecular variance (AMOVA)and computed pairwise FST with ARLEQUIN.
ARLEQUIN: http://cmpg.unibe.ch/software/arlequin35/

##Step 3. R-based population genetic analyses (R-based_analyses.R)

Packages used:

SNPRelate: https://bioconductor.org/packages/release/bioc/html/SNPRelate.html (SNPRelate.R)

HIERFSTAT: https://cran.r-project.org/web/packages/hierfstat/index.html

PopGenome: https://cran.r-project.org/web/packages/PopGenome/index.html



## Step 4. Gene flow 

Program used: 

BAYESASS: https://github.com/brannala/BA3/releases/tag/3.0.5.6

MIGRATE-N: https://popgen.sc.fsu.edu/Migrate/Info.html


## Step 5. Folded site frequency spectrum (SFS) calculation from vcf file

Program used: 

easySFS: https://github.com/isaacovercast/easySFS

## Step 6. Demographical history and phylogenetic reconstruction.

Input files and scripts in the folder "fastsimcoal"

Program used: 

Fastsimcoal2: http://cmpg.unibe.ch/software/fastsimcoal27/

GENEIOUS: https://www.geneious.com/

MAFFT: https://mafft.cbrc.jp/alignment/software/

JMODELTEST: https://evomics.org/learning/phylogenetics/jmodeltest/

BEAST2: https://www.beast2.org/

## Step 7. Future genetic variability

Program used: 

QuantiNemo: https://www2.unil.ch/popgen/softwares/quantinemo/

