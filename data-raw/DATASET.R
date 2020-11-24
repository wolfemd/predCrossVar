## code to prepare `DATASET` dataset goes here
library(tidyverse); library(magrittr); library(AlphaSimR)
founderPop = quickHaplo(nInd=100, nChr=2, segSites=50)
SP<-SimParam$new(founderPop)
# set-up 2 genetically correlated traits each with Additive and Dominance efects
G = 1.5*diag(2)-0.5 #Genetic correlation matrix
SP$addTraitAD(nQtlPerChr = 10, mean=c(0,0), meanDD=c(0,0), var=c(1,1), varDD = c(1,1), corA=G, corDD=G)

parent_pop <- newPop(founderPop) # initial generation
parent_pop %<>% setPheno(simParam = SP)

# store the pop-object
usethis::use_data(parent_pop, overwrite = TRUE)

# centimorgan scale map as a named vector
## will be converted to recomb. freq. matrix in tutorial
genmap<-getQtlMap()
m<-genmap$pos;
names(m)<-paste0(genmap$chr,"_",genmap$site)
genmap<-m; rm(m)

# store the genetic map
usethis::use_data(genmap, overwrite = TRUE)

haploMat<-pullQtlHaplo(pop)
# for now, it seems the critical internal functions `predOneCrossVarA()` and `predOneCrossVarAD()`
# require haplotypes for the same individual to be distinguished by "_HapA" and "_HapB"
# will fix in a future version
rownames(haploMat) %<>%
      gsub("_1","_HapA",.) %>%
      gsub("_2","_HapB",.)
# column IDs for the loci should match in order and ID between all files
colnames(haploMat)<-names(genmap)

# store the haplotype matrix
usethis::use_data(haploMat, overwrite = TRUE)

# effects
AddEffectsList<-list(Trait1=matrix(SP$traits[[1]]@addEff,nrow=1) %>% set_colnames(.,colnames(haploMat)),
                     Trait2=matrix(SP$traits[[2]]@addEff,nrow=1) %>% set_colnames(.,colnames(haploMat)))
DomEffectsList<-list(Trait1=matrix(SP$traits[[1]]@domEff,nrow=1) %>% set_colnames(.,colnames(haploMat)),
                     Trait2=matrix(SP$traits[[2]]@domEff,nrow=1) %>% set_colnames(.,colnames(haploMat)))

usethis::use_data(AddEffectsList, overwrite = TRUE)
usethis::use_data(DomEffectsList, overwrite = TRUE)

#save(pop,genmap,haploMat,AddEffectsList,DomEffectsList,file=here::here("data","example_data.Rdata"))

