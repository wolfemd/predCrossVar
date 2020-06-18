
#' effectsArray2list
#'
#' Converts an array of posterior samples of multi-trait marker effects to a named list (one for each trait).
#'
#' @param effectsArray According to BGLR documentation: 3D array, with dim=c(nRow,p,traits), where nRow number of MCMC samples saved, p is the number of predictors and traits is the number of traits. BGLR::Multitrait() writes a binary file to disk when saveEffects=TRUE is specified. It can be read to R with BGLR::readBinMatMultitrait().
#' @param snpIDs character vector with labels for the predictors (SNPs), numeric should work too, but untested.
#' @param traits character vector to label the traits.
#' @param nIter number of iterations used for MCMC; used internally only to exclude burn-in samples from computation
#' @param burnIn burnIn for MCMC; used internally only to exclude burn-in samples from computation
#' @param thin thin for MCMC; used internally only to exclude burn-in samples from computation
#'
#' @return list of matrices, one matrix per trait, each matrix has `nrow((nIter-burnIn)/thin)` and `ncol(length(snpIDs))`. Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @export
#'
#' @examples
#' effectsArray<-BGLR::readBinMatMultitrait(filename = "mt_effects_ETA_g_beta.bin")
#' effectsArray<-effectsArray2list(effectsArray, snpIDs, traits, nIter, burnIn, thin)
#'
effectsArray2list<-function(effectsArray, snpIDs, traits, nIter, burnIn, thin){
   # Discard burnIn
   effectsArray<-effectsArray[-c(1:burnIn/thin),,]
   # Add dimnames
   dimnames(effectsArray)[[2]]<-snpIDs
   dimnames(effectsArray)[[3]]<-traits
   # 3D arrays of effects to lists-of-matrices (by trait)
   effectsArray<-purrr::array_branch(effectsArray,3)
   return(effectsArray)
}

#' posteriorMeanVarCovarA
#'
#' For an additive only model, compute multiple types of genomic estimator for a single variance or covariance component.
#' Posterior mean variance or co-variance (PMV) and also the variance of posterior means (VPM) estimators are returned.
#' Both M1, i.e. "Method 1", which does not account for LD in the population and M2 or "Method 2", which does, are returned.
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2  string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#' @param genoVarCovarMat variance-covariance matrix of marker genotypes, i.e. measuring linkage disequilibrium
#'
#' @return tibble with estimators for one variance or co-variance parameter.
#' @export
#'
#' @examples
posteriorMeanVarCovarA<-function(Trait1,Trait2,
                                 AddEffectList,postMeanAddEffects,genoVarCovarMat){

   # Posterior Sample Variance-Covariance Matrix of Marker Effects
   postVarCovarOfAddEffects<-(1/(nrow(AddEffectList[[Trait1]])-1))*crossprod(AddEffectList[[Trait1]],AddEffectList[[Trait2]])

   # Method 1 (Unconditional, Not accounting for LD)
   ## Additive
   #### (Co)Variance of Posterior Means
   vpm_m1a<-sum(diag(genoVarCovarMat)*(postMeanAddEffects[[Trait1]]*postMeanAddEffects[[Trait2]]))
   #### Posterior Mean (Co)Variance
   pmv_m1a<-vpm_m1a+sum(diag(postVarCovarOfAddEffects)*diag(genoVarCovarMat))

   # Method 2 (Conditioned on LD)
   ## Additive
   #### (Co)Variance of Posterior Means
   vpm_m2a<-postMeanAddEffects[[Trait1]]%*%genoVarCovarMat%*%postMeanAddEffects[[Trait2]]
   #### Posterior Mean (Co)Variance
   pmv_m2a<-vpm_m2a+sum(diag(genoVarCovarMat%*%postVarCovarOfAddEffects))

   # Tidy the results
   out<-bind_rows(tibble(VarComp=c("VarA"),
                         Method=c("M1"),
                         VPM=c(vpm_m1a), # NOTE: unsure relevance of "VPM" for Method 1, suggest ignoring
                         PMV=c(pmv_m1a)), #### PMV for Method 1 matches the standard VarComps you would get from BGLR
                  tibble(VarComp=c("VarA"),
                         Method=c("M2"),
                         VPM=c(vpm_m2a),
                         PMV=c(pmv_m2a)))
   return(out)
}



#' posteriorMeanVarCovarAD
#'
#' For an additive plus dominance model, compute multiple types of genomic estimator for a single variance or covariance component.
#' Posterior mean variance or co-variance (PMV) and also the variance of posterior means (VPM) estimators are returned.
#' Both M1, i.e. "Method 1", which does not account for LD in the population and M2 or "Method 2", which does, are returned.
#' Additive and dominance variances/covariances are returned.
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2 string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param DomEffectList list of DOMINANCE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#' @param postMeanDomEffects list of named vectors (or column matrices) with the posterior mean DOMINANCE marker effects.
#' @param genoVarCovarMat variance-covariance matrix of marker genotypes, i.e. measuring linkage disequilibrium
#'
#' @return tibble with estimators for one variance or co-variance parameter including both additive and dominance components.
#' @export
#'
#' @examples
posteriorMeanVarCovarAD<-function(Trait1,Trait2,
                                  AddEffectList,DomEffectList,
                                  postMeanAddEffects,postMeanDomEffects,
                                  genoVarCovarMat){

   # Posterior Sample Variance-Covariance Matrix of Marker Effects
   postVarCovarOfAddEffects<-(1/(nrow(AddEffectList[[Trait1]])-1))*crossprod(AddEffectList[[Trait1]],AddEffectList[[Trait2]])
   postVarCovarOfDomEffects<-(1/(nrow(DomEffectList[[Trait1]])-1))*crossprod(DomEffectList[[Trait1]],DomEffectList[[Trait2]])

   # Method 1 (Unconditional, Not accounting for LD)
   ## Additive
   #### (Co)Variance of Posterior Means
   vpm_m1a<-sum(diag(genoVarCovarMat)*(postMeanAddEffects[[Trait1]]*postMeanAddEffects[[Trait2]]))
   #### Posterior Mean (Co)Variance
   pmv_m1a<-vpm_m1a+sum(diag(postVarCovarOfAddEffects)*diag(genoVarCovarMat))
   ## Dominance
   #### (Co)Variance of Posterior Means
   vpm_m1d<-sum(diag(genoVarCovarMat)^2*(postMeanDomEffects[[Trait1]]*postMeanDomEffects[[Trait2]]))
   #### Posterior Mean (Co)Variance
   pmv_m1d<-vpm_m1d+sum(diag(postVarCovarOfDomEffects)*(diag(genoVarCovarMat)^2))

   # Method 2 (Conditioned on LD)
   ## Additive
   #### (Co)Variance of Posterior Means
   vpm_m2a<-postMeanAddEffects[[Trait1]]%*%genoVarCovarMat%*%postMeanAddEffects[[Trait2]]
   #### Posterior Mean (Co)Variance
   pmv_m2a<-vpm_m2a+sum(diag(genoVarCovarMat%*%postVarCovarOfAddEffects))
   ## Dominance
   #### (Co)Variance of Posterior Means
   genoVarCovarMatSq<-genoVarCovarMat^2
   vpm_m2d<-postMeanDomEffects[[Trait1]]%*%genoVarCovarMatSq%*%postMeanDomEffects[[Trait2]]
   #### Posterior Mean (Co)Variance
   pmv_m2d<-vpm_m2d+sum(diag(genoVarCovarMatSq%*%postVarCovarOfDomEffects))

   # Tidy the results
   out<-bind_rows(tibble(VarComp=c("VarA","VarD"),
                         Method=c("M1","M1"),
                         VPM=c(vpm_m1a,vpm_m1d), # NOTE: unsure relevance of "VPM" for Method 1, suggest ignoring
                         PMV=c(pmv_m1a,pmv_m1d)), #### PMV for Method 1 matches the standard VarComps you would get from BGLR
                  tibble(VarComp=c("VarA","VarD"),
                         Method=c("M2","M2"),
                         VPM=c(vpm_m2a,vpm_m2d),
                         PMV=c(pmv_m2a,pmv_m2d)))
   return(out)
}

#' getMultiTraitPMVs_A
#'
#' For an additive only model, compute _all pairwise_ variance and co-variance parameters from the post-burn-In, thinned samples of marker effects obtained using an Bayesian multi-trait model fit to some current population.
#' Wraps `posteriorMeanVarCovarA()` over all trait variances and co-variances.
#' Posterior mean variance or co-variance (PMV) and also the variance of posterior means (VPM) estimators are returned.
#' Both M1, i.e. "Method 1", which does not account for LD in the population and M2 or "Method 2", which does, are returned.
#'
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param genoVarCovarMat variance-covariance matrix of marker genotypes, i.e. measuring linkage disequilibrium
#'
#' @return tibble with estimators for all pairwise variance and co-variance parameters.
#' @export
#'
#' @examples
getMultiTraitPMVs_A<-function(AddEffectList, genoVarCovarMat){
   traits<-names(AddEffectList)
   # Center posterior distribution of effects
   ## on posterior mean across MCMC samples
   AddEffectList %<>% map(.,~scale(.,center = T, scale = F))

   ## Get the posterior mean effects vectors
   postMeanAddEffects<-map(AddEffectList,~attr(.,which = "scaled:center"))

   # Make tibble of pairwise variance parameters to compute
   ## Include trait-trait variances, avoid duplicating covariances
   ## i.e. lower triangle including diagonals of var-covar matrix
   varcovars<-bind_rows(tibble(Trait1=traits,Trait2=traits), # trait variances
                        combn(traits,2,simplify = T) %>% # covariances
                           t(.) %>% #
                           `colnames<-`(.,c("Trait1","Trait2")) %>%
                           as_tibble)

   # Compute over each variance parameter
   varcovars %<>%
      mutate(varcomps=pmap(.,posteriorMeanVarCovarA,AddEffectList,postMeanAddEffects)) %>%
      unnest(varcomps)
   return(varcovars)
}


#' getMultiTraitPMVs_AD
#'
#' For an additive plus dominance model, compute _all pairwise_ variance and co-variance parameters from the post-burn-In, thinned samples of marker effects obtained using an Bayesian multi-trait model fit to some current population.
#' Wraps `posteriorMeanVarCovarAD()` over all trait variances and co-variances.
#' Posterior mean variance or co-variance (PMV) and also the variance of posterior means (VPM) estimators are returned.
#' Both M1, i.e. "Method 1", which does not account for LD in the population and M2 or "Method 2", which does, are returned.
#' Additive and dominance variances/covariances are returned.
#'
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param DomEffectList list of DOMINANCE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param genoVarCovarMat variance-covariance matrix of marker genotypes, i.e. measuring linkage disequilibrium
#'
#' @return tibble with estimators for all pairwise variance and co-variance parameters including both additive and dominance components.
#' @export
#'
#' @examples
getMultiTraitPMVs_AD<-function(AddEffectList, DomEffectList, genoVarCovarMat){
   traits<-names(AddEffectList)
   # Center posterior distribution of effects
   ## on posterior mean across MCMC samples
   AddEffectList %<>% map(.,~scale(.,center = T, scale = F))
   DomEffectList %<>% map(.,~scale(.,center = T, scale = F))

   ## Get the posterior mean effects vectors
   postMeanAddEffects<-map(AddEffectList,~attr(.,which = "scaled:center"))
   postMeanDomEffects<-map(DomEffectList,~attr(.,which = "scaled:center"))

   # Make tibble of pairwise variance parameters to compute
   ## Include trait-trait variances, avoid duplicating covariances
   ## i.e. lower triangle including diagonals of var-covar matrix
   varcovars<-bind_rows(tibble(Trait1=traits,Trait2=traits), # trait variances
                        combn(traits,2,simplify = T) %>% # covariances
                           t(.) %>% #
                           `colnames<-`(.,c("Trait1","Trait2")) %>%
                           as_tibble)

   # Compute over each variance parameter
   varcovars %<>%
      mutate(varcomps=pmap(.,posteriorMeanVarCovarAD,AddEffectList,DomEffectList,postMeanAddEffects,postMeanDomEffects)) %>%
      unnest(varcomps)
   return(varcovars)
}

#' predCrossMeanBVsOneTrait
#'
#' Predict the mean breeding value of each family, for a single trait, given parental allelic dosages and (posterior mean) marker effects.
#' Should from an additive-only model.
#'
#' @param Trait string, label of trait (name in list of postMeanAddEffects) to compute
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param doseMat dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/colnames to indicate SNP/ind ID
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#'
#' @return tibble with parental GEBV and the pred Mean GEBV (mean of parents) for each cross.
#' @export
#'
#' @examples
predCrossMeanBVsOneTrait<-function(Trait,CrossesToPredict,doseMat,postMeanAddEffects){
   parentGEBVs<-doseMat%*%postMeanAddEffects[[Trait]]
   predictedCrossMeanBVs<-CrossesToPredict %>%
      left_join(tibble(sireID=rownames(parentGEBVs),sireGEBV=as.numeric(parentGEBVs))) %>%
      left_join(tibble(damID=rownames(parentGEBVs),damGEBV=as.numeric(parentGEBVs))) %>%
      mutate(predMeanBV=(sireGEBV+damGEBV)/2)
   return(predictedCrossMeanBVs)
}

#' predCrossMeanGVsOneTrait
#'
#' Predict the total genetic merit of the cross based on a model with additive+dominance marker effects.
#' For each family , for a single trait, given  parental allelic dosages and (posterior mean) marker effects.
#' G = sum( pr(AA)*a-pr(aa)*a+pr(Aa)*d ), where A is counted allele in dosages
#'
#' @param Trait string, label of trait (name in list of postMeanAddEffects) to compute
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param doseMat dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/colnames to indicate SNP/ind ID
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#' @param postMeanDomEffects list of named vectors (or column matrices) with the posterior mean DOMINANCE marker effects.
#'
#' @return tibble with predicted mean GV for each family
#' @export
#'
#' @examples
predCrossMeanGVsOneTrait<-function(Trait,CrossesToPredict,doseMat,postMeanAddEffects,postMeanDomEffects){
   predictedCrossMeanGVs<-CrossesToPredict %>%
      mutate(predMeanGV=map2_dbl(sireID,damID,
                                 function(sireID,damID){
                                    p1<-doseMat[sireID,]/2
                                    p2<-doseMat[damID,]/2
                                    q1<-1-p1
                                    q2<-1-p2
                                    gfreqs<-cbind(q1*q2,p1*q2+p2*q1,p1*p2)

                                    meanG<-sum((gfreqs[,3]-gfreqs[,1])*postMeanAddEffects[[Trait]]+gfreqs[,2]*postMeanDomEffects[[Trait]])
                                    return(meanG)
                                 }))
   return(predictedCrossMeanGVs) }

#' predCrossMeansA
#'
#' From an additive only model, fit to multiple traits, predict the mean (breeding value) of each cross for each trait.
#' Corresponds to the mean GEBV of the parents, given  parental allelic dosages and (posterior mean) marker effects.
#'
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#' @param doseMat dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/colnames to indicate SNP/ind ID
#'
#' @return tibble with predicted mean BV for each trait in each family
#' @export
#'
#' @examples
predCrossMeansA<-function(CrossesToPredict,postMeanAddEffects,doseMat){
   means<-tibble(Trait=names(postMeanAddEffects))
   parents<-CrossesToPredict %$% union(sireID,damID)
   doseMat<-doseMat[parents,names(postMeanAddEffects[[1]])]
   ## Predicted Mean Total Merit
   means %<>%
      mutate(predMeanBVs=map(Trait,predCrossMeanBVsOneTrait,CrossesToPredict,doseMat,postMeanAddEffects)) %>%
      select(Trait,predMeanBVs) %>%
      unnest(predMeanBVs)
   return(means) }

#' predCrossMeansAD
#'
#' From an additive+dominance model, fit to multiple traits, predict the total genetic merit of each cross for each trait.
#' For each family , for a single trait, given  parental allelic dosages and (posterior mean) marker effects.
#' G = sum( pr(AA)*a-pr(aa)*a+pr(Aa)*d ), where A is counted allele in dosages
#'
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#' @param postMeanDomEffects list of named vectors (or column matrices) with the posterior mean DOMINANCE marker effects.
#' @param doseMat dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/colnames to indicate SNP/ind ID
#'
#' @return tibble with predicted mean GV for each trait in each family
#' @export
#'
#' @examples
predCrossMeansAD<-function(CrossesToPredict,postMeanAddEffects,postMeanDomEffects,doseMat){
   means<-tibble(Trait=names(postMeanAddEffects))
   parents<-CrossesToPredict %$% union(sireID,damID)
   doseMat<-doseMat[parents,names(postMeanAddEffects[[1]])]
   ## Predicted Mean Total Merit
   means %<>%
      mutate(predMeanGVs=map(Trait,predCrossMeanGVsOneTrait,CrossesToPredict,doseMat,postMeanAddEffects,postMeanDomEffects)) %>%
      select(Trait,predMeanGVs) %>%
      unnest(predMeanGVs)
   return(means) }

#' predOneCrossVarA
#'
#' Predict the additive variance or co-variance in one cross.
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2 string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param sireID string, Sire genotype ID. Needs to correspond to renames in haploMat
#' @param damID string, Dam genotype ID. Needs to correspond to renames in haploMat
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param predType string, "VPM" or "PMV". Default = "VPM" for variance of posterior means, this is faster but expected to be less accurate / more biased than the alternative predType=="PMV". PMV requires user to supply a variance-covariance matrix of effects estimates.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#' @param postVarCovarOfAddEffects matrix of dimension N SNP x N SNP. ADDITIVE Posterior Sample Variance-Covariance Matrix of Marker Effects Estimates.
#' @param ...
#'
#' @return tibble with predicted additive variance for one cross, one variance parameter
#' @export
#'
#' @examples
predOneCrossVarA<-function(Trait1,Trait2,sireID,damID,
                           haploMat,recombFreqMat,predType="VPM",
                           postMeanAddEffects,
                           postVarCovarOfAddEffects=NULL,...){
   starttime<-proc.time()[3]

   # Before predicting variances
   # check for and remove SNPs that
   # won't segregate, i.e. are fixed in parents
   ### hopes to save time / mem
   x<-colSums(rbind(haploMat[grep(paste0(sireID,"_"),rownames(haploMat)),],
                    haploMat[grep(paste0(damID,"_"),rownames(haploMat)),]))
   segsnps2keep<-names(x[x>0 & x<4])

   if(length(segsnps2keep)>0){
      recombFreqMat<-recombFreqMat[segsnps2keep,segsnps2keep]
      haploMat<-haploMat[sort(c(paste0(parents,"_HapA"),paste0(parents,"_HapB"))),segsnps2keep]
      postMeanAddEffects<-map(postMeanAddEffects,~.[segsnps2keep])
      if(predType=="PMV"){
         postVarCovarOfAddEffects<-postVarCovarOfAddEffects[segsnps2keep,segsnps2keep]
      }
      # sire and dam LD matrices
      sireLD<-calcGameticLD(sireID,recombFreqMat,haploMat)
      damLD<-calcGameticLD(damID,recombFreqMat,haploMat)
      progenyLD<-sireLD+damLD

      rm(recombFreqMat,haploMat,sireLD,damLD); gc()

      ## Additive
      #### (Co)Variance of Posterior Means
      vpm_m2a<-postMeanAddEffects[[Trait1]]%*%progenyLD%*%postMeanAddEffects[[Trait2]]
      #### Posterior Mean (Co)Variance
      if(predType=="PMV"){ pmv_m2a<-vpm_m2a+sum(diag(progenyLD%*%postVarCovarOfAddEffects)) }
      totcomputetime<-proc.time()[3]-starttime

      rm(progenyLD); gc()

      ## Tidy the results
      out<-tibble(VarComp=c("VarA"),
                  VPM=c(vpm_m2a),
                  PMV=ifelse(predType=="PMV",c(pmv_m2a),c(NA)),
                  Nsegsnps=c(length(segsnps2keep)),
                  totcomputetime=c(totcomputetime))
      print(paste0("Variances predicted for family: ",sireID,"x",damID,"- took ",round(totcomputetime/60,3)," mins"))
   } else {
      out<-tibble(VarComp=c("VarA"),
                  VPM=c(0),
                  PMV=c(0),
                  Nsegsnps=c(0),
                  computetime=c(0))
      print(paste0("Variances predicted for family: ",sireID,"x",damID,"- had no segregating SNPs"))
   }
   return(out)
}

#' predOneCrossVarAD
#'
#' Predict the additive and dominance variance or co-variance in one cross.
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2 string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param sireID string, Sire genotype ID. Needs to correspond to renames in haploMat
#' @param damID string, Dam genotype ID. Needs to correspond to renames in haploMat
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param predType string, "VPM" or "PMV". Default = "VPM" for variance of posterior means, this is faster but expected to be less accurate / more biased than the alternative predType=="PMV". PMV requires user to supply a variance-covariance matrix of effects estimates.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#' @param postMeanDomEffects list of named vectors (or column matrices) with the posterior mean DOMINANCE marker effects.
#' @param postVarCovarOfAddEffects matrix of dimension N SNP x N SNP. ADDITIVE Posterior Sample Variance-Covariance Matrix of Marker Effects Estimates.
#' @param postVarCovarOfDomEffects matrix of dimension N SNP x N SNP. DOMINANCE Posterior Sample Variance-Covariance Matrix of Marker Effects Estimates.
#' @param ...
#'
#' @return tibble with predicted additive and dominance variance for one cross, one variance parameter
#' @export
#'
#' @examples
predOneCrossVarAD<-function(Trait1,Trait2,sireID,damID,
                            haploMat,recombFreqMat,predType="VPM",
                            postMeanAddEffects,postMeanDomEffects,
                            postVarCovarOfAddEffects=NULL,postVarCovarOfDomEffects=NULL,...){
   starttime<-proc.time()[3]

   # Before predicting variances
   # check for and remove SNPs that
   # won't segregate, i.e. are fixed in parents
   ### hopes to save time / mem
   x<-colSums(rbind(haploMat[grep(paste0(sireID,"_"),rownames(haploMat)),],
                    haploMat[grep(paste0(damID,"_"),rownames(haploMat)),]))
   segsnps2keep<-names(x[x>0 & x<4])

   if(length(segsnps2keep)>0){
      recombFreqMat<-recombFreqMat[segsnps2keep,segsnps2keep]
      haploMat<-haploMat[sort(c(paste0(parents,"_HapA"),paste0(parents,"_HapB"))),segsnps2keep]
      postMeanAddEffects<-map(postMeanAddEffects,~.[segsnps2keep])
      postMeanDomEffects<-map(postMeanDomEffects,~.[segsnps2keep])
      if(predType=="PMV"){
         postVarCovarOfAddEffects<-postVarCovarOfAddEffects[segsnps2keep,segsnps2keep]
         postVarCovarOfDomEffects<-postVarCovarOfDomEffects[segsnps2keep,segsnps2keep]
      }
      # sire and dam LD matrices
      sireLD<-calcGameticLD(sireID,recombFreqMat,haploMat)
      damLD<-calcGameticLD(damID,recombFreqMat,haploMat)
      progenyLD<-sireLD+damLD

      rm(recombFreqMat,haploMat,sireLD,damLD); gc()

      ## Additive
      #### (Co)Variance of Posterior Means
      vpm_m2a<-postMeanAddEffects[[Trait1]]%*%progenyLD%*%postMeanAddEffects[[Trait2]]
      #### Posterior Mean (Co)Variance
      if(predType=="PMV"){ pmv_m2a<-vpm_m2a+sum(diag(progenyLD%*%postVarCovarOfAddEffects)) }
      ## Dominance
      #### (Co)Variance of Posterior Means
      progenyLDsq<-progenyLD^2
      vpm_m2d<-postMeanDomEffects[[Trait1]]%*%progenyLDsq%*%postMeanDomEffects[[Trait2]]
      #### Posterior Mean (Co)Variance
      if(predType=="PMV"){ pmv_m2d<-vpm_m2d+sum(diag(progenyLDsq%*%postVarCovarOfDomEffects)) }
      totcomputetime<-proc.time()[3]-starttime

      rm(progenyLDsq,progenyLD); gc()

      ## Tidy the results
      out<-tibble(VarComp=c("VarA","VarD"),
                  VPM=c(vpm_m2a,vpm_m2d),
                  PMV=ifelse(predType=="PMV",c(pmv_m2a,pmv_m2d),c(NA,NA)),
                  Nsegsnps=c(length(segsnps2keep),NA),
                  totcomputetime=c(totcomputetime,NA))
      print(paste0("Variances predicted for family: ",sireID,"x",damID,"- took ",round(totcomputetime/60,3)," mins"))
   } else {
      out<-tibble(VarComp=c("VarA","VarD"),
                  VPM=c(0,0),
                  PMV=c(0,0),
                  Nsegsnps=c(0,0),
                  computetime=c(0,0))
      print(paste0("Variances predicted for family: ",sireID,"x",damID,"- had no segregating SNPs"))
   }
   return(out)
}

#' predCrossVarsA
#'
#' Predict the additive variance or co-variance for a set of crosses, potentially in parallel across families.
#' Wraps predOneCrossVarA across families.
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2 string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param predType string, "VPM" or "PMV". Default = "VPM" for variance of posterior means, this is faster but expected to be less accurate / more biased than the alternative predType=="PMV". PMV requires user to supply a variance-covariance matrix of effects estimates.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param ncores If ncores set > 1 parallelizes across families, but beware it is memory intensive and options(future.globals.maxSize=___) may need to be adjusted.
#' @param ...
#'
#' @return list with two elements, "predictedfamvars" contains a tibble with all predictions for all requested families, "totcomputetime" gives the time taken to compute one var. parameter across all familes, at the given ncores.
#' @export
#'
#' @examples
predCrossVarsA<-function(Trait1,Trait2,CrossesToPredict,predType="VPM",
                         haploMat,recombFreqMat,
                         postMeanAddEffects,
                         AddEffectList=NULL,
                         ncores=1,...){

   starttime<-proc.time()[3]
   if(predType=="PMV"){
      # Posterior Sample Variance-Covariance Matrix of Marker Effects
      postVarCovarOfAddEffects<-(1/(nrow(AddEffectList[[Trait1]])-1))*crossprod(AddEffectList[[Trait1]],AddEffectList[[Trait2]])
      rm(AddEffectList); gc()
   } else {
      postVarCovarOfAddEffects<-NULL;
   }

   require(furrr); options(mc.cores=ncores); plan(multiprocess)
   predictedfamvars<-CrossesToPredict %>%
      dplyr::mutate(predVars=future_pmap(.,
                                         predOneCrossVarA,
                                         Trait1=Trait1,Trait2=Trait2,
                                         haploMat=haploMat,recombFreqMat=recombFreqMat,
                                         predType=predType,
                                         postMeanAddEffects=postMeanAddEffects,
                                         postVarCovarOfAddEffects=postVarCovarOfAddEffects))
   totcomputetime<-proc.time()[3]-starttime
   print(paste0("Done predicting fam vars. ",
                "Took ",round((totcomputetime)/60,2),
                " mins for ",nrow(predictedfamvars)," families"))
   predictedfamvars<-list(predictedfamvars=predictedfamvars,totcomputetime=totcomputetime)
   return(predictedfamvars)
}

#' predCrossVarsAD
#'
#' Predict the additive and dominance variance or co-variance for a set of crosses, potentially in parallel across families.
#' Wraps predOneCrossVarAD across families
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2 string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param predType string, "VPM" or "PMV". Default = "VPM" for variance of posterior means, this is faster but expected to be less accurate / more biased than the alternative predType=="PMV". PMV requires user to supply a variance-covariance matrix of effects estimates.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.
#' @param postMeanDomEffects list of named vectors (or column matrices) with the posterior mean DOMINANCE marker effects.
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param DomEffectList list of DOMINANCE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param ncores If ncores set > 1 parallelizes across families, but beware it is memory intensive and options(future.globals.maxSize=___) may need to be adjusted.
#' @param ...
#'
#' @return list with two elements, "predictedfamvars" contains a tibble with all predictions for all requested families, "totcomputetime" gives the time taken to compute one var. parameter across all familes, at the given ncores.
#' @export
#'
#' @examples
predCrossVarsAD<-function(Trait1,Trait2,CrossesToPredict,predType="VPM",
                          haploMat,recombFreqMat,
                          postMeanAddEffects,postMeanDomEffects,
                          AddEffectList=NULL,DomEffectList=NULL,
                          ncores=1,...){

   starttime<-proc.time()[3]
   if(predType=="PMV"){
      # Posterior Sample Variance-Covariance Matrix of Marker Effects
      postVarCovarOfAddEffects<-(1/(nrow(AddEffectList[[Trait1]])-1))*crossprod(AddEffectList[[Trait1]],AddEffectList[[Trait2]])
      postVarCovarOfDomEffects<-(1/(nrow(DomEffectList[[Trait1]])-1))*crossprod(DomEffectList[[Trait1]],DomEffectList[[Trait2]])
      rm(AddEffectList,DomEffectList); gc()
   } else {
      postVarCovarOfAddEffects<-NULL;
      postVarCovarOfDomEffects<-NULL;
   }

   require(furrr); options(mc.cores=ncores); plan(multiprocess)
   predictedfamvars<-CrossesToPredict %>%
      dplyr::mutate(predVars=future_pmap(.,
                                         predOneCrossVarAD,
                                         Trait1=Trait1,Trait2=Trait2,
                                         haploMat=haploMat,recombFreqMat=recombFreqMat,
                                         predType=predType,
                                         postMeanAddEffects=postMeanAddEffects,
                                         postMeanDomEffects=postMeanDomEffects,
                                         postVarCovarOfAddEffects=postVarCovarOfAddEffects,
                                         postVarCovarOfDomEffects=postVarCovarOfDomEffects))
   totcomputetime<-proc.time()[3]-starttime
   print(paste0("Done predicting fam vars. ",
                "Took ",round((totcomputetime)/60,2),
                " mins for ",nrow(predictedfamvars)," families"))
   predictedfamvars<-list(predictedfamvars=predictedfamvars,totcomputetime=totcomputetime)
   return(predictedfamvars)
}


#' runMtCrossVarPredsA
#'
#' Predict the additive variances _and_ covariances for a set of crosses, potentially in parallel across families.
#' Wraps predCrossVarsAD across variance parameters.
#'
#' @param outprefix string, prefix for *.rds file to be written to disk with output. DEFAULT = NULL (no disk write)
#' @param outpath string, path to disk location where files should be written. Can be left DEFAULT = NULL (no disk write)
#' @param predType string, "VPM" or "PMV". Default = "VPM" for variance of posterior means, this is faster but expected to be less accurate / more biased than the alternative predType=="PMV". PMV requires user to supply a variance-covariance matrix of effects estimates.
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param ncores If ncores set > 1 parallelizes across families, but beware it is memory intensive and options(future.globals.maxSize=___) may need to be adjusted.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runMtCrossVarPredsA<-function(outprefix=NULL,outpath=NULL,predType="VPM",
                              CrossesToPredict,AddEffectList,
                              haploMat,recombFreqMat,ncores=1,...){
   starttime<-proc.time()[3]
   traits<-names(AddEffectList)
   parents<-CrossesToPredict %$% union(sireID,damID)

   # Center posterior distribution of effects
   ## on posterior mean across MCMC samples
   AddEffectList %<>% map(.,~scale(.,center = T, scale = F))

   ## Get the posterior mean effects vectors
   postMeanAddEffects<-map(AddEffectList,~attr(.,which = "scaled:center"))

   ## Predict trait (co)variances
   varcovars<-bind_rows(tibble(Trait1=traits,Trait2=traits), # trait variances
                        combn(traits,2,simplify = T) %>% # covariances
                           t(.) %>% #
                           `colnames<-`(.,c("Trait1","Trait2")) %>%
                           as_tibble)

   haploMat<-haploMat[sort(c(paste0(parents,"_HapA"),paste0(parents,"_HapB"))),]

   if(predType!="PMV"){
      AddEffectList<-NULL;
   }

   varcovars %<>%
      mutate(varcomps=pmap(.,predCrossVarsA,CrossesToPredict=CrossesToPredict,predType=predType,
                           AddEffectList=AddEffectList,
                           haploMat=haploMat,recombFreqMat=recombFreqMat,
                           postMeanAddEffects=postMeanAddEffects,ncores=ncores))

   totcomputetime<-proc.time()[3]-starttime
   varcovars<-list(varcovars=varcovars,
                   totcomputetime=totcomputetime)
   if(!is.null(outprefix) & !is.null(outpath)){
      saveRDS(varcovars,file=here::here(outpath,paste0(outprefix,"_predVarsAndCovars.rds")))
   }
   return(varcovars)
}

#' runMtCrossVarPredsAD
#'
#' Predict the additive and dominance variances _and_ covariances for a set of crosses, potentially in parallel across families.
#' Wraps predCrossVarsAD across variance parameters.
#'
#' @param outprefix string, prefix for *.rds file to be written to disk with output. DEFAULT = NULL (no disk write)
#' @param outpath string, path to disk location where files should be written. Can be left DEFAULT = NULL (no disk write)
#' @param predType string, "VPM" or "PMV". Default = "VPM" for variance of posterior means, this is faster but expected to be less accurate / more biased than the alternative predType=="PMV". PMV requires user to supply a variance-covariance matrix of effects estimates.
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param AddEffectList list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param DomEffectList list of DOMINANCE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param ncores If ncores set > 1 parallelizes across families, but beware it is memory intensive and options(future.globals.maxSize=___) may need to be adjusted.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runMtCrossVarPredsAD<-function(outprefix=NULL,outpath=NULL,predType="VPM",
                               CrossesToPredict,AddEffectList,DomEffectList,
                               haploMat,recombFreqMat,ncores=1,...){
   starttime<-proc.time()[3]
   traits<-names(AddEffectList)
   parents<-CrossesToPredict %$% union(sireID,damID)

   # Center posterior distribution of effects
   ## on posterior mean across MCMC samples
   AddEffectList %<>% map(.,~scale(.,center = T, scale = F))
   DomEffectList %<>% map(.,~scale(.,center = T, scale = F))

   ## Get the posterior mean effects vectors
   postMeanAddEffects<-map(AddEffectList,~attr(.,which = "scaled:center"))
   postMeanDomEffects<-map(DomEffectList,~attr(.,which = "scaled:center"))

   ## Predict trait (co)variances
   varcovars<-bind_rows(tibble(Trait1=traits,Trait2=traits), # trait variances
                        combn(traits,2,simplify = T) %>% # covariances
                           t(.) %>% #
                           `colnames<-`(.,c("Trait1","Trait2")) %>%
                           as_tibble)

   haploMat<-haploMat[sort(c(paste0(parents,"_HapA"),paste0(parents,"_HapB"))),]

   if(predType!="PMV"){
      AddEffectList<-NULL;
      DomEffectList<-NULL;
   }

   varcovars %<>%
      mutate(varcomps=pmap(.,predCrossVarsAD,CrossesToPredict=CrossesToPredict,predType=predType,
                           AddEffectList=AddEffectList,DomEffectList=DomEffectList,
                           haploMat=haploMat,recombFreqMat=recombFreqMat,
                           postMeanAddEffects=postMeanAddEffects,postMeanDomEffects=postMeanDomEffects,ncores=ncores))

   totcomputetime<-proc.time()[3]-starttime
   varcovars<-list(varcovars=varcovars,
                   totcomputetime=totcomputetime)
   if(!is.null(outprefix) & !is.null(outpath)){
      saveRDS(varcovars,file=here::here(outpath,paste0(outprefix,"_predVarsAndCovars.rds")))
   }
   return(varcovars)
}



