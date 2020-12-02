
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
   out<-dplyr::bind_rows(tibble::tibble(VarComp=c("VarA"),
                                        Method=c("M1"),
                                        VPM=c(vpm_m1a), # NOTE: unsure relevance of "VPM" for Method 1, suggest ignoring
                                        PMV=c(pmv_m1a)), #### PMV for Method 1 matches the standard VarComps you would get from BGLR
                         tibble::tibble(VarComp=c("VarA"),
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
   out<-dplyr::bind_rows(tibble::tibble(VarComp=c("VarA","VarD"),
                                        Method=c("M1","M1"),
                                        VPM=c(vpm_m1a,vpm_m1d), # NOTE: unsure relevance of "VPM" for Method 1, suggest ignoring
                                        PMV=c(pmv_m1a,pmv_m1d)), #### PMV for Method 1 matches the standard VarComps you would get from BGLR
                         tibble::tibble(VarComp=c("VarA","VarD"),
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
   AddEffectList<-purrr::map(AddEffectList,~scale(.,center = T, scale = F))

   ## Get the posterior mean effects vectors
   postMeanAddEffects<-purrr::map(AddEffectList,~attr(.,which = "scaled:center"))

   # Make tibble of pairwise variance parameters to compute
   ## Include trait-trait variances, avoid duplicating covariances
   ## i.e. lower triangle including diagonals of var-covar matrix
   varcovars<-dplyr::bind_rows(tibble::tibble(Trait1=traits,Trait2=traits), # trait variances
                               combn(traits,2,simplify = T) %>% # covariances
                                  t(.) %>% #
                                  `colnames<-`(.,c("Trait1","Trait2")) %>%
                                  tibble::as_tibble)

   # Compute over each variance parameter
   varcovars<-varcovars %>%
      dplyr::mutate(varcomps=purrr::pmap(.,posteriorMeanVarCovarA,
                           AddEffectList=AddEffectList,
                           postMeanAddEffects=postMeanAddEffects,
                           genoVarCovarMat=genoVarCovarMat)) %>%
      tidyr::unnest(varcomps)
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
   AddEffectList<-purrr::map(AddEffectList,~scale(.,center = T, scale = F))
   DomEffectList<-purrr::map(DomEffectList,~scale(.,center = T, scale = F))

   ## Get the posterior mean effects vectors
   postMeanAddEffects<-purrr::map(AddEffectList,~attr(.,which = "scaled:center"))
   postMeanDomEffects<-purrr::map(DomEffectList,~attr(.,which = "scaled:center"))

   # Make tibble of pairwise variance parameters to compute
   ## Include trait-trait variances, avoid duplicating covariances
   ## i.e. lower triangle including diagonals of var-covar matrix
   varcovars<-dplyr::bind_rows(tibble(Trait1=traits,Trait2=traits), # trait variances
                               combn(traits,2,simplify = T) %>% # covariances
                                  t(.) %>% #
                                  `colnames<-`(.,c("Trait1","Trait2")) %>%
                                  tibble::as_tibble)

   # Compute over each variance parameter
   varcovars<-varcovars %>%
      dplyr::mutate(varcomps=pmap(.,posteriorMeanVarCovarAD,
                                  AddEffectList=AddEffectList,
                                  DomEffectList=DomEffectList,
                                  postMeanAddEffects=postMeanAddEffects,
                                  postMeanDomEffects=postMeanDomEffects,
                                  genoVarCovarMat=genoVarCovarMat)) %>%
      tidyr::unnest(varcomps)
   return(varcovars)
}

#' predCrossMeanBVsOneTrait
#'
#' Predict the mean breeding value of each family, for a single trait, given parental allelic dosages
#' and (posterior mean) marker effects.
#' NOTE: Marker effects should represent allele substitution effects.
#' Either by fitting an additive-only model OR an genotypic-partitioned additive+dominance model,
#' with allele sub effects computed as a+d(q-p), where q and p are allele freqs in the training pop. used.
#'
#' @param Trait string, label of trait (name in list of postMeanAddEffects) to compute
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param doseMat dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/colnames to indicate SNP/ind ID
#' @param postMeanAlleleSubEffects list of named vectors (or column matrices) with the posterior mean ALLELE SUBSTITUTION marker effects.
#'
#' @return tibble with parental GEBV and the pred Mean GEBV (mean of parents) for each cross.
#' @export
#'
#' @examples
predCrossMeanBVsOneTrait<-function(Trait,CrossesToPredict,doseMat,postMeanAlleleSubEffects){
   parentGEBVs<-doseMat%*%postMeanAlleleSubEffects[[Trait]]
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
#' Eqn. 14.6 Falconer+MacKay, and elsewhere
#' G = sum( 𝑎(𝑝 − 𝑞 − 𝑦) + 𝑑[2𝑝𝑞 + 𝑦(𝑝 − 𝑞)] )
#' a and d being the additive and dominance effects
#' p and q being the allele frequencies of one parent
#' y is the difference of freq. between the two parents
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
predCrossMeanGVsOneTrait<-function(Trait,CrossesToPredict,doseMat,
                                   postMeanAddEffects,postMeanDomEffects){
   predictedCrossMeanGVs<-CrossesToPredict %>%
      mutate(predMeanGV=map2_dbl(sireID,damID,
                                 function(sireID,damID){
                                    # Eqn 14.6 from Falconer+MacKay
                                    p1<-doseMat[sireID,]/2
                                    p2<-doseMat[damID,]/2
                                    q<-1-p1
                                    y<-p1-p2
                                    g<-postMeanAddEffects[[Trait]]*(p1-q-y) + postMeanDomEffects[[Trait]]*((2*p1*q)+y*(p1-q))
                                    meanG<-sum(g)
                                    # Equivalent, e.g. Toro & Varona 2010
                                    # q1<-1-p1
                                    # q2<-1-p2
                                    # gfreqs<-cbind(q1*q2,p1*q2+p2*q1,p1*p2)
                                    # meanG<-sum((gfreqs[,3]-gfreqs[,1])*postMeanAddEffects[[Trait]]+gfreqs[,2]*postMeanDomEffects[[Trait]])
                                    return(meanG)
                                 }))
   return(predictedCrossMeanGVs) }

#' predCrossMeanBVs
#'
#' From an additive only model, fit to multiple traits, predict the mean (breeding value) of each cross for each trait.
#' Corresponds to the mean GEBV of the parents, given  parental allelic dosages and (posterior mean) marker effects.
#'
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param postMeanAlleleSubEffects list of named vectors (or column matrices) with the posterior mean ALLELE SUBSTITUTION marker effects.
#' @param doseMat dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/colnames to indicate SNP/ind ID
#'
#' @return tibble with predicted mean BV for each trait in each family
#' @export
#'
#' @examples
predCrossMeanBVs<-function(CrossesToPredict,postMeanAlleleSubEffects,doseMat){
   means<-tibble(Trait=names(postMeanAlleleSubEffects))
   parents<-CrossesToPredict %$% union(sireID,damID)
   doseMat<-doseMat[parents,names(postMeanAlleleSubEffects[[1]])]
   ## Predicted Mean Total Merit
   means<-means %>%
      mutate(predMeanBVs=map(Trait,predCrossMeanBVsOneTrait,CrossesToPredict,doseMat,postMeanAlleleSubEffects)) %>%
      select(Trait,predMeanBVs) %>%
      unnest(predMeanBVs)
   return(means) }

#' predCrossMeanTGVs
#'
#' From an additive+dominance model, fit to multiple traits, predict the total genetic merit of each cross for each trait.
#' For each family , for a single trait, given  parental allelic dosages and (posterior mean) marker effects.
#' G = sum( 𝑎(𝑝 − 𝑞 − 𝑦) + 𝑑[2𝑝𝑞 + 𝑦(𝑝 − 𝑞)] )
#' a and d being the additive and dominance effects
#' p and q being the allele frequencies of one parent
#' y is the difference of freq. between the two parents
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
predCrossMeanTGVs<-function(CrossesToPredict,postMeanAddEffects,postMeanDomEffects,doseMat){
   means<-tibble(Trait=names(postMeanAddEffects))
   parents<-CrossesToPredict %$% union(sireID,damID)
   doseMat<-doseMat[parents,names(postMeanAddEffects[[1]])]
   ## Predicted Mean Total Merit
   means<-means %>%
      mutate(predMeanGVs=map(Trait,predCrossMeanGVsOneTrait,CrossesToPredict,doseMat,postMeanAddEffects,postMeanDomEffects)) %>%
      select(Trait,predMeanGVs) %>%
      unnest(predMeanGVs)
   return(means) }

#' Multi-trait prediction of one additive genetic variance (or covariance) among full-siblings.
#'
#' User specifies a trait variance (or trait-trait covariance) to predict for a specific pair of parents. Predicts the addtiive genetic variance (or covariance) among full-siblings of that cross.
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2 string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param sireID string, Sire genotype ID. Needs to correspond to renames in haploMat
#' @param damID string, Dam genotype ID. Needs to correspond to renames in haploMat
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param predType string, "VPM" or "PMV". Choose option "VPM" if you have REML marker effect estimates (or posterior-means from MCMC) one set of marker effect estimates per trait. Variance of posterior means is faster but the alternative predType=="PMV" is expected to be less biassed. PMV requires user to supply a (probably LARGE) variance-covariance matrix of effects estimates.
#' @param postMeanAlleleSubEffects list of named vectors (or column matrices) with the allele substitution marker effects (can posterior-mean effects from MCMC _or_ from REML, if if setting predType="PMV".
#' @param postVarCovarOfAlleleSubEffects Only if setting predType="PMV". Matrix of dimension N SNP x N SNP. ALLELE SUBSTITUTION Posterior Sample Variance-Covariance Matrix of Marker Effects Estimates.
#' @param ...
#'
#' @return tibble with predicted additive variance for one cross, one variance parameter
#' @export
#'
#' @examples
predOneCrossVarA<-function(Trait1,Trait2,sireID,damID,
                           haploMat,recombFreqMat,predType,
                           postMeanAlleleSubEffects,
                           postVarCovarOfAlleleSubEffects=NULL,...){
   starttime<-proc.time()[3]

   # Before predicting variances
   # check for and remove SNPs that
   # won't segregate, i.e. are fixed in parents
   ### hopes to save time / mem
   x<-colSums(rbind(haploMat[grep(paste0("^",sireID,"_"),rownames(haploMat)),],
                    haploMat[grep(paste0("^",damID,"_"),rownames(haploMat)),]))
   segsnps2keep<-names(x[x>0 & x<4])

   if(length(segsnps2keep)>0){
      recombFreqMat<-recombFreqMat[segsnps2keep,segsnps2keep,drop=F]
      haploMat<-haploMat[,segsnps2keep,drop=F]

      postMeanAlleleSubEffects<-purrr::map(postMeanAlleleSubEffects,~.[segsnps2keep])
      if(predType=="PMV"){
         postVarCovarOfAlleleSubEffects<-postVarCovarOfAlleleSubEffects[segsnps2keep,segsnps2keep,drop=F]
      }
      # sire and dam LD matrices
      sireLD<-calcGameticLD(sireID,recombFreqMat,haploMat)
      damLD<-calcGameticLD(damID,recombFreqMat,haploMat)
      progenyLD<-sireLD+damLD

      rm(recombFreqMat,haploMat,sireLD,damLD); gc()

      ## Additive
      #### (Co)Variance of Posterior Means
      vpm_m2a<-postMeanAlleleSubEffects[[Trait1]]%*%progenyLD%*%postMeanAlleleSubEffects[[Trait2]]
      #### Posterior Mean (Co)Variance
      if(predType=="PMV"){ pmv_m2a<-vpm_m2a+sum(diag(progenyLD%*%postVarCovarOfAlleleSubEffects)) }
      totcomputetime<-proc.time()[3]-starttime

      rm(progenyLD); gc()

      ## Tidy the results
      # out<-tibble::tibble(VarComp=c("VarA"),
      #                     VPM=c(vpm_m2a),
      #                     PMV=ifelse(predType=="PMV",c(pmv_m2a),c(NA)),
      #                     Nsegsnps=c(length(segsnps2keep)),
      #                     totcomputetime=c(totcomputetime))
      out<-tibble::tibble(VarComp="VarA",VPM=c(vpm_m2a),PMV=ifelse(predType=="PMV",pmv_m2a,NA)) %>%
         dplyr::mutate(Nsegsnps=length(segsnps2keep),
                       totcomputetime=c(totcomputetime))



      print(paste0("Variances predicted for family: ",sireID,"x",damID,"- took ",round(totcomputetime/60,3)," mins"))
   } else {
      out<-tibble::tibble(VarComp=c("VarA"),
                          VPM=c(0),
                          PMV=c(0),
                          Nsegsnps=c(0),
                          computetime=c(NA))
      print(paste0("Variances predicted for family: ",sireID,"x",damID,"- had no segregating SNPs"))
   }
   return(out)
}

#' Multi-trait prediction of one additive _and_ one dominance genetic variance (or covariance) among full-siblings.
#'
#' User specifies a trait variance (or trait-trait covariance) to predict for a specific pair of parents. Predicts the additive genetic variance (or covariance) among full-siblings of that cross.
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2 string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param sireID string, Sire genotype ID. Needs to correspond to renames in haploMat
#' @param damID string, Dam genotype ID. Needs to correspond to renames in haploMat
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param predType string, "VPM" or "PMV". Choose option "VPM" if you have REML marker effect estimates (or posterior-means from MCMC) one set of marker effect estimates per trait. Variance of posterior means is faster but the alternative predType=="PMV" is expected to be less biassed. PMV requires user to supply a (probably LARGE) variance-covariance matrix of effects estimates.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the additive marker effects (can posterior-mean effects from MCMC _or_ from REML, if setting predType="PMV".
#' @param postMeanDomEffects list of named vectors (or column matrices) with the dominance marker effects (can posterior-mean effects from MCMC _or_ from REML, if setting predType="PMV".
#' @param postVarCovarOfAddEffects Only if setting predType="PMV". Matrix of dimension N SNP x N SNP. ADDITIVE Posterior Sample Variance-Covariance Matrix of Marker Effects Estimates.
#' @param postVarCovarOfDomEffects Only if setting predType="PMV". Matrix of dimension N SNP x N SNP. DOMINANCE Posterior Sample Variance-Covariance Matrix of Marker Effects Estimates.
#' @param ...
#'
#' @return tibble with predicted additive and dominance variance for one cross, one variance parameter
#' @export
#'
#' @examples
predOneCrossVarAD<-function(Trait1,Trait2,sireID,damID,
                            haploMat,recombFreqMat,predType,
                            postMeanAddEffects,postMeanDomEffects,
                            postVarCovarOfAddEffects=NULL,postVarCovarOfDomEffects=NULL,...){
   starttime<-proc.time()[3]

   # Before predicting variances
   # check for and remove SNPs that
   # won't segregate, i.e. are fixed in parents
   ### hopes to save time / mem
   x<-colSums(rbind(haploMat[grep(paste0("^",sireID,"_"),rownames(haploMat)),],
                    haploMat[grep(paste0("^",damID,"_"),rownames(haploMat)),]))
   segsnps2keep<-names(x[x>0 & x<4])

   if(length(segsnps2keep)>0){
      recombFreqMat<-recombFreqMat[segsnps2keep,segsnps2keep,drop=F]
      haploMat<-haploMat[,segsnps2keep,drop=F]
      postMeanAddEffects<-purrr::map(postMeanAddEffects,~.[segsnps2keep])
      postMeanDomEffects<-purrr::map(postMeanDomEffects,~.[segsnps2keep])
      if(predType=="PMV"){
         postVarCovarOfAddEffects<-postVarCovarOfAddEffects[segsnps2keep,segsnps2keep,drop=F]
         postVarCovarOfDomEffects<-postVarCovarOfDomEffects[segsnps2keep,segsnps2keep,drop=F]
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
      # out<-tibble::tibble(VarComp=c("VarA","VarD"),
      #                     VPM=c(vpm_m2a,vpm_m2d),
      #                     PMV=ifelse(predType=="PMV",c(pmv_m2a,pmv_m2d),c(NA,NA)),
      #                     Nsegsnps=c(length(segsnps2keep),NA),
      #                     totcomputetime=c(totcomputetime,NA))

      out<-dplyr::bind_rows(tibble::tibble(VarComp="VarA",VPM=c(vpm_m2a),PMV=ifelse(predType=="PMV",pmv_m2a,NA)),
                            tibble::tibble(VarComp="VarD",VPM=c(vpm_m2d),PMV=ifelse(predType=="PMV",pmv_m2d,NA))) %>%
         dplyr::mutate(Nsegsnps=length(segsnps2keep),
                       totcomputetime=c(totcomputetime,NA))

      print(paste0("Variances predicted for family: ",sireID,"x",damID,"- took ",round(totcomputetime/60,3)," mins"))
   } else {
      out<-tibble::tibble(VarComp=c("VarA","VarD"),
                          VPM=c(0,0),
                          PMV=c(0,0),
                          Nsegsnps=c(0,0),
                          computetime=c(NA,NA))
      print(paste0("Variances predicted for family: ",sireID,"x",damID,"- had no segregating SNPs"))
   }
   return(out)
}

#' Multi-trait, multi-cross prediction of one additive genetic variance (or covariance) among full-siblings.
#'
#' User specifies a trait variance (or trait-trait covariance) to predict. Wrapper around `predOneCrossVarA()` for predicting multiple families. Input a data.frame of crosses to predict.
#' NOTE: Marker effects should represent allele substitution effects.
#' Either by fitting an additive-only model OR an genotypic-partitioned additive+dominance model,
#' with allele sub effects computed as a+d(q-p), where q and p are allele freqs in the training pop. used.
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2 string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param predType string, "VPM" or "PMV". Choose option "VPM" if you have REML marker effect estimates (or posterior-means from MCMC) one set of marker effect estimates per trait. Variance of posterior means is faster but the alternative predType=="PMV" is expected to be less biassed. PMV requires user to supply a (probably LARGE) variance-covariance matrix of effects estimates.
#' @param postMeanAlleleSubEffects list of named vectors (or column matrices) with the allele substitution marker effects (can posterior-mean effects from MCMC _or_ from REML, if setting predType="PMV".
#' @param AlleleSubEffectList Only if setting predType="PMV". List of of ALLELE SUBSTITUTION EFFECT matrices. One matrix per trait. Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs. If users effects are from REML or posterior-means MCMC, matrices will be of dimension 1 x N SNP. If users chose predType="PMV", each matrix will be dimension N thinned-MCMC sample x N SNP.
#' @param ncores If ncores set > 1 parallelizes across families, but beware it is memory intensive and options(future.globals.maxSize=___) may need to be adjusted.
#' @param ...
#'
#' @return list with two elements, "predictedfamvars" contains a tibble with all predictions for all requested families, "totcomputetime" gives the time taken to compute one var. parameter across all families, at the given ncores.
#' @export
#'
#' @examples
predCrossVarsA<-function(Trait1,Trait2,CrossesToPredict,predType,
                         haploMat,recombFreqMat,
                         postMeanAlleleSubEffects,
                         AlleleSubEffectList=NULL,
                         ncores,...){

   starttime<-proc.time()[3]
   if(predType=="PMV"){
      # Posterior Sample Variance-Covariance Matrix of Marker Effects
      postVarCovarOfAlleleSubEffects<-(1/(nrow(AlleleSubEffectList[[Trait1]])-1))*crossprod(AlleleSubEffectList[[Trait1]],AlleleSubEffectList[[Trait2]])
      #rm(AddEffectList); gc()
   } else {
      postVarCovarOfAlleleSubEffects<-NULL;
   }

   require(furrr); options(mc.cores=ncores); plan(multisession)
   predictedfamvars<-CrossesToPredict %>%
      dplyr::mutate(predVars=furrr::future_pmap(.,
                                                predOneCrossVarA,
                                                Trait1=Trait1,Trait2=Trait2,
                                                haploMat=haploMat,recombFreqMat=recombFreqMat,
                                                predType=predType,
                                                postMeanAlleleSubEffects=postMeanAlleleSubEffects,
                                                postVarCovarOfAlleleSubEffects=postVarCovarOfAlleleSubEffects))
   totcomputetime<-proc.time()[3]-starttime
   print(paste0("Done predicting fam vars. ",
                "Took ",round((totcomputetime)/60,2),
                " mins for ",nrow(predictedfamvars)," families"))
   predictedfamvars<-list(predictedfamvars=predictedfamvars,totcomputetime=totcomputetime)
   return(predictedfamvars)
}

#' Multi-trait, multi-cross prediction of one additive _and_ one dominance genetic variance (or covariance) among full-siblings.
#'
#' User specifies a trait variance (or trait-trait covariance) to predict. Wrapper around `predOneCrossVarAD()` for predicting multiple families. Input a data.frame of crosses to predict.
#'
#' @param Trait1 string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param Trait2 string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param predType string, "VPM" or "PMV". Choose option "VPM" if you have REML marker effect estimates (or posterior-means from MCMC) one set of marker effect estimates per trait. Variance of posterior means is faster but the alternative predType=="PMV" is expected to be less biassed. PMV requires user to supply a (probably LARGE) variance-covariance matrix of effects estimates.
#' @param postMeanAddEffects list of named vectors (or column matrices) with the additive marker effects (can posterior-mean effects from MCMC _or_ from REML, if setting predType="VPM".
#' @param postMeanDomEffects list of named vectors (or column matrices) with the dominance marker effects (can posterior-mean effects from MCMC _or_ from REML, if setting predType="VPM".
#' @param AddEffectList Only if setting predType="PMV". List of ADDITIVE effect matrices. One matrix per trait. Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs. If users effects are from REML or posterior-means MCMC, matrices will be of dimension 1 x N SNP. If users chose predType="PMV", each matrix will be dimension N thinned-MCMC sample x N SNP.
#' @param DomEffectList Only if setting predType="PMV". List of DOMINANCE effect matrices. One matrix per trait. Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs. If users effects are from REML or posterior-means MCMC, matrices will be of dimension 1 x N SNP. If users chose predType="PMV", each matrix will be dimension N thinned-MCMC sample x N SNP.
#' @param ncores If ncores set > 1 parallelizes across families, but beware it is memory intensive and options(future.globals.maxSize=___) may need to be adjusted.
#' @param ...
#'
#' @return list with two elements, "predictedfamvars" contains a tibble with all predictions for all requested families, "totcomputetime" gives the time taken to compute one var. parameter across all families, at the given ncores.
#' @export
#'
#' @examples
predCrossVarsAD<-function(Trait1,Trait2,CrossesToPredict,predType,
                          haploMat,recombFreqMat,
                          postMeanAddEffects,postMeanDomEffects,
                          AddEffectList=NULL,DomEffectList=NULL,
                          ncores,...){
   starttime<-proc.time()[3]
   if(predType=="PMV"){
      # Posterior Sample Variance-Covariance Matrix of Marker Effects
      postVarCovarOfAddEffects<-(1/(nrow(AddEffectList[[Trait1]])-1))*crossprod(AddEffectList[[Trait1]],AddEffectList[[Trait2]])
      postVarCovarOfDomEffects<-(1/(nrow(DomEffectList[[Trait1]])-1))*crossprod(DomEffectList[[Trait1]],DomEffectList[[Trait2]])
      #rm(AddEffectList,DomEffectList); gc()
   } else {
      postVarCovarOfAddEffects<-NULL;
      postVarCovarOfDomEffects<-NULL;
   }

   require(furrr); options(mc.cores=ncores); plan(multisession)
   predictedfamvars<-CrossesToPredict %>%
      dplyr::mutate(predVars=furrr::future_pmap(.,
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


#' Multi-trait, multi-cross prediction of additive genetic variance and covariances among full-siblings.
#'
#' Wrapper around `predCrossVarsA()` for predicting all trait variances and trait-trait covariances across multiple families .
#' Input a data.frame of crosses to predict.
#' NOTE: Marker effects should represent allele substitution effects.
#' Either by fitting an additive-only model OR an genotypic-partitioned additive+dominance model,
#' with allele sub effects computed as a+d(q-p), where q and p are allele freqs in the training pop. used.
#'
#' @param outprefix string, prefix for *.rds file to be written to disk with output. DEFAULT = NULL (no disk write)
#' @param outpath string, path to disk location where files should be written. Can be left DEFAULT = NULL (no disk write)
#' @param predType Default = "VPM". string, "VPM" or "PMV" Choose option "VPM" if you have REML marker effect estimates (or posterior-means from MCMC) one set of marker effect estimates per trait. Variance of posterior means is faster but the alternative predType=="PMV" is expected to be less biassed. PMV requires user to supply a (probably LARGE) variance-covariance matrix of effects estimates.
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param AlleleSubEffectList List of of ALLELE SUBSTITUTION EFFECT matrices. One matrix per trait. Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs. If users effects are from REML or posterior-means MCMC, matrices will be of dimension 1 x N SNP. If users chose predType="PMV", each matrix will be dimension N thinned-MCMC sample x N SNP.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param ncores Default = 1. If ncores set > 1 parallelizes across families, but beware it is memory intensive and options(future.globals.maxSize=___) may need to be adjusted.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runMtCrossVarPredsA<-function(outprefix=NULL,outpath=NULL,predType="VPM",
                              CrossesToPredict,AlleleSubEffectList,
                              haploMat,recombFreqMat,ncores=1,...){
   starttime<-proc.time()[3]
   traits<-names(AlleleSubEffectList)
   parents<-union(CrossesToPredict$sireID,
                  CrossesToPredict$damID)
   # Center posterior distribution of effects
   ## on posterior mean across MCMC samples
   AlleleSubEffectList<-purrr::map(AlleleSubEffectList,~scale(.,center = T, scale = F))

   ## Get the posterior mean effects vectors
   postMeanAlleleSubEffects<-purrr::map(AlleleSubEffectList,~attr(.,which = "scaled:center"))

   ## Predict trait (co)variances
   if(length(traits)>1){
      varcovars<-dplyr::bind_rows(tibble::tibble(Trait1=traits,Trait2=traits), # trait variances
                                  combn(traits,2,simplify = T) %>% # covariances
                                     t(.) %>% #
                                     `colnames<-`(.,c("Trait1","Trait2")) %>%
                                     tibble::as_tibble(.)) } else {
                                        varcovars<-tibble::tibble(Trait1=traits,Trait2=traits)
                                     }

   haploMat<-haploMat[sort(c(paste0(parents,"_HapA"),paste0(parents,"_HapB"))),]

   if(predType!="PMV"){
      AlleleSubEffectList<-NULL;
   }

   varcovars<-varcovars %>%
      dplyr::mutate(varcomps=purrr::pmap(.,predCrossVarsA,CrossesToPredict=CrossesToPredict,predType=predType,
                                         AlleleSubEffectList=AlleleSubEffectList,
                                         haploMat=haploMat,recombFreqMat=recombFreqMat,
                                         postMeanAlleleSubEffects=postMeanAlleleSubEffects,ncores=ncores))

   totcomputetime<-proc.time()[3]-starttime
   varcovars<-list(varcovars=varcovars,
                   totcomputetime=totcomputetime)
   if(!is.null(outprefix) & !is.null(outpath)){
      saveRDS(varcovars,file=here::here(outpath,paste0(outprefix,"_predVarAndCovarBVs.rds")))
   }
   return(varcovars)
}


#' Multi-trait, multi-cross prediction of additive _and_ dominance genetic variances _and_ covariances among full-siblings.
#'
#' User specifies a trait variance (or trait-trait covariance) to predict.
#' Wrapper around `predCrossVarsAD()` for predicting all trait variances and trait-trait covariances across multiple families for both additive _and_ dominance components.
#' Input a data.frame of crosses to predict.
#' NOTE: Marker effects should represent allele substitution effects.
#'
#' @param outprefix string, prefix for *.rds file to be written to disk with output. DEFAULT = NULL (no disk write)
#' @param outpath string, path to disk location where files should be written. Can be left DEFAULT = NULL (no disk write)
#' @param predType Default = "VPM". string, "VPM" or "PMV" Choose option "VPM" if you have REML marker effect estimates (or posterior-means from MCMC) one set of marker effect estimates per trait. Variance of posterior means is faster but the alternative predType=="PMV" is expected to be less biassed. PMV requires user to supply a (probably LARGE) variance-covariance matrix of effects estimates.
#' @param CrossesToPredict data.frame or tibble, col/colnames: sireID, damID. sireID and damID must both be in the haploMat.
#' @param AddEffectList List of ADDITIVE effect matrices. One matrix per trait. Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs. If users effects are from REML or posterior-means MCMC, matrices will be of dimension 1 x N SNP. If users chose predType="PMV", each matrix will be dimension N thinned-MCMC sample x N SNP.
#' @param DomEffectList List of DOMINANCE effect matrices. One matrix per trait. Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs. If users effects are from REML or posterior-means MCMC, matrices will be of dimension 1 x N SNP. If users chose predType="PMV", each matrix will be dimension N thinned-MCMC sample x N SNP.
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
   parents<-union(CrossesToPredict$sireID,
                  CrossesToPredict$damID)

   # Center posterior distribution of effects
   ## on posterior mean across MCMC samples
   AddEffectList<-purrr::map(AddEffectList,~scale(.,center = T, scale = F))
   DomEffectList<-purrr::map(DomEffectList,~scale(.,center = T, scale = F))

   ## Get the posterior mean effects vectors
   postMeanAddEffects<-purrr::map(AddEffectList,~attr(.,which = "scaled:center"))
   postMeanDomEffects<-purrr::map(DomEffectList,~attr(.,which = "scaled:center"))

   ## Predict trait (co)variances
   if(length(traits)>1){
      varcovars<-dplyr::bind_rows(tibble::tibble(Trait1=traits,Trait2=traits), # trait variances
                                  combn(traits,2,simplify = T) %>% # covariances
                                     t(.) %>% #
                                     `colnames<-`(.,c("Trait1","Trait2")) %>%
                                     tibble::as_tibble(.)) } else {
                                        varcovars<-tibble::tibble(Trait1=traits,Trait2=traits)
                                     }

   haploMat<-haploMat[sort(c(paste0(parents,"_HapA"),paste0(parents,"_HapB"))),]

   if(predType!="PMV"){
      AddEffectList<-NULL;
      DomEffectList<-NULL;
   }

   varcovars<-varcovars %>%
      dplyr::mutate(varcomps=purrr::pmap(.,predCrossVarsAD,CrossesToPredict=CrossesToPredict,predType=predType,
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



