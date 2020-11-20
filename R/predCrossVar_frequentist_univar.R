#' genomicVarA_m1
#'
#' Compute standard genomic variance term, what Lehermeier called Method 1 or "M1".
#' LD is not accounted for. SNP effects are assumed _i.i.d._.
#'
#' @param snpVarA vector of _additive_ SNP effects estimate. Ideally vector is named with SNP IDs and matches the order of related SNP matrices.
#' @param freq vector of allele frequencies. Must be in same order as snp effects
#' @return additive genomic variance component
#' @export
#'
#' @examples
genomicVarA_m1<-function(snpVarA,freq){
      sum((2*freq*(1-freq)))*snpVarA
}

#' genomicVarD_m1
#'
#' Compute standard genomic variance term, what Lehermeier called Method 1 or "M1".
#' LD is not accounted for. SNP effects are assumed _i.i.d._.
#'
#' @param snpVarD vector of _dominance_ SNP effects estimate. Ideally vector is named with SNP IDs and matches the order of related SNP matrices.
#' @param freq vector of allele frequencies. Must be in same order as snp effects
#' @return additive genomic variance component
#' @export
#'
#' @examples
genomicVarD_m1<-function(snpVarD,freq){
      sum((2*freq*(1-freq))^2)*snpVarD
}



#' getM2varcomp
#'
#' Computes the genomic variance component accounting account for linkage disequilibrium, which Lehermeier called Method 2 or "M2".
#'
#' @param effects vector of SNP effects
#' @param varcovarmat square, symmetric var-covar matrix, estimator of LD between loci
#' @param type string, "add" or "dom". If "dom" will square varcovarmat to get correct varcomp.
#'
#' @details Make sure **effects** and **varcovarmat** are in same order.
#' This might use _alot_ of RAM and take _alot_ of time as it involves very large matrix operations.
#' May be worth computing in an R session using multi-threaded BLAS
#' @return additive or dominance genomic variance estimate accounting account for linkage disequilibrium
#' @export
#'
#' @examples
#' Z<-centerDosage(M)
#' varcovarmat<-genoVarCovarMatFunc(Z)
#' VarAddM2<-getM2varcomp(effects,varcovarmat,"add")
getM2varcomp<-function(effects,varcovarmat,type){
      if(type=="dom"){ varcovarmat<-varcovarmat^2 }
      vgM2<-t(effects)%*%varcovarmat%*%effects
      return(vgM2)
}


#' predCrossVarA function
#'
#' Function to predict the additive and dominance variances in a full-sibling family based on marker effects from a Additive+Dominance genome-wide SNP-BLUP model.
#'
#' @param sireID string, Sire genotype ID. Needs to correspond to renames in haploMat
#' @param damID string, Dam genotype ID. Needs to correspond to renames in haploMat
#' @param addEffects column matrix of _additive_ SNP effects estimate with rownames == SNP_IDs and matches the order of related SNP matrices.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param ...
#'
#' @details SNP_IDs must match: names(addEffects)==names(domEffects)==colnames(haploMat)==rownames(recombFreqMat)==colnames(recombFreqMat).
#'
#' @return a tibble with values for predicted add/dom variances as well as compute time and memory usage stats
#' @export
#'
#' @examples
predCrossVarA<-function(sireID,damID,addEffects,
                        haploMat,recombFreqMat,...){
      starttime<-proc.time()[3]
      # check for and remove SNPs fixed in parents
      ### hopes to save time / mem
      x<-colSums(rbind(haploMat[grep(paste0("^",sireID,"_"),rownames(haploMat)),],
                       haploMat[grep(paste0("^",damID,"_"),rownames(haploMat)),]))
      segsnps2keep<-names(x[x>0 & x<4])
      if(length(segsnps2keep)>0){
            # drop=F safegaurds against matrix->vector conversion if only 1 seg. locus (column)
            haploMat<-haploMat[,segsnps2keep,drop=F]
            recombFreqMat<-recombFreqMat[segsnps2keep,segsnps2keep,drop=F]
            addEffects<-addEffects[segsnps2keep,,drop=F]
            # sire and dam LD matrices
            sireLD<-calcGameticLD(sireID,recombFreqMat,haploMat)
            damLD<-calcGameticLD(damID,recombFreqMat,haploMat)
            progenyLD<-sireLD+damLD
            rm(recombFreqMat,haploMat,sireLD,damLD); gc()
            # predict additive variance
            predAdd<-t(addEffects)%*%(progenyLD)%*%addEffects
            totcomputetime<-proc.time()[3]-starttime
            gc() #gcend
            # output predicted add and dom vars for current family
            out<-tibble(varA=predAdd,
                        segsnps=list(segsnps2keep),
                        computetime=totcomputetime)
            print(paste0("Variances predicted for family: ",sireID,"x",damID,"- took ",round(totcomputetime/60,3)," mins"))
      } else {
            out<-tibble(varA=0,
                        segsnps=list(),
                        computetime=0)
            print(paste0("Variances predicted for family: ",sireID,"x",damID,"- had no segregating SNPs"))
      }
      return(out)
}


#' predCrossVarAD function
#'
#' Function to predict the additive and dominance variances in a full-sibling family based on marker effects from a Additive+Dominance genome-wide SNP-BLUP model.
#'
#' @param sireID string, Sire genotype ID. Needs to correspond to renames in haploMat
#' @param damID string, Dam genotype ID. Needs to correspond to renames in haploMat
#' @param addEffects column matrix of _additive_ SNP effects estimate with rownames == SNP_IDs and matches the order of related SNP matrices.
#' @param domEffects column matrix of _dominance_ SNP effects estimate with rownames == SNP_IDs and matches the order of related SNP matrices.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param ...
#'
#' @details SNP_IDs must match: names(addEffects)==names(domEffects)==colnames(haploMat)==rownames(recombFreqMat)==colnames(recombFreqMat).
#'
#' @return a tibble with values for predicted add/dom variances as well as compute time and memory usage stats
#' @export
#'
#' @examples
predCrossVarAD<-function(sireID,damID,addEffects,domEffects,
                         haploMat,recombFreqMat,...){
      #gcbegin<-gc(reset = T)
      starttime<-proc.time()[3]
      # check for and remove SNPs fixed in parents
      ### hopes to save time / mem
      x<-colSums(rbind(haploMat[grep(paste0("^",sireID,"_"),rownames(haploMat)),],
                       haploMat[grep(paste0("^",damID,"_"),rownames(haploMat)),]))
      segsnps2keep<-names(x[x>0 & x<4])
      if(length(segsnps2keep)>0){
            # drop=F safegaurds against matrix->vector conversion if only 1 seg. locus (column)
            haploMat<-haploMat[,segsnps2keep,drop=F]
            recombFreqMat<-recombFreqMat[segsnps2keep,segsnps2keep,drop=F]
            addEffects<-addEffects[segsnps2keep,,drop=F]
            domEffects<-domEffects[segsnps2keep,,drop=F]
            # sire and dam LD matrices
            sireLD<-calcGameticLD(sireID,recombFreqMat,haploMat)
            damLD<-calcGameticLD(damID,recombFreqMat,haploMat)
            progenyLD<-sireLD+damLD
            rm(recombFreqMat,haploMat,sireLD,damLD); gc()
            # predict additive variance
            predAdd<-t(addEffects)%*%(progenyLD)%*%addEffects
            # predict dominance variance
            predDom<-t(domEffects)%*%(progenyLD)^2%*%domEffects
            totcomputetime<-proc.time()[3]-starttime
            gc() #gcend
            # output predicted add and dom vars for current family
            out<-tibble(varA=predAdd,
                        varD=predDom,
                        segsnps=list(segsnps2keep),
                        computetime=totcomputetime)
            print(paste0("Variances predicted for family: ",sireID,"x",damID,"- took ",round(totcomputetime/60,3)," mins"))
      } else {
            out<-tibble(varA=0,
                        varD=0,
                        segsnps=list(),
                        computetime=0)
            print(paste0("Variances predicted for family: ",sireID,"x",damID,"- had no segregating SNPs"))
      }
      return(out)
}


#' runCrossVarPredsA
#'
#' Function to compute predicted add + dom vars for an entire pedigree.
#' If outprefix and outpath are supplied, writes output to disk so impatient users can see results.
#'
#' @param outprefix
#' @param outpath
#' @param ped pedigree data.frame, cols: sireID, damID. sireID and damID must both be in the haploMat.
#' @param addEffects vector of _additive_ SNP effects estimate. Ideally vector is named with SNP IDs and matches the order of related SNP matrices.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param ncores If ncores set > 1 parallelizes across families, but beware it is memory intensive and options(future.globals.maxSize=___) may need to be adjusted.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runCrossVarPredsA<-function(outprefix=NULL,outpath=NULL,
                            ped,addEffects,
                            haploMat,recombFreqMat,ncores=1,...){
      require(furrr); options(mc.cores=ncores); plan(multisession)
      timestart<-proc.time()[3]
      predictedfamvars<-ped %>%
            dplyr::mutate(predVars=future_map2(sireID,damID,
                                               ~predCrossVarA(sireID=.x,damID=.y,
                                                              addEffects=addEffects,
                                                              haploMat=haploMat,
                                                              recombFreqMat=recombFreqMat))) %>%
            tidyr::unnest(predVars)
      endtime<-proc.time()[3];
      print(paste0("Done predicting fam vars. ",
                   "Took ",round((endtime-timestart)/60,2),
                   " mins for ",nrow(predictedfamvars)," families"))
      if(!is.null(outprefix) & !is.null(outpath)){
            saveRDS(predictedfamvars,
                    file=here::here(outpath,paste0(outprefix,"_predCrossVars_ModelAD.rds")))
      }
      return(predictedfamvars)
}


#' runCrossVarPredsAD
#'
#' Function to compute predicted add + dom vars for an entire pedigree.
#' If outprefix and outpath are supplied, writes output to disk so impatient users can see results.
#'
#' @param outprefix
#' @param outpath
#' @param ped pedigree data.frame, cols: sireID, damID. sireID and damID must both be in the haploMat.
#' @param addEffects vector of _additive_ SNP effects estimate. Ideally vector is named with SNP IDs and matches the order of related SNP matrices.
#' @param domEffects vector of _dominance_ SNP effects estimate. Ideally vector is named with SNP IDs and matches the order of related SNP matrices.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param ncores If ncores set > 1 parallelizes across families, but beware it is memory intensive and options(future.globals.maxSize=___) may need to be adjusted.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
runCrossVarPredsAD<-function(outprefix=NULL,outpath=NULL,
                             ped,addEffects,domEffects,
                             haploMat,recombFreqMat,ncores=1,...){
      require(furrr); options(mc.cores=ncores); plan(multisession)
      timestart<-proc.time()[3]
      predictedfamvars<-ped %>%
            dplyr::mutate(predVars=future_map2(sireID,damID,
                                               ~predCrossVarAD(sireID=.x,damID=.y,
                                                               addEffects=addEffects,
                                                               domEffects=domEffects,
                                                               haploMat=haploMat,
                                                               recombFreqMat=recombFreqMat))) %>%
            tidyr::unnest(predVars)
      endtime<-proc.time()[3];
      print(paste0("Done predicting fam vars. ",
                   "Took ",round((endtime-timestart)/60,2),
                   " mins for ",nrow(predictedfamvars)," families"))
      if(!is.null(outprefix) & !is.null(outpath)){
            saveRDS(predictedfamvars,
                    file=here::here(outpath,paste0(outprefix,"_predCrossVars_ModelAD.rds")))
      }
      return(predictedfamvars)
}