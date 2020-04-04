#' kinship function
#'
#' Function to create additive and dominance genomic relationship matrices from biallelic dosages.
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID.
#' @param type string, "add" or "dom". type="add" gives same as rrBLUP::A.mat(), i.e. Van Raden, Method 1. type="dom" gives classical parameterization according to Vitezica et al. 2013.
#'
#' @return square symmetic genomic relationship matrix
#' @export
#'
#' @examples
#' K<-kinship(M,"add")
kinship<-function(M,type){
      M<-round(M)
      freq <- colMeans(M,na.rm=T)/2
      P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
      if(type=="add"){
            Z <- M-2*P
            varD<-sum(2*freq*(1-freq))
            K <- tcrossprod(Z)/ varD
            return(K)
      }
      if(type=="dom"){
            W<-M;
            W[which(W==1)]<-2*P[which(W==1)];
            W[which(W==2)]<-(4*P[which(W==2)]-2);
            W <- W-2*(P^2)
            varD<-sum((2*freq*(1-freq))^2)
            D <- tcrossprod(W) / varD
            return(D)
      }
}

#' makeKinship function
#'
#' function to create a additive or dominance kinship matrix from a AlphaSimR pop-class object. A wrapper for the kinship() function.
#'
#' @param pop pop-class object from AlphaSimR
#' @param SP simulation parameters from AlphaSimR
#' @param type string, "add" or "dom". type="add" gives same as rrBLUP::A.mat(), i.e. Van Raden, Method 1. type="dom" gives classical parameterization according to Vitezica et al. 2013.
#' @param snpsORqtl string, "snps" (default) or "qtl". Determines which loci are used to compute the GRM.
#'
#' @return
#' @export
#'
#' @examples
#'  K<-makeKinship(pop, SP, type="add")
#'
makeKinship<-function(pop, SP, type, snpsORqtl="snps"){
      if(snpsORqtl=="snps"){
            M<-pullSnpGeno(pop=pop, simParam=SP) }
      if(snpsORqtl=="qtl"){
            M<-pullQtlGeno(pop=pop, simParam=SP) }
      grm<-kinship(M=M,type=type)
      return(grm)
}

#' centerDosage
#'
#' Centers dosage matrix, e.g. for use in whole-genome regressions like rrBLUP
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#'
#' @return
#' @export
#'
#' @examples
#' centeredM<-centerDosage(M)
centerDosage<-function(M){
      freq <- colMeans(M,na.rm=T)/2
      P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
      Z <- M-2*P
      return(Z)
}

#' dose2domDev
#'
#' Converts a dosage matrix into a matrix of centered dominance deviations
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return
#' @export
#'
#' @examples
#' domDev<-dose2domDev(M)
dose2domDev<-function(M){
      freq <- colMeans(M,na.rm=T)/2
      P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
      W<-M;
      W[which(W==1)]<-2*P[which(W==1)];
      W[which(W==2)]<-(4*P[which(W==2)]-2);
      W <- W-2*(P^2)
      return(W)
}

#' getAF
#'
#' get a vector of allele frequences from a dosage matrix
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return vector of allele frequencies, names = SNP IDs if in cols of M
#' @export
#'
#' @examples
getAF<-function(M){ colMeans(M,na.rm=T)/2 }

#' getMAF
#'
#' get a vector of _minor_ allele frequences from a dosage matrix
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return vector of _minor_ allele frequencies, names = SNP IDs if in cols of M
#' @export
#'
#' @examples
getMAF<-function(M){
      freq<-colMeans(M, na.rm=T)/2; maf<-freq;
      maf[which(maf > 0.5)]<-1-maf[which(maf > 0.5)]
      return(maf) }

#' maf_filter
#'
#' filter a dosage matrix by minor allele frequence.
#' get a vector of allele frequences from a dosage matrix
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#' @param thresh threshold value. Columns of M with maf<thresh will be removed
#' @return dosage matrix potentially with columns removed
#' @export
#'
#' @examples
maf_filter<-function(M,thresh){
      freq<-colMeans(M, na.rm=T)/2; maf<-freq;
      maf[which(maf > 0.5)]<-1-maf[which(maf > 0.5)]
      snps1<-M[,which(maf>thresh)];
      return(snps1) }

#' remove_invariant
#'
#' filter a dosage matrix, removing invariant markers. Removes e.g. cases where MAF=0.5 but all dosages == 1 (het).
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#' @param thresh threshold value. Columns of M with maf<thresh will be removed
#' @return dosage matrix potentially with columns removed
#' @export
#'
#' @examples
remove_invariant<-function(M){
      snps1<-M[ ,apply(M, 2, var) != 0]
      return(snps1) }

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

#' genoVarCovarMatFunc
#'
#' Compute the *p SNPs*-by-*p SNPs* variance-covariance matrix of SNP dosages.
#' This is an estimator of the LD between loci within a given population.
#'
#' @param Z column-centered matrix of SNP dosages. Assumes SNPs in Z were originally coded 0, 1, 2 were column centered.
#'
#' @return
#'  NOTE: this matrix is going to be big in practice.
#'  The *p SNPs*-by-*p SNPs* variance-covariance matrix of SNP dosages.
#'  may be worth computing in an R session using multi-threaded BLAS
#' @export
#'
#' @examples
#' Z<-centerDosage(M)
#' genoVarCovarMat<-genoVarCovarMatFunc(Z)
genoVarCovarMatFunc<-function(Z){
      SigmaM<-1/nrow(Z)*t(Z)%*%Z
      return(SigmaM)
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

#' genmap2recombfreq
#'
#' Compute the pairwise recombination frequencies between all loci from genetic map positions.
#'
#' @param m vector of centiMorgan-scale genetic positions. names(m) correspond to a SNP_ID. Since m potentially contains all chromosomes, sets recomb. freq. b/t chrom. to 0.5
#' @param nChr number of chromosomes
#'
#' @details names(m) must be formatted as "chr"_"id" with "chr" being integer. For example: 2_QTL1 for a locus on chr. 2.
#' May be worth computing in an R session using multi-threaded BLAS.
#' @return potentially really large matrix of pairwise recombination frequencies between loci
#' @export
#'
#' @examples
genmap2recombfreq<-function(m,nChr){
      d<-as.matrix(dist(m,upper=T,diag = T,method='manhattan'))
      c1<-0.5*(1-exp(-2*d))
      # Since m contains all chromosomes, set recomb. freq. b/t chrom. to 0.5
      for(i in 1:nChr){
            c1[grepl(paste0("^",i,"_"),rownames(c1)),!grepl(paste0("^",i,"_"),colnames(c1))]<-0.5
            c1[!grepl(paste0("^",i,"_"),rownames(c1)),grepl(paste0("^",i,"_"),colnames(c1))]<-0.5
      }
      return(c1)
}

#' calcGameticLD
#'
#' Function to compute gametic LD matrix for a single parent. Uses a matrix of recombination frequencies and the supplied haplotypes of the parent as in Bijma et al. 2020 and Lehermeier et al. 2017b (and others).
#'
#' @param parentGID string, the GID of an individual. Needs to correspond to renames in haploMat
#' @param recombFreqMat a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.
#' @param haploMat matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes
#'
#' @details Columns of haploMat and row/cols of recombFreqMat should be in same order. May be worth computing in an R session using multi-threaded BLAS.
#' @return Potentially really large matrix representing the LD between loci in gametes
#' @export
#'
#' @examples
calcGameticLD<-function(parentGID,recombFreqMat,haploMat){
      X<-haploMat[grep(paste0(parentGID,"_"),rownames(haploMat)),,drop=F]
      ### drop=F safegaurds against matrix->vector conversion if only 1 seg. locus (column)
      # Equiv versions (in case another is more efficient)
      # Version1
      D<-recombFreqMat*((0.5*crossprod(X))-tcrossprod(colMeans(X)))
      # Version2
      # p<-colMeans(X)
      # D<-recombFreqMat*((0.5*t(X)%*%X)-p%*%t(p))
      # Version3
      # D<-recombFreqMat*(crossprod(scale(X,scale = F))/2)
      return(D)
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
   x<-colSums(rbind(haploMat[grep(paste0(sireID,"_"),rownames(haploMat)),],
                    haploMat[grep(paste0(damID,"_"),rownames(haploMat)),]))
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
   x<-colSums(rbind(haploMat[grep(paste0(sireID,"_"),rownames(haploMat)),],
                    haploMat[grep(paste0(damID,"_"),rownames(haploMat)),]))
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
   require(furrr); options(mc.cores=ncores); plan(multiprocess)
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
      require(furrr); options(mc.cores=ncores); plan(multiprocess)
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

#' backsolveSNPeff
#'
#' From the GBLUP solutions and a centered SNP matrix backsolve SNP effects
#'
#' @param Z Centered marker matrix (dominance deviations must also be centered)
#' @param g The solutions (blups, i.e. GEBVs) from the GBLUP model
#'
#' @return matrix of SNP effects matching RR-BLUP / SNP-BLUP
#' @export
#'
#' @examples
#' A<-kinship(M,type="add")
#' trainingDF %<>% dplyr::mutate(ga=factor(as.character(id),
#'                                         levels=rownames(A)),
#'                               gd=ga)
#' gblup<-mmer(pheno~1,
#'             random=~vs(ga,Gu = A),
#'             weights=WT,
#'             data=trainingDF,verbose = T)
#' ga<-as.matrix(gblup$U$`u:ga`$pheno,ncol=1)
#' Za<-centerDosage(M)
#' snpeff<-backsolveSNPeff(Za,ga)
backsolveSNPeff<-function(Z,g){
      ZZt<-tcrossprod(Z)
      diag(ZZt)<-diag(ZZt)+1e-8
      bslvEffs<-crossprod(Z,solve(ZZt))%*%g
      return(bslvEffs)
}