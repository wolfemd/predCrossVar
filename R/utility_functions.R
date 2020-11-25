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
#' Converts a dosage matrix into a matrix of centered dominance deviations.
#' This sets up the "classical" (aka "Statistical") partition of additive dominance in terms of breeding values and dom. deviations.
#' See function dose2domDevGenotypic() to set-up the "genotypic" (aka "biological") partition in terms of genotypic effects.
#' See Vitezica et al. 2013.
#' Also see Varona et al. 2018. https://doi.org/10.3389/fgene.2018.00078
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


#' dose2domDevGenotypic
#'
#' Converts a dosage matrix into a matrix of centered dominance deviations.
#' This sets up the "genotypic" (aka "biological") partition of additive dominance in terms of their genotypic effects instead of in terms of breeding values or dominance deviations.
#' See function dose2domDev() to set-up the "statistical" (aka "classical") partition in terms of genotypic effects.
#' See Vitezica et al. 2013.
#' Also see Varona et al. 2018. https://doi.org/10.3389/fgene.2018.00078
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return
#' @export
#'
#' @examples
#' domDev<-dose2domDevGenotypic(M)
dose2domDevGenotypic<-function(M){
   freq <- colMeans(M,na.rm=T)/2
   P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
   W<-M; W[which(W==2)]<-0;
   W <- W-(2*P*(1-P))
   return(W)
}

#' getPropHom
#'
#' Compute the per-individual proportion homozygous.
#' For example, to use as a predictor for a directional dominance model.
#'
#' @param M dosage matrix. Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers). M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
#'
#' @return vector of per-individual proportion homozygous
#' @export
#'
#' @examples
getPropHom<-function(M){
   W<-M; W[which(W==2)]<-0;
   # f = 1 − h/N,
   # where N is the number of SNPs
   f<-1-(rowSums(W)/ncol(W))
   return(f)
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
      X<-haploMat[paste0(parentGID,c("_HapA","_HapB")),,drop=F]
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

#' crosses2predict
#'
#' Make a data.frame of all pairwise matings given a vector of parent IDs.
#' Include selfs. No reciprocal crosses, i.e. use as male == use as female.
#' Diagonal and upper-triangle of mating matrix.
#'
#' @param parents
#'
#' @return tibble, two columns, sireID and damID, all pairwise crosses (see details).
#' @export
#'
#' @examples
crosses2predict<-function(parents){
      CrossesToPredict<-matrix(NA,nrow=length(parents),ncol=length(parents))
      CrossesToPredict[upper.tri(CrossesToPredict,diag = T)]<-1
      rownames(CrossesToPredict)<-colnames(CrossesToPredict)<-parents
      CrossesToPredict<-CrossesToPredict %>%
            as.data.frame %>%
            tibble::rownames_to_column(var = "sireID") %>%
            tidyr::pivot_longer(cols = (-sireID), names_to = "damID", values_to = "keep") %>%
            dplyr::filter(keep==1) %>%
            dplyr::select(-keep)
      return(CrossesToPredict)
}

#' intensity
#'
#' Compute the standardized selection intensity from the proportion selection.
#'
#' @param propSel proportion selection
#'
#' @return
#' @export
#'
#' @examples
intensity<-function(propSel){ dnorm(qnorm(1-propSel))/propSel } # standardized selection intensity
