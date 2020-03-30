#' Title
#'
#' @param M
#' @param type
#'
#' @return
#' @export
#'
#' @examples
kinship<-function(M,type){
      # Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers)
      # M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
      # Two types of dominance matrix, as described in
      # Vitezica et al. 2013. Genetics (Genotypic and Classical)
      # Both predict identically.
      # Difference in meaning of variance components.
      # type == "Add" should match A.mat() / Van Raden, Method 2
      M<-round(M)
      freq <- colMeans(M,na.rm=T)/2
      P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
      if(type=="Add"){
            Z <- M-2*P
            varD<-sum(2*freq*(1-freq))
            K <- tcrossprod(Z)/ varD
            return(K)
      }
      if(type=="Dom"){
            W<-M;
            W[which(W==1)]<-2*P[which(W==1)];
            W[which(W==2)]<-(4*P[which(W==2)]-2);
            W <- W-2*(P^2)
            varD<-sum((2*freq*(1-freq))^2)
            D <- tcrossprod(W) / varD
            return(D)
      }
}

#' Title
#'
#' @param pop
#' @param SP
#' @param type
#'
#' @return
#' @export
#'
#' @examples
makeGRM<-function(pop, SP, type){
      grm<-kinship(M=pullSnpGeno(pop=pop, simParam=SP),
                   type=type)
      return(grm)
}

#' Title
#'
#' @param M
#'
#' @return
#' @export
#'
#' @examples
centerDosage<-function(M){
      # Modified from my "kinship" function
      # Converts a dosage matrix into a centered-dosage
      # Purpose is for whole-genome regressions like rrBLUP
      # Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers)
      # M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
      freq <- colMeans(M,na.rm=T)/2
      P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
      Z <- M-2*P
      return(Z)
}

#' Title
#'
#' @param M
#'
#' @return
#' @export
#'
#' @examples
dose2domDev<-function(M){
      # Modified from my "kinship" function
      # Converts a dosage matrix into a matrix of centered dominance deviations
      # Purpose is for whole-genome regressions like rrBLUP
      # Assumes SNPs in M coded 0, 1, 2 (requires rounding dosages to integers)
      # M is Nind x Mrow, numeric matrix, with row/columanes to indicate SNP/ind ID
      # Vitezica et al. 2013. Genetics ("Classical")
      freq <- colMeans(M,na.rm=T)/2
      P <- matrix(rep(freq,nrow(M)),byrow=T,ncol=ncol(M))
      W<-M;
      W[which(W==1)]<-2*P[which(W==1)];
      W[which(W==2)]<-(4*P[which(W==2)]-2);
      W <- W-2*(P^2)
      return(W)
}