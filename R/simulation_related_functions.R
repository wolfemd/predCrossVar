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
