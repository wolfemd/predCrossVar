% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predCrossVar_frequentist_univar.R
\name{predCrossVarAD}
\alias{predCrossVarAD}
\title{Single trait prediction of the additive _and_ dominance genetic variance in a one cross}
\usage{
predCrossVarAD(
  sireID,
  damID,
  addEffects,
  domEffects,
  haploMat,
  recombFreqMat,
  ...
)
}
\arguments{
\item{sireID}{string, Sire genotype ID. Needs to correspond to renames in haploMat}

\item{damID}{string, Dam genotype ID. Needs to correspond to renames in haploMat}

\item{addEffects}{column matrix of _additive_ SNP effects estimate with rownames == SNP_IDs and matches the order of related SNP matrices.}

\item{domEffects}{column matrix of _dominance_ SNP effects estimate with rownames == SNP_IDs and matches the order of related SNP matrices.}

\item{haploMat}{matrix of phased haplotypes, 2 rows per sample, cols = loci, {0,1}, rownames assumed to contain GIDs with a suffix, separated by "_" to distinguish haplotypes}

\item{recombFreqMat}{a square symmetric matrix with values = (1-2*c1), where c1=matrix of expected recomb. frequencies. The choice to do 1-2c1 outside the function was made for computation efficiency; every operation on a big matrix takes time.}

\item{...}{}
}
\value{
a tibble with values for predicted add/dom variances as well as compute time and memory usage stats
}
\description{
Function to predict the additive genetic _and_ dominance variances among full-siblings of a single, user-specified cross.
}
\details{
SNP_IDs must match: names(addEffects)==names(domEffects)==colnames(haploMat)==rownames(recombFreqMat)==colnames(recombFreqMat).
}
