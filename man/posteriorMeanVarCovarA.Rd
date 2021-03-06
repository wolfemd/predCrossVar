% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predCrossVar.R
\name{posteriorMeanVarCovarA}
\alias{posteriorMeanVarCovarA}
\title{posteriorMeanVarCovarA}
\usage{
posteriorMeanVarCovarA(
  Trait1,
  Trait2,
  AddEffectList,
  postMeanAddEffects,
  genoVarCovarMat
)
}
\arguments{
\item{Trait1}{string, label for Trait1. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.}

\item{Trait2}{string, label for Trait2. When Trait1==Trait2 computes the genomic variance of the trait, when Trait1!=Trait2 computes the genomic covariance between traits.}

\item{AddEffectList}{list of ADDITIVE effect matrices, one matrix per trait, Each element of the list is named with a string identifying the trait and the colnames of each matrix are labelled with snpIDs.}

\item{postMeanAddEffects}{list of named vectors (or column matrices) with the posterior mean ADDITIVE marker effects.}

\item{genoVarCovarMat}{variance-covariance matrix of marker genotypes, i.e. measuring linkage disequilibrium}
}
\value{
tibble with estimators for one variance or co-variance parameter.
}
\description{
For an additive only model, compute multiple types of genomic estimator for a single variance or covariance component.
Posterior mean variance or co-variance (PMV) and also the variance of posterior means (VPM) estimators are returned.
Both M1, i.e. "Method 1", which does not account for LD in the population and M2 or "Method 2", which does, are returned.
}
