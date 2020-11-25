# predCrossVar

<!-- badges: start -->

<!-- badges: end -->

The **predCrossVar** package

contains a complete set of functions for the prediction of additive and dominance genetic variances and *co*-variances among full-siblings based on parents. **predCrossVar** enables the prediction of genetic variance on multi-trait selection indices. Built for diploid organisms with phased, chromosome- or linkage-group ordered biallelic marker data, and a centimorgan-scale genetic map.

## Installation

You can install **predCrossVar** [My GitHub](https://www.github.com/wolfemd/) with:

``` {.r}
devtools::install_github("wolfemd/predCrossVar", ref = 'master') 
```

## HIGHLIGHTS

-   Allows for parents to be of arbitrary heterozygosity/homozygosity (outbred or inbred)
-   Predicts the **additive** *and* **dominance** genetic variances in the $F_1$
-   Predicts genetic **variances** *and* ***co***-**variances**. Enables prediction of genetic variance on an multi-trait **selection index**.
-   Handles simple (one trait, one cross) predictions, but built for complex (multi-trait, many crosses) scenarios.
-   Single estimate of marker effects from REML or MCMC (posterior mean effects --\> predicts "variance of posterior means") supported
-   Posterior Mean Variance (PMV) *also* supported: the estimator of Lehermeier et al. 2017b which computes the predicted variance across a sample of marker effects, e.g. the thinned MCMC samples, usually stored on disk. For the multi-trait case, a multivariate Bayesian model is required as only a marker effects for each trait must be computed on the same Gibbs chain.

## Core functions

The package helps automate prediction for both simple (one trait, one cross) and complex (multi-trait, many crosses) scenarios.

| Single-trait Functions | REML or $MCMC^{VPM}$                                       |
|------------------------|------------------------------------------------------------|
| `predCrossVarA()`      | Predicts one cross variance                                |
| `runCrossVarPredsA()`  | Wraps around `predCrossVarA()` to predict multiple crosses |

Equivalent functions for an additive-plus-dominance models are: `predCrossVarAD()` --\> `runCrossVarPredsAD()` No function as of now to do dominance separately.

Note that these functions were developed early, and the multi-trait functions below should supersede these in function.

| Multi-trait Functions | REML, $MCMC^{VPM}$, $MCMC^{PMV}$ | 
|-----------------------|----------------------------------|
|`predOneCrossVarA()`   | Predicts one cross, one variance (or covariance), additive-only model |
| `predCrossVarsA()`    | wraps around `predOneCrossVarA()`. Predicts multiple crosses (potentially in parallel/multicore), one variance (or covariance), additive-only model | 
| `runMtCrossVarPredsA()` | wraps around `predCrossVarsA()`. Predicts the variances *and* covariances in the multi-trait case for a set of crosses, additive-only model. | 

Equivalent functions for an additive-plus-dominance models are: `predOneCrossVarAD()` --\> `predCrossVarsAD()` --\> `runMtCrossVarPredsAD()`. No function as of now to do only dominance.

There are also functions for the much less computationally challenging prediction of cross means (`predCrossMeanBVsOneTrait()` --\> `predCrossMeanBVs()`) for convenience.


