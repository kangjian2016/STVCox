# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#'@importFrom Rcpp evalCpp
#'@useDynLib STVCox, .registration=TRUE
NULL

h_beta_t0_Rcpp <- function(tau, beta_t, alph) {
    .Call(`_STVCox_h_beta_t0_Rcpp`, tau, beta_t, alph)
}

#'@export
lplk_Rcpp <- function(par, alph, time, delta, z, Bs) {
    .Call(`_STVCox_lplk_Rcpp`, par, alph, time, delta, z, Bs)
}

#'@export
lplk_Rcpp_penalty <- function(par, alph, time, delta, z, Bs, rho) {
    .Call(`_STVCox_lplk_Rcpp_penalty`, par, alph, time, delta, z, Bs, rho)
}

lplk_Rcpp_min_beta <- function(par, alph, time, delta, z, Bs) {
    .Call(`_STVCox_lplk_Rcpp_min_beta`, par, alph, time, delta, z, Bs)
}

