# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

algoritmo_probit <- function(y, nit, beta, f, sigma2, alpha, C0 = 100000, a = 0.001, b = 0.001, sdpropbeta = 0.05, sdpropbeta2 = 0.075, sdpropf = 0.4, sdpropsigma2 = 0.05, display_progress = TRUE) {
    .Call('_fatcat_algoritmo_probit', PACKAGE = 'fatcat', y, nit, beta, f, sigma2, alpha, C0, a, b, sdpropbeta, sdpropbeta2, sdpropf, sdpropsigma2, display_progress)
}

