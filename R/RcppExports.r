# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

minFunc <- function(x, daten, dataNA) {
    .Call('spass_minFunc', PACKAGE = 'spass', x, daten, dataNA)
}

minFuncMult <- function(x, daten, dataNA, n) {
    .Call('spass_minFuncMult', PACKAGE = 'spass', x, daten, dataNA, n)
}

minFuncBlinded <- function(x, daten, dataNA, n, delta) {
    .Call('spass_minFuncBlinded', PACKAGE = 'spass', x, daten, dataNA, n, delta)
}

mlFirst <- function(y, groupE, groupC, nE, nC, tpE, tpC, type) {
    .Call('spass_mlFirst', PACKAGE = 'spass', y, groupE, groupC, nE, nC, tpE, tpC, type)
}

mlFirstOneGroup <- function(y, groupC, nC, tpC, type) {
    .Call('spass_mlFirstOneGroup', PACKAGE = 'spass', y, groupC, nC, tpC, type)
}

mlSecond <- function(rho, y, groupE, groupC, nE, nC, tpE, tpC, type) {
    .Call('spass_mlSecond', PACKAGE = 'spass', rho, y, groupE, groupC, nE, nC, tpE, tpC, type)
}

mlSecondOneGroup <- function(rho, y, groupC, nC, tpC, type) {
    .Call('spass_mlSecondOneGroup', PACKAGE = 'spass', rho, y, groupC, nC, tpC, type)
}

mlFirstGrad <- function(y, groupE, groupC, nE, nC, tpE, tpC, type) {
    .Call('spass_mlFirstGrad', PACKAGE = 'spass', y, groupE, groupC, nE, nC, tpE, tpC, type)
}

mlFirstHObs <- function(y, groupE, groupC, nE, nC, tpE, tpC, type) {
    .Call('spass_mlFirstHObs', PACKAGE = 'spass', y, groupE, groupC, nE, nC, tpE, tpC, type)
}

mlFirstHExp <- function(y, kf, tp, type) {
    .Call('spass_mlFirstHExp', PACKAGE = 'spass', y, kf, tp, type)
}

mlFirstJObs <- function(y, groupE, groupC, nE, nC, tpE, tpC, type) {
    .Call('spass_mlFirstJObs', PACKAGE = 'spass', y, groupE, groupC, nE, nC, tpE, tpC, type)
}

mlFirstJExp <- function(y, rho, kf, tp, type, approx = 50L) {
    .Call('spass_mlFirstJExp', PACKAGE = 'spass', y, rho, kf, tp, type, approx)
}

mlFirstBlinded <- function(y, group, n, tp, type, theta, k) {
    .Call('spass_mlFirstBlinded', PACKAGE = 'spass', y, group, n, tp, type, theta, k)
}

mlSecondBlinded <- function(rho, y, group, n, tp, type, k) {
    .Call('spass_mlSecondBlinded', PACKAGE = 'spass', rho, y, group, n, tp, type, k)
}

