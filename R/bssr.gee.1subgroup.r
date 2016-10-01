#' @title Blinded Sample Size Recalculation for longitudinal data in a One Subgroup Design
#' @description Given data from an Internal Pilot Study (IPS), \code{bssr.GEE.1subgroup} given the reestimated nuisance parameteres are calculated. \code{bssr.gee.1subgroup} is a wraper for \code{n.gee.1subgroup} because the reestimation of the variances can be highly dependable on the user and should be done seperatly. see "detail" for more information on that.
#'
#' @param alpha    level (type I error) to which the hypothesis is tested.
#' @param beta     type II error (power=1-beta) to which an alternative should be proven.
#' @param delta    vector of regression coefficients values which shall be proven, c(allcomers, only subpopulation).
#' @param k        sample size allocation factor between groups: see 'Details'.
#' @param tail     which type of test is used, e.g. which quartile und H0 is calculated
#' @param estsigma reestimated vector of assymptotic standard deviations.
#' @param tau      ration between F/S and S
#'
#' @details
#' This function provides a simple warper for \code{n.gee.1subgroup} where instead of initial assumptions blind estimated nuisance parameter inserted.
#' For information see \code{n.gee.1subgroup}.
#' alternative \code{delta} with a specified power 1-\code{beta} when testing the global null hypothesis \eqn{H_0: \beta_3^F=\beta_3^S=0} to level \code{alpha} is calculated.
#'
#' The data matrix \code{data} should have as many columns as observed timepoints: first column first observed timepoint. As of now the timepoints must be equispaced to calculate the correct intra-subject covariance-matrix. Entries can be \code{NA}. See \code{r.gee.1subgroup.r} for more information.
#'
#' For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the
#' sample size allocation factor, i.e. \eqn{k = n_T/n_C}.
#'
#'
#' @return \code{bssr.gee.1subgroup} returns a list containing the recalculated required sample size within the control group and treatment group along with all relevant parameters. Use \code{\link{summary.bssrest}} for a structured overview.
#'
#' @source \code{bssr.gee.1subgroup} uses code contributed by Roland Gerard Gera.
#'
#' @seealso \code{\link{n.gee.1subgroup}} for sample size calculation prior to the trial, \code{r.gee.1subgroup} how list data should look like and \code{estimcov} how the reestimation of nuisance parameters works. See \eqn{sim.gee} for an enxample for an initial sample size estimation and reestimation to see the functions working in junction.
#'
#' @examples
#' estimate<-bssr.gee.1subgroup(alpha=0.05,beta=0.2,delta=c(0.1,0.1),estsigma=c(0.8,0.4),tau=0.4, k=1)
#' summary(estimate)
#' @export

bssr.gee.1subgroup <- function(alpha, tail="both",beta=NULL, delta, estsigma, tau=0.5, k = 1) { 
  if(alpha >=1 || alpha <=0) stop("Wrong Type-I-Error specefied")
  if(beta >=1 || beta <=0) stop("Wrong Type-II-Error specefied")
  
  result = n.gee.1subgroup(delta=delta, 
                           sigma=estsigma,
                           alpha=alpha,
                           tail=tail,
                           beta=beta, 
                           tau=tau , 
                           k = k)
  
  class(result)<-"bssrest"
  return(result)
}
