% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bssr.gee.1subgroup.r
\name{bssr.gee.1subgroup}
\alias{bssr.gee.1subgroup}
\title{Blinded Sample Size Recalculation for longitudinal data in a One Subgroup Design}
\source{
\code{bssr.gee.1subgroup} uses code contributed by Roland Gerard Gera.
}
\usage{
bssr.gee.1subgroup(
  alpha,
  tail = "both",
  beta = NULL,
  delta,
  estsigma,
  tau = 0.5,
  k = 1
)
}
\arguments{
\item{alpha}{level (type I error) to which the hypothesis is tested.}

\item{tail}{which type of test is used, e.g. which quartile und H0 is calculated.}

\item{beta}{type II error (power=1-beta) to which an alternative should be proven.}

\item{delta}{vector of estimated treatment effect in overall and sub population, c(overall population, only subpopulation).}

\item{estsigma}{vector of re-estimated standard deviations, c(full population, subpopulation). See 'Details'.}

\item{tau}{ratio between complementary F/S and sub-population S.}

\item{k}{treatment allocation factor between groups: see 'Details'.}
}
\value{
\code{bssr.gee.1subgroup} returns a list containing the recalculated sample sizes along with all relevant parameters. Use \code{\link{summary.bssrest}} for a structured overview.
}
\description{
Given re-estimations from an Internal Pilot Study (IPS), \code{bssr.GEE.1subgroup} re-estimates required sample size given the re-estimated nuisance parameters are given. \code{bssr.gee.1subgroup} is a wrapper for \code{n.gee.1subgroup} where the re-estimation of the variances can be highly dependable on the user and should be supplied separately. see "detail" for more information.
}
\details{
This function provides a simple warped for \code{n.gee.1subgroup} where instead of initial assumptions, reestimated nuisance parameter are used.
For more information see \code{n.gee.1subgroup}.
Required samplesize to test alternative \code{delta} with specified power 1-\code{beta} when testing the global null hypothesis \eqn{H_0: \beta_3^F=\beta_3^S=0} to level \code{alpha} is estimated. When testing outcomes have variance \code{estsigma}.

For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group respectively, the argument \code{k} is the
sample size allocation factor, i.e. \eqn{k = n_T/n_C} and \code{tau} represents the ratio of the sub-population.
}
\examples{
estimate<-bssr.gee.1subgroup(alpha=0.05,beta=0.2,delta=c(0.1,0.1),estsigma=c(0.8,0.4),tau=0.4, k=1)
summary(estimate)
}
\seealso{
\code{\link{n.gee.1subgroup}} for sample size calculation prior to a trial and \code{estimcov} how the re-estimate nuisance parameters. See \eqn{sim.gee} for a working example for an initial sample size estimation and a re-estimation mid trial.
}
