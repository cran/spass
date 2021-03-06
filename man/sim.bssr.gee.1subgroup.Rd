% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.bssr.gee.1subgroup.r
\name{sim.bssr.gee.1subgroup}
\alias{sim.bssr.gee.1subgroup}
\title{Simulation of a longitudinal one subgroup design with internal pilot Study}
\source{
\code{sim.bssr.gee.1subgroup} uses code contributed by Roland Gerard Gera.
}
\usage{
sim.bssr.gee.1subgroup(
  nsim = 1000,
  alpha = 0.05,
  tail = "both",
  beta = 0.2,
  delta = c(0.1, 0.1),
  vdelta = c(0.1, 0.1),
  sigma_pop = c(3, 3),
  vsigma_pop = c(3, 3),
  tau = 0.5,
  rho = 0.25,
  vrho = 0.25,
  theta = 1,
  vtheta = 1,
  Time = 0:5,
  rec.at = 0.5,
  k = 1,
  model = 1,
  V = diag(rep(1, length(Time))),
  OD = 0,
  vdropout = rep(0, length(Time)),
  missingtype = "none",
  vmissingtype = "none",
  seed = 2015
)
}
\arguments{
\item{nsim}{number of simulation runs.}

\item{alpha}{level (type I error) to which the hypothesis is tested.}

\item{tail}{which type of test is used, e.g. which quartile und H0 is calculated}

\item{beta}{type II error (power=1-beta) to which an alternative should be proven.}

\item{delta}{vector of true treatment effects, c(overall population, inside subgroup).}

\item{vdelta}{vector of treatment effects to be proven, c(overall population, inside subgroup).}

\item{sigma_pop}{vector of true standard deviations of the treatment effects, c(overall population, subgroup).}

\item{vsigma_pop}{vector of assumed standard deviations, c(overall population, inside subgroup).}

\item{tau}{subgroup prevalence.}

\item{rho}{true correlation coefficient between two adjacent timepoints}

\item{vrho}{initial expectation of the correlation coefficient between two adjacent timepoints}

\item{theta}{true correlation absorption coefficient if timepoints are farther apart}

\item{vtheta}{expected correlation absorption coefficient if timepoints are farther apart}

\item{Time}{vector of measured timepoints}

\item{rec.at}{blinded sample size review is performed after \code{rec.at}*\eqn{100\%} subjects of the initial sample size calculation.}

\item{k}{sample size allocation factor between groups: see 'Details'.}

\item{model}{which of the two often revered statistical models should be used?: see 'Details'.}

\item{V}{working covariance matrix.}

\item{OD}{overall dropout measured at last timepoint}

\item{vdropout}{vector of expected dropouts per timepoint if missingness is to be expected}

\item{missingtype}{true missingtype underlying the missingness}

\item{vmissingtype}{initial assumptions about the missingtype underlying the missingness}

\item{seed}{set seed value for the simulations to compare results.}
}
\value{
\code{sim.bssr.1subgroup} returns a data.frame containing the mean and variance of recalculated sample sizes within the control group and treatment group respectively and the achieved simulated power along with all relevant parameters.
}
\description{
Given estimates of the treatment effects to be proven, the variances, and the prevalence,
\code{sim.bssr.gee.1subgroup} calculates an initial sample size and performs a blinded sample size recalculation
after a pre-specified number of subjects have been enrolled. Each observation is simulated and a final analysis executed.
Several variations are included, such as different approximations or sample size allocation.
}
\details{
This function combines sample size estimation, blinded sample size re-estimation and analysis in a design with a subgroup within a full population where we want to test for treatment effects between a control and a treatment group.
The required sample size for the control and treatment group to prove an existing
alternative \code{delta} with a specified power 1-\code{beta} when testing the global null hypothesis \eqn{H_0: \Delta_F=\Delta_S=0} to level \code{alpha} is calculated prior to the study and then recalculated in an internal pilot study.

For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the
sample size allocation factor, i.e. \eqn{k = n_T/n_C}.
}
\examples{
sim.bssr.gee.1subgroup(nsim = 5,missingtype = "intermittened")

}
\seealso{
\code{sim.bssr.gee.1subgroup} makes use of \code{\link{n.gee.1subgroup}}, \code{\link{bssr.gee.1subgroup}}, and \code{\link{r.gee.1subgroup}}.
}
