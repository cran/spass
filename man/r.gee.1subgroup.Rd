% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r.gee.1subgroup.r
\name{r.gee.1subgroup}
\alias{r.gee.1subgroup}
\title{Generate dataset of normal distributed repeated observations in a one subgroup design}
\source{
\code{r.gee.1subgroup} uses code contributed by Roland Gerard Gera
}
\usage{
r.gee.1subgroup(n, reg, sigma, rho, theta, tau, k, Time, OD)
}
\arguments{
\item{n}{overall sample size for the overall population}

\item{reg}{list containing coefficients \deqn{\beta_0} to \deqn{\beta_0} for complementary population, \code{reg[[1]]} and subpopulation, \code{reg[[2]]}: see 'Details'.}

\item{sigma}{vector with standard deviations for generated observations c(complementary population, subpopulation).}

\item{rho}{variable used together with \code{theta} to describe correlation between two adjacent timepoints: see 'Details'.}

\item{theta}{variable used together with \code{rho} to describe correlation between two adjacent timepoints: see 'Details'.}

\item{tau}{subgroup prevalence.}

\item{k}{sample size allocation factor between treatment groups: see 'Details'.}

\item{Time}{list of timepoints \eqn{t} that have to be generated: see 'Details'.}

\item{OD}{percentage of observed overall dropout at last timepoint: see 'Details'.}
}
\value{
\code{r.gee.1subgroup} returns a list with 7 different entries. Every Matrix rows are the simulated subjects and the columns are the observed time points.

The first list element is a vector containing subject ids.
The second element contains a matrix with the outcomes of a subject with row being the subjects and columns being the measuring-timepoints
Elements 3 to 5 return matrices with the information of which patients have baseline-measurements, which patients belong to treatment and which to control and what are the observed timepoints for each patient respectively.
The sixth entry returns a matrix which contains the residuals of each measurement.
The seventh entry returns the sub-population identification.
}
\description{
\code{r.gee.1subgroup} generates data for a design with one subgroup within a full population. Each baseline-observation is normal distributed with mean \deqn{\beta_0} in placebo group and \deqn{\beta_0+\beta_1} in treatment group.
Measurements after baseline have mean \deqn{\beta_0+\beta_2*t} in placebo group and \deqn{\beta_0+\beta_1+\beta_2*t+\beta_3*t} in treatment group where \deqn{t} is the measurement time. Whether the effect can be found solely in the subgroup or additionally a certain amount outside of the subgroup can be specified as well as a potential different covariance-structure within subgroup and in the complementary subgroup.
}
\details{
For \code{reg}\code{list}(c(\eqn{\beta_0^F\S,\beta_1^F\S,\beta_2^F\S,\beta_3^F\S}), c(\eqn{\beta_0^S,\beta_1^S,\beta_2^S,\beta_3^S})) and variances \code{sigma}=(\eqn{\sigma_F\S, \sigma_S}) function \code{r.gee.1subgroup} generates data given correlation-variables \eqn{\rho} and \eqn{\theta} as follows (and let t=0 be the baseline measurement):

Placebo group - complementary population \eqn{y_{it}=N(\beta_0+\beta_2*t,\sigma_F\S)},
Placebo group - within subgroup \eqn{y_{it}=N(\beta_0+\beta_2*t,\sigma_S)},
Treatment group - complementary population \eqn{y_{it}=N(\beta_0+\beta_1+\beta_2*t+\beta_3*t,\sigma_F\S)},
Treatment group - within subgroup \eqn{y_{it}=N(\beta_0+\beta_1+\beta_2*t+\beta_3*t,\sigma_S)}.
Correlation between measurements - \eqn{corr(\epsilon_it,\epsilon_io)=\rho^{(t-o)^\theta}}

 Argument \code{k} is the sample size allocation factor, i.e. the ratio between control and treatment. Let \eqn{n_C} and \eqn{n_T} denote sample sizes of control and treatment groups respectively, then \eqn{k = n_T/n_C}.

Argument \code{Time} is the vector denoting all measuring-times, i. e. every value for \eqn{t}.

Argument \code{OD} sets the overall dropout rate observed at the last timepoint. For \code{OD}=0.5, 50 percent of all observation had a dropout event at some point. If a subject experienced a dropout the starting time of the dropout is equally distributed over all timepoints.
}
\examples{

set.seed(2015)
dataset<-r.gee.1subgroup(n=200, reg=list(c(0,0,0,0.1),c(0,0,0,0.1)), sigma=c(3,2.5),
tau=0.5, rho=0.25, theta=1, k=1.5, Time=c(0:5), OD=0)
dataset
}
