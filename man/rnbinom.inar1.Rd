% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnbinom.inar1.r
\name{rnbinom.inar1}
\alias{rnbinom.inar1}
\title{Generate Time Series with Negative Binomial Distribution and Autoregressive Correlation Structure of Order One: NB-INAR(1)}
\source{
\code{rnbinom.inar1} computes a reparametrization of the NB-INAR(1) model by \emph{McKenzie 1986} using code contributed by Thomas Asendorf.
}
\usage{
rnbinom.inar1(n, size, mu, rho, tp)
}
\arguments{
\item{n}{number of observations. If length(n) > 1, the length is taken to be the number required.}

\item{size}{dispersion parameter (the shape parameter of the gamma mixing distribution). Must be strictly positive, need not be integer.}

\item{mu}{parametrization via mean: see 'Details'.}

\item{rho}{correlation coefficient of the underlying autoregressive correlation structure. Must be between 0 and 1.}

\item{tp}{number of observed time points.}
}
\value{
\code{rnbinom.inar1} returns a matrix of dimension \code{n} x \code{tp} with marginal negative binomial
distribution with mean \code{mu} and dispersion parameter \code{size}, and an autoregressive correlation structure
between time points.
}
\description{
\code{rnbinom.inar1} generates one or more independent time series following the NB-INAR(1) model. The generated data has negative binomial marginal distribution and an autoregressive covariance structure.
}
\details{
The generated marginal negative binomial distribution with mean \code{mu} = \eqn{\mu} and \code{size} = \eqn{\eta} has density
\deqn{(\mu/(\mu+\eta))^x \Gamma(x + \eta)/(\Gamma(x+1)\Gamma(\eta)) (\eta/(\mu+\eta))^\eta}
for \eqn{0 < \mu}, \eqn{0 < \eta} and \eqn{x=0, 1, 2, ...}.

Within the NB-INAR(1) model, the correlation between two time points \eqn{t} and \eqn{s} for \code{rho} = \eqn{\rho} is given through
\deqn{\rho^|t-s|}
for \eqn{0 \le \rho \le 1}.
}
\examples{
set.seed(8)
random<-rnbinom.inar1(n=1000, size=0.6, mu=2, rho=0.8, tp=6)
cor(random)

#Check the marginal distribution of time point 3
plot(table(random[,3])/1000, xlab="Probability", ylab="Observation")
lines(0:26, dnbinom(0:26, mu=2, size=0.6), col="red")
legend("topright",legend=c("Theoretical Marginal Distribution", "Observed Distribution"), 
col=c("red", "black"), lty=1, lwd=c(1,2))

}
\references{
McKenzie Ed (1986), Autoregressive Moving-Average Processes with Negative-Binomial and Geometric Marginal Distributions. \emph{Advances in Applied Probability} Vol. 18, No. 3, pp. 679-705.
}
