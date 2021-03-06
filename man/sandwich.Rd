% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sandwich.r
\name{sandwich}
\alias{sandwich}
\title{Calculate the robust covariance estimator for GEE given an}
\source{
\code{sandwich} computes the asymptotic sandwich covariance estimator and uses code contributed by Roland Gerard Gera.
}
\usage{
sandwich(
  yCov,
  D,
  V,
  correctionmatrix,
  missing = rep(0, dim(yCov)[[2]]),
  missingtype = c("none", "monotone", "intermittened")
)
}
\arguments{
\item{yCov}{\eqn{yCov} matrix containing either an estimation for the covariance between timepoints or an empirical covariance matrix itsel. see 'Details'.}

\item{D}{\eqn{D} denotes the mean matrix of all entries of \eqn{\Delta\mu_i/\delta\beta}, where this is the average over all \code{i} patiens. see 'Details'.}

\item{V}{\eqn{V} denotes the working covariance matrix.  see 'Details'.}

\item{correctionmatrix}{As of this version this matrix is needed to correct some calculations. see 'Details' to see for more details and how to correctly select matrices.}

\item{missing}{vector which denotes the probability to experience a dropout at each timepoint. If \code{missingtype} is \code{"none"} then all entries are 0.}

\item{missingtype}{String which describes the type of missingness occuring in the data. \code{none} if no missingnes occured, \code{"monotone"} if missing was monotone and \code{"intermittened"} if the missingness was independent across all timepoints.}
}
\value{
\code{sandwich} returns the robust covariance estimator of regression coefficients which are implicitly defined by \code{D}.
}
\description{
\code{sandwich} calculates the covariance structure between timepoints given matrices \code{yCov}, \code{D},\code{V} and \code{correctionmatrix}. This is done to be able to account for missingness in the Data.
}
\details{
\code{yCov} is either empirical or the estimated covariance-matrix between timepoints which is needed to calculate the sandwich estimator. This matrix can either be generated by estimating the empirical covariance matrix using existing data or by using function \code{gen_cov_cor} to calculate a estimation for the covariance.

\code{D} denotes the estimation of \eqn{n^-1* \sum_i^N \Delta\mu_i/\delta\beta}, which means that \eqn{D=E(D_i)}. As of yet this has the unfortunate side effect that E(D_i%*%D_i) is not equal to E(D_i)%*%E(D_i). Entries that differ between those two formulas are equated using \code{correctionsmatrix}.
}
\examples{
#Let's assume we wish to calculate the robust variance estimator for  equation
#\eqn{y_{it}=\beta_0+\beta_1*I_{treat}+\beta_2*t+\beta_3*I _{treat}*t+\epsilon_{it}}.
#Furthermore we use the identitiy matrix as the working covariance matrix.
#The chance to get treatment is 60 percent and the observed timerange ranges from 0:5.

  ycov = gen_cov_cor(var = 3,rho = 0.25,theta = 1,Time = 0:5,cov = TRUE)
  D = matrix(c(1,0.6,0,0,
               1,0.6,1,0.6,
               1,0.6,2,1.2,
               1,0.6,3,1.8,
               1,0.6,4,2.4,
               1,0.6,5,3.0),nrow=4)

 D=t(D)
 V=diag(1,length(0:5))
 #We correct entries where E(D_i \%*\% D_i) is unequal to E(D_i)\%*\%E(D_i) (D \%*\% D).
 correctionmatrix=matrix(c(1,1,1,1,1,1/0.6,1,1/0.6,1,1,1,1,1,1/0.6,1,1/0.6),nrow=4)
 missingtype = "none"

 robust=sandwich(yCov=ycov,D=D,V=V,missingtype=missingtype,correctionmatrix=correctionmatrix)
 robust

}
\references{
Liang Kung-Yee, Zeger Scott L. (1986); Jung Sin-Ho, Ahn Chul (2003); Wachtlin Daniel Kieser Meinhard (2013)
}
