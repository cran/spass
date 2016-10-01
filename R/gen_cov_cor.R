#' @title Generation of covariance- or correlation-matrices
#' 
#' @description Generate covariance- or correlation-matrices given the parameters \code{var}, \code{rho}, \code{theta} for the covariance structure, \code{Time} for the observed timepoints and \code{cov=TRUE} if a covariance or \code{cov=FALSE} if a correlation-matrix is to be generated.
#'
#' @param var       variance at each timepoint 
#' @param rho       correlation between two adjacent timepoints 1 timeunit appart
#' @param theta     variable specifying the type of the correlation structure: see 'Details'
#' @param Time      list with time measures which are used to generate the covariance- or correlation-structure: see 'Details'
#' @param cov       TRUE/FALSE statement which determines if a covariance- or a correlation-matrix is generated.
#'
#'
#' @details The function \code{gen_cov_cor} is used to generate either a covariance- or a correlation-matrix. Given vector \code{Time} and parameters \code{var}, \code{rho} and \code{theta} the covariance or correlation between two timepoints is described: The way the correlation between two timepoints is described by
#' 
#' \eqn{cov(Time[i],Time[j])=var*(rho^(abs(Time[i]-Time[j])^theta))}
#'
#'for covariance and 
#'
#'\eqn{corr(Time[i],Time[j])=rho^(abs(Time[i]-Time[j])^theta)}
#'
#'for correlation. 
#'[The above sentence would be better is we write the following sentence:]
#'[[The following two equations are used to calculate
#'the covariance and the correlation between two timepoints, respectively:
#' cov(Time[i],Time[j])=var*(rho^(abs(Time[i]-Time[j])^theta))
#' corr(Time[i],Time[j])=rho^(abs(Time[i]-Time[j])^theta) ]]
#'
#'Which of the two formulas is used depends on \code{cov}
#' @return \code{gen_cov_cor} returns the covariance- or correlation-matrix, depending if \code{cov} was \code{TRUE} or \code{FALSE}.
#'
#' @source \code{gen_cov_cor} uses code contributed by Roland Gerard Gera
#'
#'  @seealso \code{\link{r.gee.1subgroup}} for information on the generated longitudinal data and \code{\link{n.gee.1subgroup}} for the calculation of
#' initial sample sizes for longitudinal GEE-models and \code{\link{bssr.gee.1subgroup}} for blinded
#' sample size reestimation within a trial See \code{\link{estimcov}} for more information on the used minimization algorithms.
#'
#' @examples
#' #Generate a covariance-matrix with measurements at Baseline and at timeunit 1,1.5,2 
#' #and 5 (hours,day,months,years, etc.) 
#'
#' covar<-gen_cov_cor(var=3,rho=0.25,theta=1,Time=c(0,1,1.5,2,5),cov=TRUE)
#' covar
#' 
#' #Generate a correlation-matrix
#'
#' corr<-gen_cov_cor(rho=0.25,theta=1,Time=c(0,1,1.5,2,5),cov=FALSE)
#' corr
#' @export


gen_cov_cor <- function(var=1,rho,theta,Time,cov=TRUE){
  
  
  tp=length(Time)
  # Generate an empty mold for your matrix
  answer=matrix(0,ncol=tp,nrow=tp)
  
  # generate covariance
  if (cov){
    for (i in 1:tp){
      for (j in 1:tp){
        # prohibits that 0^0=1; we define her 0^0=0
        if (i-j==0 && theta==0) hoch=0
        else hoch=abs(Time[i]-Time[j])^theta
        answer[i,j]=var*(rho^(hoch))
      }
    }
    
  # generate correlation
  } else {
    for (i in 1:tp){
      for (j in 1:tp){
        # prohibits that 0^0=1; we define her 0^0=0
        if (i-j==0 && theta==0) hoch=0
        else hoch=abs(Time[i]-Time[j])^theta
        answer[i,j]=rho^(hoch)
      }
    }
  }
  
  return(answer)
  
}