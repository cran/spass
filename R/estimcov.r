#' @title Estimation of variance, intra-subject-correlation and dropout 
#' 
#' @description \code{estimcov} estimates variance, intra-subject-correlation and dropout given empirical data.
#'
#' @param data         list of gathered data. The list must be consistent with the gernerated data of \code{\link{r.gee.1subgroup}}
#' @param Time         list with observed time points: see 'Details'
#' @param Startvalues  startvalues for the paramteres \code{var},\code{rho} and \code{theta}
#' @param stepwidth    vector of stepwidths which the optimisation-function should use
#' @param maxiter      value setting maximal amount of iterations for the optimisation algorithm
#' @param lower        lower bound for \code{var},\code{rho} and \code{theta}
#' @param upper        upper bound for \code{var},\code{rho} and \code{theta}
#'
#' @details Function \code{estimcov} fits a covariance-matrix with parameters \code{var},\code{rho} and \code{theta} (see \code{\link{gen_cov_cor}} for matrix generation) to an empirical covarince-matrix provided by \code{data}.
#'
#'
#' @return \code{estimcov} returns a list with two vectors. The first entry consists of a vector with estimations for c(\code{var},\code{rho},\code{theta})
#' while the second entry contains a vector, describing the empirical dropout-chance per timepoint.
#'
#' @source \code{estimcov} uses code contributed Roland Gerard Gera.
#'
#'
#' @seealso \code{\link{r.gee.1subgroup}} for information on the generated longitudinal data and \code{\link{n.gee.1subgroup}} for the calculation of
#' initial sample sizes for longitudinal GEE-models and \code{\link{bssr.gee.1subgroup}} for blinded
#' sample size reestimation within a trial. See \code{\link{gen_cov_cor}} for more information on the generation of covariance matrices.
#'
#' @examples
#' #Generate data from longitudinal-model
#' set.seed(2015)
#' dataset<-r.gee.1subgroup(n=300, reg=list(c(0,0,0,0.1),c(0,0,0,0.1)), sigma=c(3,2.5), tau=0.5, 
#' rho=0.25, theta=1, k=1.5, Time=c(0:5), OD=0.2)
#'
#' estimate<-estimcov(data=dataset,Time=c(0:5))
#' estimate
#' @export

# just testing!

estimcov <- function(data,
                       Time, 
                       Startvalues = c(3,0.5,1),
                       stepwidth=c(1e-3,1e-3,1e-3),
                       maxiter=10000,
                       lower = c(0.0001,0.0001,0.0001), 
                       upper = c(Inf,5,3)){
  
  # starting values for optimasation
  estimation = Startvalues 
  
  # get empiric covariance matrix from data
  CovMat=cov(data$y,use="pairwise.complete.obs")
  
  estimation = optim(par=estimation,
                     fn = optim.function,
                     CovMat=CovMat,
                     Time=Time,
                     method = "L-BFGS-B",
                     lower = lower, 
                     upper=  upper,
                     control = list(maxit=maxiter,ndeps=stepwidth))$par
  
  if (estimation[2]<0) estimation[2]=0
  if (estimation[3]<0) estimation[3]=0.0001
  
  #### estimation of dropout
  
  # How many supervised times exist
  m = dim(data$y)[2]
  
  # empty vector for dropout
  drop=c()
  
  # Estimate missing subjects
  for (i in 1:m){
    drop[i] = sum(is.na(data$y[,i]))/length(data$y[,i]) 
  }
  
  #integrate estimated parameters and estimated dropout in the list "answer" 
  
  answer=list(estimation=estimation,drop=drop)
  return(answer)
}

# function which is to be optimated
optim.function <- function(para,CovMat,Time){
  
  if (para[1]<=0) para[1]=0.0001
  if (para[2]<0) para[2]=0.0001
  fit = sum(abs(CovMat-gen_cov_cor(var=para[1],rho=para[2],theta=para[3],Time=Time)))
  return(fit)
}