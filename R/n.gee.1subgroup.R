#' @title Sample Size estimation for longitudinal GEE Models when testing 1 coefficient 
#' @description \code{n.gee.1subgroup} calculates the required sample size for proving a desired alternative when testing regression coefficients in the full or subpopulation. See 'Details' for more information.
#'
#' @param alpha   level (type I error) to which the hypothesis is tested.
#' @param beta    type II error (power=1-beta) to which an alternative should be proven.
#' @param delta   vector of regression coefficients values which shall be proven, c(allcomers, subpopulation).
#' @param sigma   vector of assymptotic standard diviation of regressors, c(full population, subpopulation).see 'Details'
#' @param tau     subgroup prevalence.
#' @param k       sample size allocation factor between control and treatment: see 'Details'.
#' @param nmax    maximum total sample size.
#' @param npow    calculates power of a test if \code{npow} is a sample size
#' @param tail    which type of test is used, e.g. which quartile und H0 is calculated
#' 
#' @details
#' This function performs sample size estimation in a design with a subgroup nested within a full population. To calculate the required sample size when testing only one regressor (e.g. effect of treatment*time) one needs to input the expected value of the regressor under alternative, \code{delta}, and the expected asymptotic variance of that regressor, \code{sigma}. 
#' The power for the global null hypothesis is given by 1-\code{beta} and \code{alpha} specifies the false positve level for rejecting \eqn{H_0: \Delta_F=\Delta_S=0} to level \code{alpha}, where in our context \eqn{\Delta_F and \Delta_S} normaly represent regressioncoefficents and \eqn{\sigma^2} their variance.
#'
#' For sample sizes \eqn{n_C} and \eqn{n_T} of the control and treatment group, respectively, the argument \code{k} is the
#' sample size allocation factor, i.e. \eqn{k = n_T/n_C}.
#'
#'
#' @return \code{n.gee.1subgroup} returns the required sample size within the control group and treatment group.
#'
#' @source \code{n.gee.1subgroup} uses code contributed by Roland Gerard Gera.
#'
#' @seealso \code{\link{bssr.1subgroup}} for blinded sample size reestimation within a running trial and \code{\link{sandwich}} for estimating asymptotic covarianc mtrices in GEE models.
#'
#' @examples
#' #Calculate required sample size to correctly reject with
#' #80% probability when testing the global Nullhypothesis H_0: Delta_F=Delta_S = 0
#' #assuming the coefficient in and outside of the subgroup is Delta=c(0.1,0,1) with a 
#' #subgroup-prevalence of tau=0.4.
#' #The assymptotic variances in and outside of the subgroup are unequal, sigma=c(0.8,0.4).
#'
#' estimate<-n.gee.1subgroup(alpha=0.05,beta=0.2,delta=c(0.1,0.1),sigma=c(0.8,0.4),tau=0.4, k=1)
#' summary(estimate)
#' 
#' #Now we want to estimate the power our study would have, 
#' #if we know the effect in and outside the subgroup, as 
#' #well as asymptotic variance of the regressors. Here we 
#' #estimate that only 300 Patiens total can be recruited.
#' #All other parameters are the same as those above. 
#' 
#'
#' n.gee.1subgroup(alpha=0.05,delta=c(0.1,0.1),sigma=c(0.8,0.4),tau=0.4, k=1, npow=300)
#' @export

# Berechnung der Sample.Size ----------------------------------------------
n.gee.1subgroup <- function(alpha, 
                            tail="both",
                            beta=NULL, 
                            delta, 
                            sigma, 
                            tau = 0.5 ,
                            k = 1,
                            npow=NULL,
                            nmax=Inf){
  
  if((is.null(beta) & is.null(npow)) | (!is.null(beta) & !is.null(npow))){
    stop("Either beta OR npow have to be specified")
  }
  
  # Load required packages
  #mvtnorm <- require(mvtnorm, quietly = TRUE)
  
  # If normal distribution should be used
  Correl = matrix(c(1,sqrt(tau),sqrt(tau),1),nrow=2)
  
  #calkulate Critical values from Null
  CritVal = qmvnorm(1-alpha, mean=c(0,0), corr=Correl, tail=tail)$quantile
  
  # Over all Samplesize
  n=min(1000,nmax)
  n_s=round(n*tau)
  
  # if the power is to be returned, and not the sample size
  if (!is.null(npow)) {
    power=1-pmvnorm(mean=c(sqrt(npow)*delta[1]/sigma[1], sqrt(npow*tau)*delta[2]/sigma[2]), corr=Correl,lower=c(-Inf,-Inf),upper=c(CritVal,CritVal))[1]
    names(power)<-c("The power with given sample size is:")
    return(power)
  }
  
  # Initial power for n=nmax
  power=1-pmvnorm(mean=c(sqrt(n)*delta[1]/sigma[1],
                                  sqrt(n_s)*delta[2]/sigma[2]),
                           corr=Correl,lower=c(-Inf,-Inf),upper=c(CritVal,CritVal))[1]
  
  # First reduce Power until you are below 0.8
  while(power > 0.80){
    n=n-1
    n_s=round(n*tau)
    # New power
    power=1-pmvnorm(mean=c(sqrt(n)*delta[1]/sigma[1],
                                    sqrt(n_s)*delta[2]/sigma[2]),
                             corr=Correl,lower=c(-Inf,-Inf),upper=c(CritVal,CritVal))[1]
    
  }
  
  # then increas Power until you are above 0.8 the first time
  while(power<0.80){
    if (n>nmax){
      print("!!Warning!!")
      cat("calculated total sample size exceeds nmax = ", nmax, ";\n")
      cat("expected power using nmax subjects: ", round(power,2))
      power=0.8
    } else{
      n=n+1
      n_s=round(n*tau)
      
      power=1-pmvnorm(mean=c(sqrt(n)*delta[1]/sigma[1],
                                      sqrt(n_s)*delta[2]/sigma[2]),
                               corr=Correl,lower=c(-Inf,-Inf),upper=c(CritVal,CritVal))[1]
    }
  }
  # Save your solutions
  n_cont=round(n/(k+1))
  n<-c(n_cont,n-n_cont)
  
  names(n)<-c("Control", "Treatment")
  result<-list(n=n, alpha=alpha, beta=beta, delta=delta, sigma_reg=sigma, tau=tau, k=k, model="normal.gee.1subgroup")
  class(result)<-"ssest"
  return(result)
}