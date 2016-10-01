#' @title Estimate the Robust covariance estimator for GEE (weigthed GEE if missing occures) of Regressor parameters
#' 
#' @description \code{sandwich} calculates the asymptotic regression covariance structure given matrices \code{yCov}, \code{D},\code{V}, \code{correctionmatrix} for further analyses and is a morte generalized, but also more complex version as \code{sandwich2}.
#'
#' @param yCov              \eqn{yCov} is the empirical or estimated Covariancematrix that we cen get from the outcomes. see 'Details'.
#' @param D                 \eqn{D} es the mean Matrix of all entries of  \eqn{\Delta\mu_i/\delta\beta}, where is the average over all \code{i} patiens. see 'Details'.
#' @param V                 \eqn{V} is the Working covariance matrix.  see 'Details'.
#' @param missing           a vector which describes the probabilety to experience a dropout at all observed timepoints. if missing is \code{"none"} then it is treated as if all entries are 0
#' @param missingtype       describes the type of missing that occured in tzhe data. Possebiletys range from \code{none} if there is no missing, to \code{"monotone"} if missing is monotone, aka dropout, and lastly \code{"intermittened"} if the missingness is independent across ale timepoints
#' @param correctionmatrix  a correctionamtrix that will correct misstakes. see 'Details' to see what these misstakes are und how to select correction matrices. see 'Details'.
#'
#' @details \code{yCov} is the either empirical or estimated intra-subject covariancematrix which is needed to calculate the sandwich (robust) covariance estimator. This matrix can either be achieved by estimating the empirical intra-subject covariance out of data or by using \code{gen_cov_cor} to calculate a estimation for the covariance.
#' 
#' \code{D} is the estimation of \eqn{n^-1* \sum_i^N \Delta\mu_i/\delta\beta}, so \eqn{D=E(D_i)}. But this is also source of an error which has to be corrected by \code{correctionsmatrix}. The error emerges when we calculate the "Bread" and "Meat" of the sandwichestimators.  Exemplary on the "Bread" we need to calculate \eqn{E(D_i \times V \times D_i)} wich is however enequal to \eqn{E(D_i)^t \times V \times E(D_i)} which we ARE calculating. \code{correctionmatrix} is now used to correct made misstakes so that \eqn{E(D_i)^t \times V \times E(D_i)}*\code{correctionsmatrix}=\eqn{E(D_i \times V \times D_i)}, which is still a point which we will improve on further itterations of the algorithm.
#'
#'
#' @return \code{sandwich} returns the sandwich (robust) covariance estimator of regression coefficients which are impicently defined by \code{D}.
#'
#' @source \code{sandwich} computes the asymptotic sandwich covariance estimator and uses code contributed by Roland Gerard Gera.
#'
#' @references Liang Kung-Yee, Zeger Scott L. (1986); Jung Sin-Ho, Ahn Chul (2003); Wachtlin Daniel Kieser Meinhard (2013)
#'
#' @examples
#' #Lets assume we wish to calculate the robust variance estimator for the equation 
#' #\eqn{y_{ij}=\beta_0+\beta_1*I_{treat}+\beta_2*j+\beta_3*I _{treat}*j+\epsilon_{ij}}. 
#' #Furthermore we use the identitiy matrix as working covariance matrix. 
#' #We compare the results with the same estimation made by \code{sandwich2} to show the 
#' #same results. The cance to get randomized to treatment is 60 percent and we observe 
#' #the timerange 0:5.
#'    
#'   ycov = gen_cov_cor(var = 3,rho = 0.25,theta = 1,Time = 0:5,cov = TRUE)
#'   D = matrix(c(1,0.6,0,0,
#'                1,0.6,1,0.6,    
#'                1,0.6,2,1.2,  
#'                1,0.6,3,1.8,    
#'                1,0.6,4,2.4,  
#'                1,0.6,5,3.0),nrow=4)
#'                
#'  D=t(D)            
#'  V=diag(1,length(0:5))
#'  #We correct entries where E(D_i %*% D_i) is unequal to E(D_i)%*%E(D_i) (D %*% D). 
#'  correctionmatrix=matrix(c(1,1,1,1,1,1/0.6,1,1/0.6,1,1,1,1,1,1/0.6,1,1/0.6),nrow=4)
#'  missingtype = "none"
#'  
#'  robust=sandwich(yCov=ycov,D=D,V=V,missingtype=missingtype,correctionmatrix=correctionmatrix)
#'  robust
#'  
#'  #To see if that is correct we can verify it with function sandwich2, which is usable for 
#'  #this particular model with: 
#'  robust2=sandwich2(sigma = c(3,3),rho = 0.25,theta = 1,k = 1.5,Time = 0:5,
#'  dropout = rep(0,6),Model = 1)
#'  robust2[[1]]
#'  
#'  # We can also test this with the the Model: 
#'  #\eqn{y_{ij}=\beta_0+\beta_2*j+\beta_3*I _{treat}*j+\epsilon_{ij}} which leads to 
#'  D = matrix(c(1,0,0,
#'                1,1,0.6,    
#'                1,2,1.2,  
#'                1,3,1.8,    
#'                1,4,2.4,  
#'                1,5,3.0),nrow=3)
#'  D=t(D)            
#'  V=diag(1,length(0:5))
#'  #We correct entries where E(D_i %*% D_i) is unequal to E(D_i)%*%E(D_i) (D %*% D). 
#'  correctionmatrix =matrix(c(1,1,1, 1,1,1, 1,1,1/0.6),nrow=3)
#'  missingtype = "none"
#'  
#'  robust=sandwich(yCov=ycov,D=D,V=V,missingtype=missingtype,correctionmatrix=correctionmatrix)
#'  robust
#'  robust2=sandwich2(sigma = c(3,3),rho = 0.25,theta = 1,k = 1.5,Time = 0:5,
#'  dropout = rep(0,6),Model = 2)
#'  robust2[[1]]
#'  
#' @export

sandwich <- function(yCov , D, V, correctionmatrix,missing=rep(0,dim(yCov)[[2]]), missingtype=c("none","monotone","intermittened") ){
#  require(MASS)
  
  if(missingtype=="none"){
    remain_i = rep(1,dim(yCov)[[2]])
    remain_ij=remain_i %*% t(remain_i) #p_{i,j}= p_i*p_j
  }
  if(missingtype=="intermittened"){
    # calculate percentage of remaining subjects per timepoint
    remain_i = 1-missing
    # calculate probabilety to remain in study at time i AND j
    remain_ij=remain_i %*% t(remain_i) #p_{i,j}= p_i*p_j
    diag(remain_ij)=remain_i # bzw p_{i,i}=\p_i
  } 
  if(missingtype=="monotone"){
    # calculate percentage of remaining subjects per timepoint
    remain_i = 1-missing
    # calculate probabilety to remain in study at time i AND j
    remain_ij=diag(rep(0,dim(yCov)[[1]]))
    for (i in 1:length(missing)){
      for (j in 1:length(missing)){
        # alternativly calculate min(remain_i[i],remain_i[j])]
        remain_ij[i,j]=remain_i[max(i,j)]
      }
    }
  }

  Bread=t(D)%*%diag(remain_i)%*%ginv(V)%*%D*correctionmatrix
  invbread=ginv(Bread)
  
  Meat=(t(D)%*%ginv(V)%*%(remain_ij*yCov)%*%ginv(V)%*%D)*correctionmatrix
  
  robust=invbread%*%Meat%*%invbread
  return(robust)
}