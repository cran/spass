#' @title Generate a dataset of normally distributed repeated measures in a one subgroup setting
#' @description \code{r.gee.1subgroup} generates data of a population which is comprised of a subgroup and the complementary subgroup. The generated longitudinal data needs the specification of the correlation (\eqn{\rho}),the correlation structure (\eqn{\theta}) and the number of repeated measurements. The intra-subject correlation is defined via \eqn{corr(y_{ij},y_{io})=\rho^{(j-o)^\theta}} for the correlation between timepoints \eqn{i} and \eqn{o}.The outcomes are generated as follows:
#' \deqn{y_{ij}=\beta_0+\beta_1*I_{treat}+\beta_2*j+\beta_3*I _{treat}*j+\epsilon_{ij}}
#'  with \eqn{i} being the subject index and \eqn{j} being the time index.  The regression coefficients and outcome-variance for subpopulation and complementary population can be defined seperatly .
#'
#' @param n        overall sample size that is generated
#' @param reg      a list containing regression coefficients for complementary population, \code{reg[[1]]} and subpopulation, \code{reg[[2]]}: see 'Details'
#' @param sigma    vector of standard deviations for \eqn{\epsilon_{ij}}, c(complementary population, subpopulation)
#' @param rho      correlation between two adjacent timepoints 1 timeunit appart
#' @param theta    variable specifying the type of the correlation structure: see 'Details'
#' @param tau      prevalence of the subgroup in the full population.
#' @param k        sample size allocation factor between control and treatment: see 'Details'.
#' @param Time     list with the time-values which are taken by \eqn{j}: see 'Details'
#' @param OD       overall dropout observed at the last timepoint in percent: see 'Details'
#'
#'
#'
#' @details
#'  Given the coefficients \code{reg}=\code{list}(c(\eqn{\beta_0^F\S,\beta_1^F\S,\beta_2^F\S,\beta_3^F\S}), c(\eqn{\beta_0^S,\beta_1^S,\beta_2^S,\beta_3^S}))
#'  and the outcome-variance \code{sigma}=(\eqn{\sigma_F\S, \sigma_S})
#' function \code{r.gee.1subgroup} generates data with intra-subject correlation defined by variables \eqn{\rho} and \eqn{\theta} as follows:
#'
#' Placebo group - complementary population \eqn{y_{ij}=\beta_0+\beta_2*j+N(0,\sigma_F\S)},
#' Placebo group - within subgroup \eqn{y_{ij}=\beta_0+\beta_2*j+N(0,\sigma_S)},       
#' Treatment group - complementary population \eqn{y_{ij}=\beta_0+\beta_1+\beta_2*j+\beta_3*j+N(0,\sigma_F\S)},      
#' Treatment group - within subgroup \eqn{y_{ij}=\beta_0+\beta_1+\beta_2*j+\beta_3*j+N(0,\sigma_S)}.
#' 
#' The intra-subject correlation is included by correlating the error terms \eqn{\epsilon_{ij}}. The formula which describes the correlation between two timepoints is \eqn{corr(\epsilon_ij,\epsilon_io)=\rho^{(j-o)^\theta}}. If for example \eqn{\theta=0} the correlation is compound symmetric. With \eqn{\theta=0} the data is AR(1) correlated.
#' 
#'  Argument \code{k} is the sample size allocation factor, i.e. the ratio between control and treatment. Let \eqn{n_C} and \eqn{n_T} denote the sample sizes of of the control and treatment group, respectively, then \eqn{k = n_T/n_C}.
#' Argument \code{Time} is a vector which are the measurment times, i. e. all the timepoint where a measurement was taken. For \code{Time}=0:5 measurments at baseline, and at timepoints 1,2,3,4 and 5 where taken.
#' 
#' Argument \code{OD} sets the overall dropout rate at the last timepoint. For \code{OD}=0.5 50 percent of all observation had an dropout event. If a subjact has a dropout the chance for that dropout is equally distributet over all time points. 
#'
#'
#' @return \code{r.gee.1subgroup} returns a list with 7 diffrent matrices.In every Matrix the rows are the simulated subjects and the columns are the observed time points.
#' 
#' The first matrix contains the id's of the subject. The id's range from 1 to N.
#' The second are are the outcomes of a subject, \eqn{y_ij}, and so the dependent variable in most analysis. The outcome for \eqn{y_ij} can be found in row i at the corresponding collumn for j. 
#' Matrix 3 to 5 are the values for the independent variables \code{Baseline}, \code{Gr} and \code{Time}. All entries of \code{Baseline} are 1 and as such the baseline of  control pations is defined by \eqn{\beta_0}. The enries of \code{Gr} corresponds to coefficient \eqn{\beta_1}, \code{Time} to coefficient \eqn{\beta_2} and the result of \code{Gr}*\code{Time} to coefficient \eqn{\beta_3}.
#' The sixth matrix contains the \code{error}-terms to preserve the abilety to look tat them later.
#' The last matrix provides the invoramtion if an observation comes from an subjoct of the subpopulation or the complementary population.
#'
#' @source \code{r.gee.1subgroup} uses code contributed by Roland Gerard Gera
#'
#' @examples
#'
#' set.seed(2015)
#' dataset<-r.gee.1subgroup(n=200, reg=list(c(0,0,0,0.1),c(0,0,0,0.1)), sigma=c(3,2.5), 
#' tau=0.5, rho=0.25, theta=1, k=1.5, Time=c(0:5), OD=0)
#' dataset
#' @export


r.gee.1subgroup<-function(n, reg, sigma, rho ,theta , tau, k, Time, OD){
  
  size_S  = round(n*tau)
  size_SC = n-size_S
  
  #  if no data from the subpopulation is generated
  if (size_S==0) {
    
    Datenliste_S = c()
    
  } else {
    
    Datenliste_S <-Datagen(Varianz = sigma[2] , 
                           rho = rho , 
                           theta = theta , 
                           k = k,
                           Koeffizienten = reg[[2]] , 
                           n = size_S , 
                           Time = Time, 
                           OverallDropout = OD, 
                           Typ="S")
    
  }
  # Falls keine Daten aus der komplement?ren Subgruppe generiert werden m?ssen
  if ( size_SC == 0 ) {
    
    Datenliste_SC <-c()
    
  } else {
    
    Datenliste_SC <- Datagen(Varianz = sigma[1] , 
                             rho = rho , 
                             theta = theta , 
                             Koeffizienten = reg[[1]] , 
                             k = k,  
                             n = size_SC , 
                             Time = Time , 
                             OverallDropout = OD,
                             Typ="SC")
    #  Falls Daten aus S erzeugt wurden, incrementiere die ID aus SC, um diese Anzupassen
    if (size_S!=0) Datenliste_SC$id<-Datenliste_SC$id+max(Datenliste_S$id)
  }
  
  
  Datenliste <- list()
  # Fuege die Daten zusammen
  for (i in 1:length(Datenliste_SC)){
    
    Datenliste[[i]] <-rbind(Datenliste_S[[i]] , Datenliste_SC[[i]])
    
  }
  
  # Gebe den Listeneintr?gen die Namen aus SC
  names(Datenliste)=names(Datenliste_S)
  
  return(Datenliste)
}


# Function which generates the List


# Generierung von Simulationsdaten, H1 True und H1 False ------------------
Datagen <- function(Varianz, rho , theta , Koeffizienten , k , n , Time , OverallDropout , Typ){
  #require(MASS)
  kp=1/(k+1)
  TimePoints=length(Time)
  # Erstelle eine "Maske", fuer NA um dropouts zu erzeugen
  Remain.matrix=matrix(1 , nrow = n , ncol = TimePoints)
  
  # only if droppout is != 0
  if(OverallDropout != 0){
    # calculate which patients expereince a dropout
    
    drops=sample(1:n , OverallDropout*n)
    
    # the range in which dropouts can occure
    range=2:length(Time)
    
    for(i in 1:length(drops)){
      
      # draw where at random the missing begins
      k=sample(range,1)
      #calculate matrix where 1 is the part where no dropout occured and NA where dropout exists
      Remain.matrix[drops[i],k:length(Time)]=NA
      
    }
  }
  
  ## Bestimme die Korrelation, die zwischen den Daten bestehen soll -----
  Korrelation <- gen_cov_cor(var=Varianz,
                             rho=rho,
                             theta=theta, 
                             Time=Time,
                             cov=TRUE)
  
  
  # Berechne die einzelnen Parameter der Formel -----------------------------
  id = matrix(1 , ncol=TimePoints , nrow=n) * 1:n
  
  error = mvrnorm(n , mu=matrix(0,ncol=TimePoints) , Sigma=Korrelation) * Remain.matrix
  
  Baseline = matrix(1 , nrow=n , ncol = TimePoints)
  
  Gr = rbind(matrix(0 , nrow = floor(n*kp)   , ncol = TimePoints),
             matrix(1 , nrow = n-floor(n*kp) , ncol = TimePoints))
  
  time = t( matrix(rep(Time,n) , nrow = TimePoints) )
  
  Gr_x_Time = time*Gr
  
  y = Baseline*Koeffizienten[1]+
    Gr*Koeffizienten[2]+
    time*Koeffizienten[3]+
    Gr_x_Time*Koeffizienten[4]+
    error
  
  population = matrix(Typ , ncol = TimePoints , nrow = n)
  
  # Fuege alle einzelne "Messdaten zusammen" ------------
  study = list(id = id ,
               y=y , 
               Baseline=Baseline, 
               Gr=Gr,
               Time=time,
               error=error, 
               population = population)
  class(study)="repdata"
  return(study)
}