
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;


double min(double x, double y){
  if(x<y){
    return(x);
  }else{
    return(y);
  }
}

double uber(double x, double y){
  return(tgamma(x+1)/(tgamma(y+1)*tgamma(x-y+1)));
}

double dnbinom(int k, double mu, double eta){
  return(tgamma(k+eta)/(tgamma(k+1)*tgamma(eta))*pow((mu/(eta+mu)),k)*pow((eta/(eta+mu)),eta));
}

double beta(double p, double q){
  return(tgamma(p)*tgamma(q)/tgamma(p+q));
}

double dnbinomCond(int j, int k, double mu, double size, double rho){
  int y;
  double sumh;

  sumh=0;
  for(y=0;y<(min(j,k)+1);y++){
    sumh += uber(j, y)*beta(rho*size+y, (1-rho)*size+j-y)/beta(rho*size, (1-rho)*size)*pow(size/(size+mu), (1-rho)*size)*tgamma(k+(1-rho)*size-y)/tgamma(k-y+1)/tgamma((1-rho)*size)*pow(mu/(size+mu), k-y);
  }
  return sumh;
}

// [[Rcpp::export]]
double minFunc(NumericVector x, NumericVector daten, int dataNA){

  int t;
  double logL, mu, size, rho;

  logL=0;
  mu = x[0];
  size=x[1];
  rho=x[2];

  logL += log(dnbinom(daten[0], mu, size));
  for(t=0;t<(dataNA-1);t++){
    logL += log(dnbinomCond(daten[t], daten[t+1], mu, size, rho));
  }

  return -logL;
}

// [[Rcpp::export]]
double minFuncMult(NumericVector x, NumericMatrix daten, NumericVector dataNA, int n){

  int j,t;
  double logL, mu, size, rho;

  logL=0;
  mu = x[0];
  size=x[1];
  rho=x[2];

  for(j=0;j<n;j++){
    logL += log(dnbinom(daten(j,0), mu, size));
    for(t=0;t<(dataNA[j]-1);t++){
      logL += log(dnbinomCond(daten(j,t), daten(j,t+1), mu, size, rho));
    }

  }

  return -logL;
}

// [[Rcpp::export]]
double minFuncBlinded(NumericVector x, NumericMatrix daten, NumericVector dataNA, NumericVector n, double delta){
  int j, t, nAll;
  double logL, muC, muE, size, rho, k;

  logL=0;
  nAll = n[0]+n[1];
  k=n[1]/n[0];
  muE=x[0]*(1+k)/(k+1/delta);
  muC=x[0]*(1+k)/(1+k*delta);
  size=x[1];
  rho=x[2];

  for(j=0;j<nAll;j++){
    logL += log(1/(1+k)*(k*dnbinom(daten(j,0), muE, size)+dnbinom(daten(j,0), muC, size)));
    for(t=0;t<(dataNA[j]-1);t++){
      logL += log(1/(1+k)*(k*dnbinomCond(daten(j,t), daten(j,t+1), muE, size, rho) + dnbinomCond(daten(j,t), daten(j,t+1), muC, size, rho)));
    }
  }
  return -logL;
}

