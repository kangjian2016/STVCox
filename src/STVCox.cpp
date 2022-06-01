#include <RcppArmadillo.h>
#include <math.h>

// alph: was lambda

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

#define pi 3.14159265


//'@importFrom Rcpp evalCpp
//'@useDynLib STVCox, .registration=TRUE


// [[Rcpp::export]]
double h_beta_t0_Rcpp(double tau, double beta_t, double alph){
  double a1=beta_t-alph;
  double b1=a1/tau;
  double a2=beta_t+alph;
  double b2=a2/tau;
  double bt=0.5*(1+2/pi*atan(b1))*a1+0.5*(1-2/pi*atan(b2))*a2;
  return bt;
}

//'@export
// [[Rcpp::export]]
double lplk_Rcpp(arma::vec par, arma::vec alph, arma::vec time, 
                 arma::vec delta, arma::mat z, arma::mat Bs){
  int n = time.n_elem;
  int p = z.n_cols;
  int qn = Bs.n_cols;
  double tau=0.01;
  double llk=0;
  arma::vec betai(p);
  arma::vec thetaj(qn);
  for(int i = 0; i < n; i++){
    if(delta(i)==1){
      for(int j = 0; j<p; j++){
        int first = j*qn;
        int last = first + qn-1;
        thetaj = par.subvec(first,last);
        double beta_t = arma::accu(Bs.row(i).t()%thetaj);
        betai(j) = h_beta_t0_Rcpp(tau, beta_t, alph(j));
        
      }
      arma::vec zbeta(n-i);
      zbeta=z.rows(i,n-1)*betai;
      double zbeta_max = zbeta.max();
      double esum=0;
      for(int k=0;k<zbeta.n_elem;k++){
        esum+=exp(zbeta(k)-zbeta_max);
      }
      llk = llk+arma::accu(z.row(i).t()%betai)-log(esum)-zbeta_max;
    }//if(delta(i)==1)
  }//for i loop
  return -llk;
}


//'@export
// [[Rcpp::export]]
double lplk_Rcpp_penalty(arma::vec par, arma::vec alph, arma::vec time, 
                         arma::vec delta, arma::mat z, arma::mat Bs, double rho){
  int n = time.n_elem;
  int p = z.n_cols;
  int qn = Bs.n_cols;
  double tau=0.01;
  double llk=0;
  double penalty = 0;
  arma::vec betai(p);
  arma::vec thetaj(qn);
  for(int j = 0; j < p; j++){
    int first = j*qn;
    int last = first + qn-1;
    thetaj = par.subvec(first,last);
    penalty = penalty + mean(thetaj%thetaj);
  }
  
  for(int i = 0; i < n; i++){
    if(delta(i)==1){
      for(int j = 0; j<p; j++){
        int first = j*qn;
        int last = first + qn-1;
        thetaj = par.subvec(first,last);
        double beta_t = arma::accu(Bs.row(i).t()%thetaj);
        betai(j) = h_beta_t0_Rcpp(tau, beta_t, alph(j));
        
      }
      arma::vec zbeta(n-i);
      zbeta=z.rows(i,n-1)*betai;
      double zbeta_max = zbeta.max();
      double esum=0;
      for(int k=0;k<zbeta.n_elem;k++){
        esum+=exp(zbeta(k)-zbeta_max);
      }
      llk = llk+arma::accu(z.row(i).t()%betai)-log(esum)-zbeta_max;
    }//if(delta(i)==1)
  }//for i loop
  llk = llk - rho*penalty;
  return -llk;
}


// [[Rcpp::export]]
double lplk_Rcpp_min_beta(arma::vec par, arma::vec alph, arma::vec time, 
                          arma::vec delta, arma::mat z, arma::mat Bs){
  int n = time.n_elem;
  int p = z.n_cols;
  int qn = Bs.n_cols;
  double tau=0.01;
  double llk=0;
  arma::vec betai(p);
  arma::vec thetaj(qn);
  double min_beta = 0;
  for(int i = 0; i < n; i++){
    if(delta(i)==1){
      for(int j = 0; j<p; j++){
        int first = j*qn;
        int last = first + qn-1;
        thetaj = par.subvec(first,last);
        double beta_t = arma::accu(Bs.row(i).t()%thetaj);
        betai(j) = h_beta_t0_Rcpp(tau, beta_t, alph(j));
        if(abs(beta_t) < alph(j)){
          min_beta += abs(beta_t);
        }
      }
      arma::vec zbeta(n-i);
      zbeta=z.rows(i,n-1)*betai;
      double esum=0;
      for(int k=0;k<zbeta.n_elem;k++){
        esum+=exp(zbeta(k));
      }
      llk = llk+arma::accu(z.row(i).t()%betai)-log(esum);
    }//if(delta(i)==1)
  }//for i loop
  double total = -llk + min_beta;
  return total;
}


