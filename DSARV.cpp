#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Eigen::VectorXd r_dens(const double &x, const Eigen::VectorXd &sd){
  int n = sd.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i] = R::dnorm(x, 0, sqrt(sd[i]), FALSE);
  }
  return(res);
}

// [[Rcpp::export]]
Eigen::VectorXd workNat(const Eigen::VectorXd &para_tilde,
                        const bool &LEVIER){
  Eigen::VectorXd para=para_tilde;
  
  para[0]=1/(1+exp(para_tilde[0]));
  para[1]=1/(1+exp(para_tilde[1]));
  para[2]=para_tilde[2];
  para[3]=exp(para_tilde[3]);
  if(LEVIER){
    para(4)=exp(para_tilde(4)); //l
    para(5)=1/(1+exp(para_tilde(5))); //theta
  }
  return para;
}

// [[Rcpp::export]]
Eigen::VectorXd natWork(const Eigen::VectorXd &para,
                        const bool &LEVIER){
  Eigen::VectorXd para_tilde=para;
  
  para_tilde[0]=log((1/para[0])-1);
  para_tilde[1]=log((1/para[1])-1);
  para_tilde[2]=para[2];
  para_tilde[3]=log(para[3]);
  if(LEVIER){
    para_tilde[4]=log(para[4]); //l
    para_tilde[5]=log((1/para[5])-1); //theta
  }
  return para_tilde;
}

// [[Rcpp::export]]
Eigen::VectorXd volatilityVector(const Eigen::VectorXd &para, const int &N){
  Eigen::VectorXd sigma(N);
  for(int i(0);i<N;i++) sigma[i]=exp(para[2]+((i+1)*para[3]));
  return sigma;
}

// [[Rcpp::export]]
Eigen::VectorXd probapi(const double &omega, const int &N){
  NumericVector probaPi=dbinom(seq(0,N-1),N-1,omega);
  return as<Eigen::VectorXd>(probaPi);
}

// [[Rcpp::export]]
Eigen::MatrixXd P(const double &phi, const double &omega, const int &N){
  return phi*(Eigen::MatrixXd::Identity(N,N))+(1-phi)*(Eigen::VectorXd::Ones(N))*(probapi(omega,N)).transpose();
}

// [[Rcpp::export]]
Rcpp::List levierVolatility(const NumericVector &ech, const int &Nl, const Eigen::VectorXd &para){
  int t=ech.size();
  NumericVector Levier=wrap(Eigen::VectorXd::Ones(t));
  NumericVector li(Nl);
  double levier;
  for(int i=0;i<Nl;i++) li[i]=para[4]*pow(para[5],i);
  for(int t=Nl;t<ech.size();t++){
    levier=1;
    for(int i=0; i<Nl;i++) levier=levier*(1+(li[i]*(-(ech[t-i-1]))*(ech[t-i-1]<0)/(sqrt(Levier[t-i-1]))));
    Levier[t]=levier;
  }
  levier=1;
  for(int i=0; i<Nl;i++) levier=levier*(1+(li[i]*(-(ech[t-i-1]))*(ech[t-i-1]<0)/(sqrt(Levier[t-i-1]))));
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("Levier") = Levier,
      Rcpp::Named("levier") = levier
    );
  return output;
} 

// [[Rcpp::export]]
double logLik(const NumericVector &ech,
              const Eigen::Map<Eigen::VectorXd> &para_tilde,
              const bool &LEVIER,
              const int &N,
              const int &Nl=70){
  Eigen::VectorXd para=workNat(para_tilde,LEVIER);
  int n(ech.size());
  Eigen::VectorXd sigma = volatilityVector(para,N);
  Eigen::VectorXd aj;
  Eigen::VectorXd p0=probapi(para[0],N);
  double a(0);
  double lik(0);
  Eigen::VectorXd w;
  Eigen::MatrixXd matP = P(para[1],para[0],N);
  Eigen::VectorXd Res;
  aj=r_dens(ech(0),sigma);
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;
  
  if(!LEVIER){
    for(int i(1); i<n;i++){
      Res=r_dens(ech(i),sigma);
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }
  
  if(LEVIER){
    Rcpp::List L=levierVolatility(ech,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    for(int i(1); i<n;i++){
      sigma1=sigma*Levier[i];
      Res=r_dens(ech(i),sigma1);
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }
  return (-lik);
}

// [[Rcpp::export]]
Rcpp::List logLik2(const NumericVector &ech,
              const Eigen::Map<Eigen::VectorXd> &para_tilde,
              const bool &LEVIER,
              const int &N,
              const double &r=0,
              const int &Nl=70){
  Eigen::VectorXd para=workNat(para_tilde,LEVIER);
  int n(ech.size());
  Eigen::VectorXd sigma = volatilityVector(para,N);
  Eigen::VectorXd aj;
  Eigen::VectorXd p0=probapi(para[0],N);
  double a(0);
  double lik(0);
  double pred_lik(0);
  Eigen::VectorXd w;
  Eigen::MatrixXd matP = P(para[1],para[0],N);
  Eigen::VectorXd Res;
  aj=r_dens(ech(0),sigma);
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;
  
  if(!LEVIER){
    for(int i(1); i<n;i++){
      Res=r_dens(ech(i),sigma);
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
    Res=r_dens(r,sigma);
    aj=(((w.transpose())*matP).array())*(Res.transpose().array());
    a=aj.sum();
    pred_lik=log(a);
  }
  
  if(LEVIER){
    Rcpp::List L=levierVolatility(ech,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    for(int i(1); i<n;i++){
      sigma1=sigma*Levier[i];
      Res=r_dens(ech(i),sigma1);
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
    double levier=L["levier"];
    sigma1=sigma*levier;
    Res=r_dens(r,sigma1);
    aj=(((w.transpose())*matP).array())*(Res.transpose().array());
    a=aj.sum();
    pred_lik=log(a);
  }
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("loglik") = -lik,
      Rcpp::Named("Pred_loglik") = pred_lik
    );
  return output;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

# # DSARV

# N <- 3
# logR<-c(1.2,1.5,-2.3,-0.5,0.7,-0.8,0.8,0.8,0.7,-0.5,-1.3)
# p0<-rep((1/N^2),N^2)
# para<-c(0.182,0.943,-0.269,0.982) # Essai DSARV simple loglinear
# para_tilde<-natWork(para,FALSE)
# para_tilde
# 
# workNat(para_tilde,FALSE)
# 
# volatilityVector(para,N)
# 
# probapi(para[1],N)
# 
# P(para[2],para[1],N)
# 
# logLik(ech=donne,para_tilde=natWork(para,FALSE),FALSE,N)
# s<-nlm(logLik,ech=donne,para_tilde,N=N,LEVIER=FALSE)
# 
# para<-c(workNat(s$estimate,FALSE),1,0.98) # Essai DSARV simple loglinear
# para_tilde<-natWork(para,TRUE)
# 
# levierVolatility(logR,Nl=3,para)
# logLik(ech=logR,para_tilde=para_tilde,TRUE,N,Nl=3)
# 
# nlm(logLik,ech=donne,para_tilde,N=N,LEVIER=TRUE,Nl=70)

#microbenchmark(loglik(logR,para0_tildeDSARV,"normal"),logLikDSARV(logR,para0_tildeDSARV,p0,"normal","linear",NDSARV),times=100)

***/