#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppEigen.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Eigen::VectorXd workNat(const Eigen::Map<Eigen::VectorXd> &para_tilde,
                   const bool &LEVIER){
  Eigen::VectorXd para=para_tilde;

  para(0)=1/(1+exp(para_tilde(0)));      //omega
  para(1)=1/(1+exp(para_tilde(1)));      //a
  para(2)=1+exp(para_tilde(2));      //b
  para(3)=exp(para_tilde(3));            //sigma
  para(4)=1/(1+exp(para_tilde(4)));    //v_0

  para(5)=para_tilde(5);    //xi
  para(6)=para_tilde(6);    //varphi
  para(7)=para_tilde(7);   //delta_1
  para(8)=para_tilde(8);    //delta_2
  para(9)=exp(para_tilde(9));//shape
  if(LEVIER){
    para(10)=exp(para_tilde(10)); //l
    para(11)=1/(1+exp(para_tilde(11))); //theta
  }
  return para;
}

// [[Rcpp::export]]
Eigen::VectorXd natWork(const Eigen::Map<Eigen::VectorXd> &para,
                              const bool &LEVIER){
  Eigen::VectorXd para_tilde=para;

  para_tilde[0]=log((1/para[0])-1);      //omega

  para_tilde[1]=log((1/para[1])-1);      //a
  para_tilde[2]=log(para[2]-1);      //b

  para_tilde[3]=log(para[3]);            //sigma
  para_tilde[4]=log((1/para[4])-1);    //v_0

  para_tilde[5]=para[5];    //xi
  para_tilde[6]=para[6];    //varphi
  para_tilde[7]=para[7];    //delta_1
  para_tilde[8]=para[8];    //delta_2
  para_tilde[9]=log(para[9]);//shape
  if(LEVIER){
    para_tilde[10]=log(para[10]); //l
    para_tilde[11]=log((1/para[11])-1); //theta
  }
  return para_tilde;
}

// [[Rcpp::export]]
NumericVector proba_pi(const double &omega, const int &N_D){
  NumericVector probaPi=dbinom(seq(0,N_D-1),N_D-1,omega);
  return probaPi;
}

// [[Rcpp::export]]
Eigen::VectorXd volatilityVector(const Eigen::VectorXd &para, const int &N_D, const int &N){
  Eigen::VectorXd sigma=Eigen::VectorXd::Ones(1);
  Eigen::VectorXd sigma_i(N_D);
  Eigen::VectorXd probaPi=as<Eigen::VectorXd>(proba_pi(para[0],N_D));
  double e_i;
  for(int k=0; k<N_D;k++) sigma_i[k]=para[4]*pow(((2-para[4])/para[4]),k);
  e_i=(probaPi.transpose())*sigma_i;
  sigma_i=sigma_i/e_i;
  for(int i(0);i<N;i++ ) sigma=kroneckerProduct(sigma,sigma_i).eval();
  sigma=para[3]*sigma;
  return sigma;
}

// [[Rcpp::export]]
Eigen::VectorXd probapi(const double &omega, const int &N_D, const int &N){
  Eigen::VectorXd probaPi=as<Eigen::VectorXd>(proba_pi(omega,N_D));
  Eigen::VectorXd proba=as<Eigen::VectorXd>(proba_pi(omega,N_D));
  if(N==1) return(probaPi);
  for(int i(1);i<N;i++) probaPi=kroneckerProduct(probaPi,proba).eval();
  return(probaPi);
}

// [[Rcpp::export]]
Eigen::MatrixXd transitionMatrix(const double &phi, const double &omega, const int &N_D){
  return phi*(Eigen::MatrixXd::Identity(N_D,N_D))+(1-phi)*(Eigen::VectorXd::Ones(N_D))*(as<Eigen::VectorXd>(proba_pi(omega,N_D))).transpose();
}

// [[Rcpp::export]]
mat P(const Eigen::VectorXd &para, const int &N_D,const int &N){
  colvec phi(N);
  phi[0]=para[1];
  for(int k=1; k<N; k++) phi[k]=pow(para[1],(pow(para[2],k)));
  mat P=as<mat>(wrap(transitionMatrix(phi[0],para[0],N_D)));
  if(N==1) return(P);
  for(int i(1);i<N;i++) P=kron(P,as<mat>(wrap(transitionMatrix(phi[i],para[0],N_D)))).eval();
  return(P);
}

// [[Rcpp::export]]
Rcpp::List levierVolatility(const NumericVector &ech, const int &Nl, const Eigen::VectorXd &para){
  int t=ech.size();
  NumericVector Levier=wrap(Eigen::VectorXd::Ones(t));
  NumericVector li(Nl);
  double levier;
  for(int i=0;i<Nl;i++) li[i]=para[10]*pow(para[11],i);
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
Eigen::VectorXd r_dens(const double &x, const Eigen::VectorXd &sd){
  int n = sd.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i] = R::dnorm(x, 0, sqrt(sd[i]), FALSE);
  }
  return(res);
}

// [[Rcpp::export]]
Eigen::VectorXd rv_dens(const double &x,const double &shape,const Eigen::VectorXd &mu_rv){
  int n = mu_rv.size();
  Eigen::VectorXd res(n);
  for(int i = 0; i < n; i++) {
    res[i]=R::dlnorm(x,mu_rv[i],sqrt(shape),FALSE);
  }
  return(res);
}

// [[Rcpp::export]]
double logLik(const Eigen::MatrixXd &ech,
                    const Eigen::Map<Eigen::VectorXd> &para_tilde,
                    const bool &LEVIER,
                    const int &N_D,
                    const int &N,
                    const int &Nl=70){
  Eigen::VectorXd para=workNat(para_tilde,LEVIER);
  int n=ech.rows();
  double xi=para[5];
  double varphi=para[6];
  double delta1=para[7];
  double delta2=para[8];
  double shape=para[9];
  Eigen::VectorXd sigma = volatilityVector(para,N_D,N);
  Eigen::VectorXd aj;
  Eigen::VectorXd p0=probapi(para[0],N_D,N);
  double a(0);
  double lik(0);
  Eigen::VectorXd w;
  double x=ech(0,1);
  double y=ech(0,0);
  Eigen::VectorXd temp=x*(sigma.array().sqrt().inverse());
  Eigen::VectorXd mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
  Eigen::MatrixXd matP = as<Eigen::MatrixXd>(wrap(P(para,N_D,N)));
  Eigen::VectorXd Res;
  aj=r_dens(x,sigma).array()*rv_dens(y,shape,mu_rv).array();
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;
  
  if(!LEVIER){
    for(int i(1); i<n;i++){
      x=ech(i,1);
      y=ech(i,0);
      temp=x*(sigma.array().sqrt().inverse());
      mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
      Res=r_dens(x,sigma).array()*rv_dens(y,shape,mu_rv).array();
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }
  if(LEVIER){
    NumericVector ech1=Rcpp::wrap(ech.col(1));
    Rcpp::List L=levierVolatility(ech1,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    for(int i(1); i<n;i++){
      x=ech(i,1);
      y=ech(i,0);
      sigma1=sigma*Levier[i];
      temp=x*(sigma1.array().sqrt().inverse());
      mu_rv=xi+varphi*(sigma1.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
      Res=r_dens(x,sigma1).array()*rv_dens(y,shape,mu_rv).array();
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;
    }
  }
  return (-lik);
}

// [[Rcpp::export]]
Rcpp::List logLik2(const Eigen::MatrixXd &ech,
                   const Eigen::Map<Eigen::VectorXd> &para_tilde,
                   const bool &LEVIER,
                   const int &N_D,
                   const int &N,
                   const double &r=0,
                   const int &t=2,
                   const int &Nl=70){
  Eigen::VectorXd para=workNat(para_tilde,LEVIER);
  int n=ech.rows();
  double xi=para[5];
  double varphi=para[6];
  double delta1=para[7];
  double delta2=para[8];
  double shape=para[9];
  Eigen::VectorXd sigma = volatilityVector(para,N_D,N);
  Eigen::VectorXd aj;
  Eigen::VectorXd ajm;
  Eigen::VectorXd p0=probapi(para[0],N_D,N);
  double a(0);
  double lik(0);
  double likm(0);
  double pred_lik(0);
  Eigen::VectorXd w;
  Eigen::VectorXd w_hat;
  double x=ech(0,1);
  double y=ech(0,0);
  Eigen::VectorXd temp=x*(sigma.array().sqrt().inverse());
  Eigen::VectorXd mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
  Eigen::MatrixXd matP = as<Eigen::MatrixXd>(wrap(P(para,N_D,N)));
  Eigen::VectorXd Res;
  Eigen::VectorXd Resm;
  aj=r_dens(x,sigma).array()*rv_dens(y,shape,mu_rv).array();
  aj=p0.array()*aj.array();
  a=aj.sum();
  lik=log(a);
  w = aj/a;

  ajm=r_dens(x,sigma);
  ajm=p0.array()*ajm.array();
  likm=log(ajm.sum());

  if(!LEVIER){
    for(int i(1); i<n;i++){
      x=ech(i,1);
      y=ech(i,0);
      temp=x*(sigma.array().sqrt().inverse());
      mu_rv=xi+varphi*(sigma.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
      Res=r_dens(x,sigma).array()*rv_dens(y,shape,mu_rv).array();
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;

      if(i==(t-1)) w_hat=w;

      Resm=r_dens(x,sigma);
      ajm=(((w.transpose())*matP).array())*(Resm.transpose().array());
      likm+=log(ajm.sum());
    }
    Res=r_dens(r,sigma);
    aj=(((w.transpose())*matP).array())*(Res.transpose().array());
    a=aj.sum();
    pred_lik=log(a);
  }
  if(LEVIER){
    NumericVector ech1=Rcpp::wrap(ech.col(1));
    Rcpp::List L=levierVolatility(ech1,Nl,para);
    NumericVector Levier=L["Levier"];
    Eigen::VectorXd sigma1;
    for(int i(1); i<n;i++){
      x=ech(i,1);
      y=ech(i,0);
      sigma1=sigma*Levier[i];
      temp=x*(sigma1.array().sqrt().inverse());
      mu_rv=xi+varphi*(sigma1.array().log())+delta1*temp.array()+delta2*(temp.array().pow(2)-1);
      Res=r_dens(x,sigma1).array()*rv_dens(y,shape,mu_rv).array();
      aj=(((w.transpose())*matP).array())*(Res.transpose().array());
      a=aj.sum();
      lik+=log(a);
      w=aj/a;

      if(i==(t-1)) w_hat=w;

      Resm=r_dens(x,sigma1);
      ajm=(((w.transpose())*matP).array())*(Resm.transpose().array());
      likm+=log(ajm.sum());
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
      Rcpp::Named("Marg_loglik") = likm,
      Rcpp::Named("Pred_loglik") = pred_lik,
      Rcpp::Named("w_hat") = w_hat
    );
  return output;
}

// [[Rcpp::export]]
NumericMatrix levierVolatilityMat(const NumericMatrix &ech, const NumericMatrix Levier, const int &Nl, const NumericVector &para){
  NumericVector li(Nl);
  colvec levier=as<colvec>(wrap(Eigen::VectorXd::Ones(ech.rows())));
  for(int i=0;i<Nl;i++) li[i]=para[10]*pow(para[11],i);
  NumericVector y1;
  NumericVector y2;
  LogicalVector x;
  int t=ech.cols();
  for(int i=0; i<Nl;i++) {
    y1 = -ech(_,t-i-1);
    x=(y1>0);
    y2=1/sqrt(Levier(_,t-i-1));
    levier=levier%(1+(li[i]*as<colvec>(y1)%as<colvec>(x)%as<colvec>(y2)));
  }
  NumericVector l=wrap(levier);
  NumericMatrix Levier2=cbind(Levier,l);
  return Levier2;
}

// [[Rcpp::export]]
arma::rowvec colMeansRcpp(NumericMatrix x){
  arma::mat X = arma::mat(x.begin(), x.nrow(), x.ncol(), false); 
  return arma::mean(X, 0); 
}

// [[Rcpp::export]]
NumericMatrix Pow(const NumericMatrix &x,const double &n){
  Eigen::MatrixXd y = as<Eigen::MatrixXd>(x); 
  return wrap(y.array().pow(n)); 
}

// [[Rcpp::export]]
Rcpp::List R_hat(const int &H,
                    const Eigen::Map<Eigen::VectorXd> &ech,
                    const Eigen::Map<Eigen::MatrixXi> &MC_sim,
                    const Eigen::Map<Eigen::MatrixXd> &z_t,
                    NumericMatrix Levier,
                    const NumericVector &sig,
                    const NumericVector &para,
                    const int &N=3,
                    const int &Nl=70){
  NumericVector temp;
  NumericVector temp1;
  NumericVector temp2;
  NumericMatrix rt_sim(z_t.rows(),ech.size());
  rt_sim=wrap(kroneckerProduct(Eigen::VectorXd::Ones(z_t.rows()),ech.transpose()));
  for(int h=0; h<H;h++){
    Levier=levierVolatilityMat(rt_sim,Levier,Nl,para);
    IntegerVector idx=wrap(MC_sim.col(h).array()-1);
    temp1=sig[idx];
    temp2=wrap(z_t.col(h));
    temp=sqrt(Levier(_,Levier.ncol()-1)*temp1)*temp2;
    rt_sim=cbind(rt_sim,temp);
  }
  NumericMatrix rt2 = Pow(rt_sim,2);
  NumericMatrix rvt=Pow(rt2,para[6])*exp(para[5]+0.5*para[9]+0.5*pow(para[7],2)-para[8])/sqrt(1-2*para[8]);
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("rt_sim") = rt_sim,
      Rcpp::Named("rt2") = colMeansRcpp(as<NumericMatrix>(rt2)),
      Rcpp::Named("rvt") = colMeansRcpp(as<NumericMatrix>(rvt))
    );
  return output;
}

// [[Rcpp::export]]
Rcpp::List f_sim(const int &H,
                      const Eigen::Map<Eigen::VectorXd> &sig,
                      const Eigen::Map<Eigen::VectorXd> &pi_0,
                      const Eigen::Map<Eigen::MatrixXd> &matP,
                      const double &varphi,
                      const double &xi,
                      const double &shape,
                      const double &delta1,
                      const double &delta2){
  Eigen::MatrixXd temp = matP;
  NumericVector temp2_r(H);
  NumericVector temp2_rv(H);
  NumericVector temp3_r =  wrap((pi_0.transpose())*sig);
  Eigen::VectorXd sig2=sig.array().pow(varphi);
  NumericVector temp3_rv =  wrap((pi_0.transpose())*sig2);
  temp2_r[0] = temp3_r[0];
  temp2_rv[0] = temp3_rv[0]*(exp(xi+0.5*shape+0.5*pow(delta1,2)-delta2)/sqrt(1-2*delta2));
  for(int h=1;h<H;h++){
    temp3_r = wrap((pi_0.transpose() * temp)*sig);
    temp3_rv = wrap((pi_0.transpose() * temp)*sig2);
    temp2_r[h] = temp3_r[0];
    temp2_rv[h] = temp3_rv[0]*(exp(xi+0.5*shape+0.5*pow(delta1,2)-delta2)/sqrt(1-2*delta2));
    temp = temp * matP;
  }
  Rcpp::List output =
    Rcpp::List::create(
      Rcpp::Named("rt2") = temp2_r,
      Rcpp::Named("rvt") = temp2_rv
    );
  return output;
}



/*** R

# para<-c(0.67044,  0.99481,  5.58986,  2.33798,  0.80352, -0.31636,  0.90798, -0.03042,  0.07726,  0.14792)
# vars<-c("omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape")
# ech<-c(0.7,0.5,-0.3,-1.75,-0.7,0.5,-0.5,0.8,0.3,-0.7,0.8)
# ech<-cbind(rv=c(0.15,1.5,0.7,0.8,1.9,0.8,0.3,0.54,0.76,0.98,0.01),r=ech)
# ech1<-donnees[[1]]
# ech<-ech1[1:5012,]
# N_D<-10
# N<-3
# LEVIER<-FALSE
# if(LEVIER) {
#   para<-c(para,1.75,0.89)
#   vars<-c(vars,"l","theta")
# }
# names(para)<-vars
# 
# para_tilde<-natWork(para,LEVIER)
# 
# logLik(ech,para_tilde,LEVIER,N_D,N,Nl)
# 
# LEVIER<-TRUE
# Nl<-70
# if(LEVIER) {
#   para<-c(para,1.75,0.89)
#   vars<-c(vars,"l","theta")
# }
# names(para)<-vars
# 
# levierVolatility(ech,Nl,para)
# para_tilde<-natWork(para,LEVIER)
# workNat(para_tilde,LEVIER)
# 
# logLik(ech,para_tilde,LEVIER,N_D,N,Nl=Nl)
# 
# n<-100
# H <- 100    #nb of years simulated
# 
# Levier<-rep(1,n)%*%t(levierVolatility(ech[,2],Nl,para)$`Levier`)
# echMat<-rep(1,n)%*%t(ech[,2])
# 
# l<-logLik2(ech,para_tilde,LEVIER,N_D,N,Nl=Nl,r=0)
# 
# pi_0 <- l$w_hat
# sig <- volatilityVector(para,N_D,N)
# matP<-P(para,N_D,N)
# 
# #simulation
# set.seed(1984)
# 
# library(mhsmm)
# MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(H, n)),H,n,byrow=FALSE)) #simulation of Markov chain
# z_t<-matrix(rnorm(n*H),nrow=n,ncol=H)
# r_t2<-R_hat(H,ech[,2],MC_sim,z_t,Levier,sig,para,N,Nl)
# rvt_sim <- r_t2$rvt
# rt2_sim  <- r_t2$rt2


*/
