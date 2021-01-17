#============================================================================
#           JMDSARV (N_D=N, N=nDSARV)
#============================================================================


#-------------------------------------------------------------------------------
#Data input
#-------------------------------------------------------------------------------
if(!require(matrixcalc)){install.packages("matrixcalc")}; library(matrixcalc)
if(!require(Rcpp)){install.packages("Rcpp")}; library(Rcpp)
if(!require(RcppArmadillo)){install.packages("RcppArmadillo")}; library(RcppArmadillo)
if(!require(RcppEigen)){install.packages("RcppEigen")}; library(RcppEigen)

data_source <- "OxfordMan" #"OxfordMan" or "realized_library" or "11series"
# ".AEX" ".AORD" ".BFX" ".BSESN" ".BVLG" ".BVSP" ".DJI" ".FCHI" ".FTMIB" ".FTSE" ".GDAXI"
# ".GSPTSE" ".HSI" ".IBEX" ".IXIC" ".KS11" ".KSE" ".MXX" ".N225" ".NSEI" ".OMXC20" ".OMXHPI"
# ".OMXSPI" ".OSEAX" ".RUT" ".SMSI" ".SPX" ".SSEC" ".SSMI" ".STI" ".STOXX50E"

index_set<-c(".AEX", ".AORD", ".BFX", ".BSESN", ".BVLG", ".BVSP", ".DJI", ".FCHI", ".FTMIB", ".FTSE", ".GDAXI",
             ".GSPTSE", ".HSI", ".IBEX", ".IXIC", ".KS11", ".KSE", ".MXX", ".N225", ".NSEI", ".OMXC20", ".OMXHPI",
             ".OMXSPI", ".OSEAX", ".RUT", ".SMSI", ".SPX", ".SSEC", ".SSMI", ".STI", ".STOXX50E")

# model_extern<-rbind(expand.grid(model=c("MDSV"),r_dis=c("normal"),rv_dis=c("lognormal"),N_D=c(10),N=c(3),
#                          LEVIER=c(FALSE,TRUE),index=c(".FTSE")),
#                     expand.grid(model=c("MDSV"),r_dis=c("normal"),rv_dis=c("lognormal"),N_D=c(5),N=c(4),
#                                 LEVIER=c(FALSE,TRUE),index=c(".FTSE")),
#                     expand.grid(model=c("MDSV"),r_dis=c("normal"),rv_dis=c("lognormal"),N_D=c(2),N=c(10),
#                                 LEVIER=c(FALSE,TRUE),index=c(".FTSE")))

model_extern<-rbind(expand.grid(model=c("MDSV"),r_dis=c("normal"),rv_dis=c("lognormal"),N_D=c(2),N=c(10),
                                LEVIER=c(FALSE,TRUE),index=c(".IXIC")))

para_names<-function(LEVIER){
  vars.names<-c("omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape")
  if(LEVIER) vars.names<-c(vars.names,"l","theta")
  return(vars.names)
}

ctrl <- list(TOL=1e-15, trace=0)

use.kernel <- TRUE #realized kernel instead of realized variance

start.date <- as.Date("2000-01-03") 
end.date   <- as.Date("2019-08-31")
#end.date   <- as.Date("2015-09-23")
#weekdays(start.date); weekdays(end.date);

DMS <- TRUE
if(DMS){
  path.R    <- "/home/maoudek/Rsim/Article1/"
  path.Data <- "/home/maoudek/Rsim/Article1/"
  path.Work <- "/home/maoudek/Rsim/Article1/MDSV/Forecast/"
} else {
  path.R    <- "C:/Users/abdou/Dropbox/Abdoul/Realized volatility/Code"
  path.Work <- "C:/Users/abdou/Dropbox/Abdoul/These/Article1/Code/MDSV"
  path.Data <- "C:/Users/abdou/Dropbox/Abdoul/These/Article1/Data"
}

# rel_tol  <- 1e-6 #optimization parameter

setwd(path.R)
source("realized_functions.r") #load functions for realized volatility

donnees<-list() 
setwd(path.Data)
for(i in 1:length(index_set)){
  index<-index_set[i]
  source("realized_data_traitement.r") #load and organize realized volatility data
  
  RV<-cbind(rv=100*100*Data$RV$rv,r=100*Data$RV$r)
  names(RV)<-as.Date(rownames(Data$RV))
  
  assign(index_set[i],RV)
  donnees[[i]]<-get(index_set[i])
}

#=======================================================

setwd(path.Work)
sourceCpp('MDSVBinomial.cpp')

t0<-9*252 #(9 ans)
H_range<-c(1,5,10,25,50,75,100)

for(i in 1:nrow(model_extern)){
  index<-as.character(model_extern[i,"index"])
  donne<-donnees[[which(index_set==index)]]
  donne[,"r"]<-donne[,"r"]-mean(donne[,"r"])
  
  date_debut<-as.Date(names(donne)[length(donne[,"r"])-t0])
  date_fin<-as.Date(names(donne)[length(donne[,"r"])])
  
  r_dis<-as.character(model_extern[i,"r_dis"])
  rv_dis<-as.character(model_extern[i,"rv_dis"])
  N_D<-as.numeric(model_extern[i,"N_D"])
  N<-as.numeric(model_extern[i,"N"])
  LEVIER<-as.logical(model_extern[i,"LEVIER"])
  
  filename <- paste("MDSVBinomial_N_",N,"_N_D_",N_D,"_", index, "_LEVIER_",LEVIER,"_",r_dis,"_",rv_dis,"_", date_debut, "_", date_fin, sep="")
  
  model<-expand.grid(date = names(donne)[(length(donne[,"r"])-t0):(length(donne[,"r"]))],r_dis=r_dis,rv_dis=rv_dis,N_D=N_D,N=N,Levier=LEVIER)
  
  vars<-c(paste0("var_R_for",H_range),paste0("var_RV_for",H_range),paste0("var_R_tru",H_range),paste0("var_RV_tru",H_range),
          'predict_loglik','Marg_loglik','loglik', 'AIC', 'BIC',
          c("omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape","l","theta"))
  
  model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
  colnames(model_add) <- vars
  model <- cbind(model, model_add)
  
  para<-c(0.59,0.95,1.5,2,0.72,0.2,0.4,-0.1,0.03,3)
  vars<-c("omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape")
  if(LEVIER) {
    para<-c(para,0.7,0.8)
    vars<-c(vars,"l","theta")
  }
  
  # para<-c(0.67044,  0.99481,  5.58986,  2.33798,  0.80352, -0.31636,  0.90798, -0.03042,  0.07726,  0.14792)
  # names(para)<-vars
  para_tilde<-natWork(para,LEVIER)
  # logLik(ech,para0_tilde,N,LEVIER,nDSARV)
  for_R_var <- rep(0,length(H_range))
  for_RV_var <- rep(0,length(H_range))
  RV_var <- rep(0,length(H_range))
  R_var <- rep(0,length(H_range))
  opt<-NULL
  
  for(t in 0:(t0-1)){
    ech<-donne[1:(length(donne[,"r"])-t0+t),]
    ech_rv<-donne[(length(donne[,"r"])-t0+t+1):length(donne[,"rv"]),"rv"]
    ech_r<-donne[(length(donne[,"r"])-t0+t+1):length(donne[,"rv"]),"r"]
    
    update_date<-seq(0,t0,length.out = 3*12+1)
    
    if(t %in% update_date){
      if(!is.null(opt)) para_tilde<-opt$pars
      opt<-try(solnp(pars=para_tilde,fun=logLik,ech=ech,N_D=N_D,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
      if(is(opt,"try-error")) {
        next
      }
      
      para<-workNat(opt$pars,LEVIER)
      names(para)<-para_names(LEVIER)
    }
    
    model[t+1,vars] <- round(para[vars],5)
    l<-logLik2(ech,opt$pars,LEVIER,N_D,N,t=nrow(ech),r=ech_r[1])
    model[t+1,"loglik"]<--as.numeric(opt$values[length(opt$values)])
    model[t+1,"Marg_loglik"]<-l$Marg_loglik
    model[t+1,"predict_loglik"]<-l$Pred_loglik
    model[t+1,'AIC'] <- model[(t+1),"loglik"]-length(para)
    model[t+1,'BIC'] <- model[(t+1),"loglik"]-length(para)*log(length(ech[,"r"]))/2
    
    #simulation
    set.seed(1984)
    
    n <- 10000 #number of simulations
    H <- max(H_range)    #nb of years simulated
    
    pi_0 <- l$w_hat
    sig <- volatilityVector(para,N_D,N)
    matP<-P(para,N_D,N)
    
    MC_sim <- t(matrix(sim.mc(pi_0, matP, rep(H, n)),H,n,byrow=FALSE)) #simulation of Markov chain
    
    w<-para["omega"]
    a<-para["a"]
    b<-para["b"]
    sigma<-para["sigma"] 
    v0<-para["v0"]
    xi<-para["xi"]
    varphi<-para["varphi"]  
    delta1<-para["delta1"]   
    delta2<-para["delta2"]
    shape<-para["shape"]
    l<-para["l"]  
    theta_l<-para["theta"]
    Nl<-70
    
    z_t<-matrix(rnorm(n*H),nrow=n,ncol=H)
    
    if(LEVIER){
      Levier<-rep(1,n)%*%t(levierVolatility(ech[,2],Nl,para)$`Levier`)
      sim<-R_hat(H,ech[,2],MC_sim,z_t,Levier,sig,para,N,Nl)
      rt2_sim<-sim$`rt2`
      rvt_sim<-sim$`rvt`
      rt2_sim <- rt2_sim[,(ncol(rt2_sim)-H+1):ncol(rt2_sim)]
      rvt_sim  <- rvt_sim[,(ncol(rvt_sim)-H+1):ncol(rvt_sim)]
    } else {
      sim<-f_sim(H,sig,pi_0,matP,varphi,xi,shape,delta1,delta2)
      rt2_sim<-sim$`rt2`
      rvt_sim<-sim$`rvt`
      # for(j in 1:H){
      #    rt2_sim_essai[j]<-sum(sig*(pi_0%*%matrix.power(matP, j-1)))
      #    rvt_sim_essai[j]<-sum((sig^varphi)*(pi_0%*%matrix.power(matP, j-1)))*exp(xi+0.5*shape+0.5*delta1^2-delta2)/sqrt(1-2*delta2)
      # }
      
      #rt2_sim<-colMeans(matrix(sig[(MC_sim)],n,H)) (Monte Carlo)
      #rvt_sim<-colMeans(matrix((sig[(MC_sim)])^varphi,n,H))*exp(xi+0.5*shape+0.5*delta1^2-delta2)/sqrt(1-2*delta2) (Monte Carlo)
    }
    
    for(k in 1:length(H_range)) for_RV_var[k] <- sum(rvt_sim[1:H_range[k]])
    names(for_RV_var) <- paste0("var_RV_for",H_range)
    model[(t+1),colnames(model) %in% names(for_RV_var)]<-for_RV_var
    
    for(k in 1:length(H_range)) for_R_var[k] <- sum(rt2_sim[1:H_range[k]])
    names(for_R_var) <- paste0("var_R_for",H_range)
    model[(t+1),colnames(model) %in% names(for_R_var)]<-for_R_var
    
    for(k in 1:length(H_range)) RV_var[k] <- sum(ech_rv[1:H_range[k]])
    names(RV_var) <- paste0("var_RV_tru",H_range)
    model[(t+1),colnames(model) %in% names(RV_var)]<-RV_var
    
    for(k in 1:length(H_range)) R_var[k] <- sum((ech_r[1:H_range[k]])^2)
    names(R_var) <- paste0("var_R_tru",H_range)
    model[(t+1),colnames(model) %in% names(R_var)]<-R_var
    
    write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
    print(paste("===",round(100*(t+1)/nrow(model),2) , "%" ,"====", "===",round(100*i/nrow(model_extern),2) , "%" ,"===="))
  }
  
  # write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
  print(paste("===",round(100*i/nrow(model_extern),2) , "%" ,"===="))
}

