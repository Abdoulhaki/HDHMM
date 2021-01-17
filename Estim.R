#-------------------------------------------------------------------------------
#Data input
#-------------------------------------------------------------------------------
if(!require(matrixcalc)){install.packages("matrixcalc")}; library(matrixcalc)
if(!require(Rcpp)){install.packages("Rcpp")}; library(Rcpp)
if(!require(RcppArmadillo)){install.packages("RcppArmadillo")}; library(RcppArmadillo)
if(!require(RcppEigen)){install.packages("RcppEigen")}; library(RcppEigen)
if(!require(DEoptim)){install.packages("DEoptim")}; library(DEoptim)

data_source <- "OxfordMan" #"OxfordMan" or "realized_library" or "11series"
# ".AEX" ".AORD" ".BFX" ".BSESN" ".BVLG" ".BVSP" ".DJI" ".FCHI" ".FTMIB" ".FTSE" ".GDAXI"
# ".GSPTSE" ".HSI" ".IBEX" ".IXIC" ".KS11" ".KSE" ".MXX" ".N225" ".NSEI" ".OMXC20" ".OMXHPI"
# ".OMXSPI" ".OSEAX" ".RUT" ".SMSI" ".SPX" ".SSEC" ".SSMI" ".STI" ".STOXX50E"

current.index <- c(".SPX",".IXIC",".FTSE")

# model<-rbind(expand.grid(model=c("MDSV"),r_dis=c("normal"),
#                          rv_dis=c("lognormal"),
#                          N_D=c(10),
#                          N=c(3),
#                          LEVIER=c(FALSE,TRUE),index=current.index),
#              expand.grid(model=c("MDSV"),r_dis=c("normal"),
#                          rv_dis=c("lognormal"),
#                          N_D=c(5),
#                          N=c(4),
#                          LEVIER=c(FALSE,TRUE),index=current.index),
#              expand.grid(model=c("MDSV"),r_dis=c("normal"),
#                          rv_dis=c("lognormal"),
#                          N_D=c(2),
#                          N=c(7:10),
#                          LEVIER=c(FALSE,TRUE),index=current.index))
model<-rbind(expand.grid(model=c("MDSV"),r_dis=c("normal"),
                         rv_dis=c("lognormal"),
                         N_D=c(1024),
                         N=c(1),
                         LEVIER=c(TRUE),index=current.index))

vars <- c('N_T','Np',"omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape","l","theta",'Marg_loglik','loglik', 'AIC', 'BIC', 'times')
model_add <- matrix(0, nrow=nrow(model), ncol=length(vars))
colnames(model_add) <- vars
model <- cbind(model, model_add)

ctrl <- list(TOL=1e-15, trace=0)

para_names<-function(LEVIER){
  vars.names<-c("omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape")
  if(LEVIER) vars.names<-c(vars.names,"l","theta")
  return(vars.names)
}

use.kernel <- TRUE #realized kernel instead of realized variance

start.date <- as.Date("2000-01-03") 
end.date   <- as.Date("2019-08-31")
#end.date   <- as.Date("2015-09-23")
#weekdays(start.date); weekdays(end.date);

DMS <- TRUE
if(DMS){
  path.R    <- "/home/maoudek/Rsim/Article1"
  path.Data <- "/home/maoudek/Rsim/Article1"
  path.Work <- "/home/maoudek/Rsim/Article1/MDSV/Essai"
} else {
  path.R    <- "C:/Users/abdou/Dropbox/Abdoul/Realized volatility/Code"
  path.Work <- "C:/Users/abdou/Dropbox/Abdoul/These/Package"
  path.Data <- "C:/Users/abdou/Dropbox/Abdoul/These/Article1/Data"
}

# rel_tol  <- 1e-6 #optimization parameter

setwd(path.R)
source("realized_functions.r") #load functions for realized volatility

index_set<-c(".AEX", ".AORD", ".BFX", ".BSESN", ".BVLG", ".BVSP", ".DJI", ".FCHI", ".FTMIB", ".FTSE", ".GDAXI",
             ".GSPTSE", ".HSI", ".IBEX", ".IXIC", ".KS11", ".KSE", ".MXX", ".N225", ".NSEI", ".OMXC20", ".OMXHPI",
             ".OMXSPI", ".OSEAX", ".RUT", ".SMSI", ".SPX", ".SSEC", ".SSMI", ".STI", ".STOXX50E")

index_set<-current.index

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

filename <- paste("Estim_MDSVBinomial_N_1_ND_1024_",paste(current.index,collapse = "_"),"_",start.date,"_",end.date, sep="")

for(i in 1:nrow(model)){
  T1<-Sys.time()
  index<-as.character(model[i,"index"])
  donne<-donnees[[which(index_set==index)]]
  donne[,"r"]<-donne[,"r"]-mean(donne[,"r"])
  
  N_D<-as.numeric(model[i,"N_D"])
  N<-as.numeric(model[i,"N"])
  r_dis<-as.character(model[i,"r_dis"])
  rv_dis<-as.character(model[i,"rv_dis"])
  LEVIER<-as.numeric(model[i,"LEVIER"])
  para<-c(0.75,0.99, 5.0,	1.1,	0.72,	0.5,	0.8,	-0.09,	0.09,	0.2)
  
  if(LEVIER) {
    para<-c(para,0.2 ,0.99)
  }
  para_tilde<-natWork(para,LEVIER)
  opt<-try(solnp(pars=para_tilde,fun=logLik,ech=donne,N_D=N_D,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
  
  vars<-c("omega","a","b","sigma","v0","xi","varphi","delta1","delta2","shape")
  if(LEVIER) {
    vars<-c(vars,"l","theta")
  }
  
  params<-workNat(opt$pars,LEVIER)
  names(params)<-para_names(LEVIER)
  model[i,colnames(model) %in% names(params)] <- round(params[vars],5)
  
  l<-logLik2(donne,opt$pars,LEVIER,N_D,N,t=nrow(donne))
  model[i,"loglik"]<--as.numeric(opt$values[length(opt$values)])
  model[i,"Marg_loglik"]<-l$Marg_loglik
  model[i,'N_T']<-nrow(donne)
  model[i,'Np']<-length(params)
  model[i,'AIC'] <- model[i,"loglik"]-model[i,'Np']
  model[i,'BIC'] <- model[i,"loglik"]-(model[i,'Np'])*log(model[i,'N_T'])/2
  model[i,'times'] <-difftime(Sys.time(),T1,units = "min")
  write.csv(model, paste(filename,"csv",sep="."), row.names=FALSE)
  print(paste("===",round(100*i/nrow(model),2) , "%" ,"==== LOGLIK = ", model[i,"loglik"], "==="))
}


