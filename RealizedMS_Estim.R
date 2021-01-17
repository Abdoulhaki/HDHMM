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

current.index <- c(".SPX", ".IXIC",".FTSE", ".STOXX50E")
# current.index <- c(".STOXX50E")

# model<-rbind(expand.grid(model=c("MDSV"),r_dis=c("normal"),
#                          rv_dis=c("lognormal"),
#                          N_D=c(2:10),
#                          N=c(1:2),
#                          LEVIER=c(FALSE,TRUE),index=current.index))

para_init<-function(LEVIER,K,Model_type){
  para<-c(5.01)
  if(Model_type=="logN") para<-c(c(0.12,1.12,-0.005,0.043),para)
  if(K==2){
    Vol<-c(0.02,3.5)
    P<-matrix(c(0.99,0.01,0.01,0.99),K,K)
  }
  if(K==3){
    Vol<-c(1,2,3)
    P<-matrix(c(0.9,0.05,0.05,0.05,0.9,0.05,0.05,0.05,0.9),K,K)
  }
  if(K==4){
    Vol<-c(0.02,0.05,0.5,1.5)
    P<-matrix(c(0.7,0.1,0.1,0.1,0.1,0.7,0.1,0.1,0.1,0.1,0.7,0.1,0.1,0.1,0.1,0.7),K,K)
  }
  if(LEVIER) para<-c(para,2.75,0.52)
  
  return(list(para=para,P=P,Vol=Vol))
}

f.opt<-function(pas,K,LEVIER,ech,Model_type){
  para1<-para2<-para3<-para_init(LEVIER,K,Model_type)$para
  P1<-P2<-P3<-para_init(LEVIER,K,Model_type)$P
  Vol1<-Vol2<-Vol3<-para_init(LEVIER,K,Model_type)$Vol
  logN<-1
  if(Model_type=="logN") logN<-5
  
  para_tilde<-natWork(para1,Vol1,P1,LEVIER,K,Model_type)
  opt<-try(solnp(pars=para_tilde,fun=logLik,ech=ech,LEVIER=LEVIER,K=K,Model_type=Model_type,Nl=70,control=ctrl),silent=T)
 
  for(i in 1:10){
    if(Model_type=="logN"){
      para1[1]<-para1[1]+pas
      para1[2]<-para1[2]+pas
      para1[3]<-para1[3]+pas
      para1[4]<-para1[4]+pas
    }
    para1[logN]<-para1[logN]+pas
    if(LEVIER){
      para1[logN+1]<-para1[logN+1]+10*pas
      para1[logN+2]<-min(para1[logN+2]+pas/10,0.999)
    }
    P1[diag(K)]<-min(P1[1,1]+pas/10,1)
    P1[!diag(K)]<-(1-P1[1,1])/(K-1)
    Vol1<-Vol1+pas
    
    para_tilde<-natWork(para1,Vol1,P1,LEVIER,K,Model_type)
    opt1<-try(solnp(pars=para_tilde,fun=logLik,ech=ech,LEVIER=LEVIER,K=K,Model_type=Model_type,Nl=70,control=ctrl),silent=T)
  
    if(class(opt) == "try-error"){
      if(!(class(opt1) == "try-error")){
        opt<-opt1
      }
    }else{
      if(!(class(opt1) == "try-error")){
        if(-as.numeric(opt1$values[length(opt1$values)])>-as.numeric(opt$values[length(opt$values)]))
          opt<-opt1
      }
    }
    
    if(Model_type=="logN"){
      para2[1]<-para2[1]-pas
      para2[2]<-para2[2]-pas
      para2[3]<-para2[3]-pas
      para2[4]<-para2[4]-pas
    }
    para2[logN]<-pmax(0.001,para2[logN]-pas)
    if(LEVIER){
      para2[logN+1]<-max(0.001,para2[logN+1]-10*pas)
      para2[logN+2]<-max(para2[logN+2]-pas/10,0.001)
    }
    P2[diag(K)]<-max(P2[1,1]-pas/10,0)
    P2[!diag(K)]<-(1-P2[1,1])/(K-1)
    Vol2<-pmax(0.001,Vol2-pas)
    
    para_tilde<-natWork(para2,Vol2,P2,LEVIER,K,Model_type)
    opt1<-try(solnp(pars=para_tilde,fun=logLik,ech=ech,LEVIER=LEVIER,K=K,Model_type=Model_type,Nl=70,control=ctrl),silent=T)
    
    if(class(opt) == "try-error"){
      if(!(class(opt1) == "try-error")){
        opt<-opt1
      }
    }else{
      if(!(class(opt1) == "try-error")){
        if(-as.numeric(opt1$values[length(opt1$values)])>-as.numeric(opt$values[length(opt$values)]))
          opt<-opt1
      }
    }
    
    if(Model_type=="logN"){
      para3[1]<-para3[1]-pas
      para3[2]<-para3[2]-pas
      para3[3]<-para3[3]-pas
      para3[4]<-para3[4]-pas
    }
    para3[logN]<-max(0.001,para3[logN]+pas)
    if(LEVIER){
      para3[logN+1]<-para3[logN+1]+10*pas
      para3[logN+2]<-max(para3[logN+2]-pas/20,0.001)
    }
    P3[diag(K)]<-max(P3[1,1]-pas/10,0)
    P3[!diag(K)]<-(1-P3[1,1])/(K-1)
    Vol3<-pmax(0.001,Vol3+pas)
    
    para_tilde<-natWork(para3,Vol3,P3,LEVIER,K,Model_type)
    opt1<-try(solnp(pars=para_tilde,fun=logLik,ech=ech,LEVIER=LEVIER,K=K,Model_type=Model_type,Nl=70,control=ctrl),silent=T)
    
    if(class(opt) == "try-error"){
      if(!(class(opt1) == "try-error")){
        opt<-opt1
      }
    }else{
      if(!(class(opt1) == "try-error")){
        if(-as.numeric(opt1$values[length(opt1$values)])>-as.numeric(opt$values[length(opt$values)]))
          opt<-opt1
      }
    }
  }
  opt
}

para_name<-function(params,K,LEVIER,Model_type){
  if(Model_type=="logN") {
    vars.names<-c("xi","varphi","delta1","delta2","shape")
    }else {vars.names<-c("shape")}
  if(LEVIER) vars.names<-c(vars.names,'l','theta')
  
  names(params$para)<-vars.names
  names(params$Vol)<-paste('Vol',1:K)
  return(params)
}

ctrl <- list(TOL=1e-15, trace=0)

use.kernel <- TRUE #realized kernel instead of realized variance

start.date <- as.Date("2000-01-03") 
end.date   <- as.Date("2020-06-06")
#end.date   <- as.Date("2015-09-23")
#weekdays(start.date); weekdays(end.date);

DMS <- FALSE
if(DMS){
  path.R    <- "/home/maoudek/Rsim/Article1"
  path.Data <- "/home/maoudek/Rsim/Article1"
  path.Work <- "/home/maoudek/Rsim/Article1/MDSV/Essai"
} else {
  path.R    <- "C:/Users/abdou/Dropbox/Abdoul/Realized volatility/Code"
  path.Work <- "C:/Users/abdou/Dropbox/Abdoul/These/Package/RealMS"
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
sourceCpp(("RealizedMS.cpp"))

filename <- paste("RealMS","_",start.date,"_",end.date, sep="")

model_file<-expand.grid(model="RealMS",K=c(2,4),Model_type=c('InvG','logN'),LEVIER=c(FALSE,TRUE),index=current.index)
vars <- c('N_T','Np',"Marg_loglik",'loglik', 'AIC', 'BIC', 'times','xi','varphi','delta1','delta2','shape', 'l', 'theta', paste('Vol',1:4), paste('P',1:4))
model_add <- matrix(0, nrow=nrow(model_file), ncol=length(vars))
colnames(model_add) <- vars
model_file <- cbind(model_file, model_add)

for(i in 1:nrow(model_file)){
  T1<-Sys.time()
  index<-as.character(model_file[i,"index"])
  donne<-donnees[[which(index_set==index)]]
  donne[,"r"]<-donne[,"r"]-mean(donne[,"r"])
  
  K<-as.numeric(model_file[i,"K"])
  LEVIER<-as.logical(model_file[i,"LEVIER"])
  Model_type<-as.character(model_file[i,"Model_type"])
  pas<-0.5
  
  opt<-f.opt(pas,K,LEVIER,donne,Model_type = Model_type)
  # opt<-try(solnp(pars=para_tilde,fun=logLik,ech=donne,LEVIER=LEVIER,K=K,Nl=70,control=ctrl),silent=T)

  params<-workNat(opt$pars,LEVIER,K,Model_type)
  if(Model_type=="logN") {
    nam<-c("xi","varphi","delta1","delta2","shape")
  }else{ nam<-'shape'}
  if(LEVIER) nam<-c(nam, 'l','theta')
  model_file[i,colnames(model_file) %in% nam] <- params$para
  model_file[i,colnames(model_file) %in% paste('Vol',1:K)] <- params$Vol
  model_file[i,colnames(model_file) %in% paste('P',1:K)] <- diag(params$P)
  l<-logLik2(donne,opt$pars,LEVIER=LEVIER,K=K,t=nrow(donne))
  model_file[i,"Marg_loglik"]<-l$Marg_loglik
  model_file[i,"loglik"]<--as.numeric(opt$values[length(opt$values)])
  model_file[i,'N_T']<-nrow(donne)
  model_file[i,'Np']<-length(opt$pars)
  model_file[i,'AIC'] <- model_file[i,"loglik"]-model_file[i,'Np']
  model_file[i,'BIC'] <- model_file[i,"loglik"]-(model_file[i,'Np'])*log(model_file[i,'N_T'])/2
  model_file[i,'times'] <-difftime(Sys.time(),T1,units = "min")
  print(paste("===",round(100*i/nrow(model_file),2) , "%" ,"==== LOGLIK = ", model_file[i,"loglik"], "===")) 
}

write.csv(model_file, paste(filename,"csv",sep="."), row.names=FALSE)


