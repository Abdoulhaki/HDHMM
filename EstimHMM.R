#-------------------------------------------------------------------------------
#Data input
#-------------------------------------------------------------------------------
if(!require(matrixcalc)){install.packages("matrixcalc")}; library(matrixcalc)
if(!require(Rcpp)){install.packages("Rcpp")}; library(Rcpp)
if(!require(RcppArmadillo)){install.packages("RcppArmadillo")}; library(RcppArmadillo)
if(!require(RcppEigen)){install.packages("RcppEigen")}; library(RcppEigen)
if(!require(parallel)){install.packages("parallel")}; library(parallel)
if(!require(doSNOW)){install.packages("doSNOW")}; library(doSNOW)

data_source <- "OxfordMan" #"OxfordMan" or "realized_library" or "11series"
# ".AEX" ".AORD" ".BFX" ".BSESN" ".BVLG" ".BVSP" ".DJI" ".FCHI" ".FTMIB" ".FTSE" ".GDAXI"
# ".GSPTSE" ".HSI" ".IBEX" ".IXIC" ".KS11" ".KSE" ".MXX" ".N225" ".NSEI" ".OMXC20" ".OMXHPI"
# ".OMXSPI" ".OSEAX" ".RUT" ".SMSI" ".SPX" ".SSEC" ".SSMI" ".STI" ".STOXX50E"

current.index <- c(".SPX",".IXIC",".FTSE",".STOXX50E")

# model<-rbind(expand.grid(model=c("MDSV"),r_dis=c("normal"),
#                          rv_dis=c("lognormal"),
#                          N_D=c(2:10),
#                          N=c(1:2),
#                          LEVIER=c(FALSE,TRUE),index=current.index))

para_names<-function(model,LEVIER){
  if(model=="MSM") vars.names<-c("m0","sigma","b","gamma")
  if(model=="FHMV") vars.names<-c("sigma","c1","theta_c","p","m1","theta_m","q")
  if(model=="FHMV_rv") vars.names<-c("sigma","c1","theta_c","p","m1","theta_m","q","shape")
  if(model=="DSARV") vars.names<-c("omega","phi","delta","gamma")
  if(LEVIER) vars.names<-c(vars.names,"l","theta")
  return(vars.names)
}

para_init<-function(model,LEVIER){
  if(model=="MSM") para<-c(0.8,1.26,5,0.05)
  if(model=="FHMV") para<-c(0.4908,4.01718,0.35581,0.99615,8.42876,0.8,0.83108)
  if(model=="DSARV") para<-c(0.31,  0.96, -2.71,  1.53)
  if(LEVIER) para<-c(para,1.89,0.89)
  return(para)
}

ctrl <- list(TOL=1e-15, trace=0)

use.kernel <- TRUE #realized kernel instead of realized variance

start.date <- as.Date("2000-01-03") 
end.date   <- as.Date("2020-06-06")
#end.date   <- as.Date("2019-08-31")
#weekdays(start.date); weekdays(end.date);

DMS <- FALSE
if(DMS){
  path.R    <- "/home/maoudek/Rsim/Article1"
  path.Data <- "/home/maoudek/Rsim/Article1"
  path.Work <- "/home/maoudek/Rsim/Article1/MDSV/Essai"
} else {
  path.R    <- "C:/Users/abdou/Dropbox/Abdoul/Realized volatility/Code"
  path.Work <- "C:/Users/abdou/Dropbox/Abdoul/These/Package"
  path.Data <- "C:/Users/abdou/Dropbox/Abdoul/These/Article1/Data"
}

DellPC <- TRUE
if(DellPC){
  path.R    <- "C:/Users/DellPC/Dropbox/Abdoul/Realized volatility/Code"
  path.Work <- "C:/Users/DellPC/Dropbox/Abdoul/These/Package"
  path.Data <- "C:/Users/DellPC/Dropbox/Abdoul/These/Article1/Data"
}

# rel_tol  <- 1e-6 #optimization parameter

setwd(path.R)
source("realized_functions.r") #load functions for realized volatility

index_set<-c(".AEX", ".AORD", ".BFX", ".BSESN", ".BVLG", ".BVSP", ".DJI", ".FCHI", ".FTMIB", ".FTSE", ".GDAXI",
             ".GSPTSE", ".HSI", ".IBEX", ".IXIC", ".KS11", ".KSE", ".MXX", ".N225", ".NSEI", ".OMXC20", ".OMXHPI",
             ".OMXSPI", ".OSEAX", ".RUT", ".SMSI", ".SPX", ".SSEC", ".SSMI", ".STI", ".STOXX50E")

# index_set<-current.index
current.index<-index_set

serie<-"rv"


if(serie=="r"){
  donnees<-list() 
  setwd(path.Data)
  for(i in 1:length(index_set)){
    index<-index_set[i]
    source("realized_data_traitement.r") #load and organize realized volatility data
    
    RV<-cbind(r=100*Data$RV$r)
    names(RV)<-as.Date(rownames(Data$RV))
    
    assign(index_set[i],RV)
    donnees[[i]]<-get(index_set[i]) 
  }
  
  #=======================================================
  
  setwd(paste0(path.Work,"/HMM"))
  
  for(model in c("DSARV","FHMV")){
    sourceCpp(paste0(model,".cpp"))
    
    filename <- paste(model,"_LEVIER","_FALSE_",start.date,"_",end.date, sep="")
    
    model_file<-expand.grid(model=model,N=c(10),LEVIER=c(FALSE,TRUE),index=current.index)
    vars <- c('N_T','Np',para_names(model,TRUE),'loglik', 'AIC', 'BIC', 'times')
    model_add <- matrix(0, nrow=nrow(model_file), ncol=length(vars))
    colnames(model_add) <- vars
    model_file <- cbind(model_file, model_add)
    
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    
    n_model <-nrow(model_file)
    registerDoSNOW(cl)
    pb <- txtProgressBar(max = n_model, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n_model)
    opts <- list(progress = progress)
    
    Y<-foreach(i=1:n_model, .combine=rbind, .export=c("solnp"), .packages = c("Rcpp","RcppArmadillo","RcppEigen"),
               .noexport = c("natWork","workNat","logLik"), .options.snow = opts) %dopar% {
                 
                 sourceCpp(paste0(model,".cpp")) 
                 T1<-Sys.time()
                 index<-as.character(model_file[i,"index"])
                 donne<-donnees[[which(index_set==index)]]
                 donne<-donne-mean(donne)
                 
                 N<-as.numeric(model_file[i,"N"])
                 LEVIER<-as.logical(model_file[i,"LEVIER"])
                 
                 para<-para_init(model,LEVIER)
                 
                 para_tilde<-natWork(para,LEVIER)
                 opt<-try(solnp(pars=para_tilde,fun=logLik,ech=donne,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
                 
                 vars<-para_names(model,LEVIER)
                 
                 params<-workNat(opt$pars,LEVIER)
                 names(params)<-para_names(model,LEVIER)
                 model_file[i,colnames(model_file) %in% names(params)] <- round(params[vars],5)
                 model_file[i,"loglik"]<--as.numeric(opt$values[length(opt$values)])
                 model_file[i,'N_T']<-nrow(donne)
                 model_file[i,'Np']<-length(params)
                 model_file[i,'AIC'] <- model_file[i,"loglik"]-model_file[i,'Np']
                 model_file[i,'BIC'] <- model_file[i,"loglik"]-(model_file[i,'Np'])*log(model_file[i,'N_T'])/2
                 model_file[i,'times'] <-difftime(Sys.time(),T1,units = "min")
                 
                 model_file[i,]
               }
    
    write.csv(Y, paste(filename,"csv",sep="."), row.names=FALSE)
    
    close(pb)
    # stopCluster(cl)
    on.exit(stopCluster(cl))
  }
} else if(serie=="rv"){
  donnees<-list() 
  setwd(path.Data)
  for(i in 1:length(index_set)){
    index<-index_set[i]
    source("realized_data_traitement.r") #load and organize realized volatility data
    
    RV<-cbind(RV=10000*Data$RV$rv,r=100*Data$RV$r)
    names(RV)<-as.Date(rownames(Data$RV))
    
    assign(index_set[i],RV)
    donnees[[i]]<-get(index_set[i]) 
  }
  
  #=======================================================
  
  setwd(paste0(path.Work,"/HMM"))
  
  model<-"FHMV_rv"
  sourceCpp("FHMV_rv.cpp")
  
  filename <- paste("FHMV_rv",start.date,"_",end.date, sep="")
  
  model_file<-expand.grid(model=model,N=c(10),LEVIER=c(FALSE,TRUE),index=current.index)
  vars <- c('N_T','Np',para_names(model,TRUE),'loglik', 'AIC', 'BIC', 'times')
  model_add <- matrix(0, nrow=nrow(model_file), ncol=length(vars))
  colnames(model_add) <- vars
  model_file <- cbind(model_file, model_add)
  
  #setup parallel backend to use many processors
  cores=detectCores()
  cl <- makeCluster(cores[1]-1) #not to overload your computer
  
  n_model <-nrow(model_file)
  registerDoSNOW(cl)
  pb <- txtProgressBar(max = n_model, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n_model)
  opts <- list(progress = progress)

  Y<-foreach(i=1:n_model, .combine=rbind, .export=c("solnp"), .packages = c("Rcpp","RcppArmadillo","RcppEigen"),
             .noexport = c("natWork","workNat","logLik"), .options.snow = opts) %dopar% {

               sourceCpp("FHMV_rv.cpp")
  # for(i in 1:n_model){
               T1<-Sys.time()
               index<-as.character(model_file[i,"index"])
               donne<-donnees[[which(index_set==index)]]
               donne[,"r"]<-donne[,"r"]-mean(donne[,"r"])
               
               N<-as.numeric(model_file[i,"N"])
               LEVIER<-as.logical(model_file[i,"LEVIER"])
               
               para<-c(0.4908,4.01718,0.35581,0.99615,8.42876,0.81,0.83108,4.0)
               if(LEVIER) para<-c(para,1.89,0.89)
               
               para_tilde<-natWork(para,LEVIER)
               opt<-try(solnp(pars=para_tilde,fun=logLik,ech=donne,LEVIER=LEVIER,N=N,Nl=70,control=ctrl),silent=T)
               
               vars<-para_names(model,LEVIER)
               
               params<-workNat(opt$pars,LEVIER)
               names(params)<-para_names(model,LEVIER)
               model_file[i,colnames(model_file) %in% names(params)] <- round(params[vars],5)
               model_file[i,"loglik"]<--as.numeric(opt$values[length(opt$values)])
               model_file[i,'N_T']<-nrow(donne)
               model_file[i,'Np']<-length(params)
               model_file[i,'AIC'] <- model_file[i,"loglik"]-model_file[i,'Np']
               model_file[i,'BIC'] <- model_file[i,"loglik"]-(model_file[i,'Np'])*log(model_file[i,'N_T'])/2
               model_file[i,'times'] <-difftime(Sys.time(),T1,units = "min")
               # print(paste0(i/n_model," %"))

               model_file[i,]
  }
  
  write.csv(Y, paste(filename,"csv",sep="."), row.names=FALSE)
  
  close(pb)
  # stopCluster(cl)
  on.exit(stopCluster(cl))
}
