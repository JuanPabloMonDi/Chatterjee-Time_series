### Load required packages ------------------------

# List of packages to check
required_packages <- c(
  "ppcor", "mgcv", "gplm", "FOCI", "dplyr", "tidyr", 
  "ggplot2", "XICOR", "forecast", "minerva", "tsDyn", "glmnet","rugarch","astsa","Rcpp"
)

# Identify which packages are missing
missing_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]

# Install missing packages if there are any
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
  message("Installed packages: ", paste(missing_packages, collapse = ", "))
} else {
  message("All packages are already installed!")
}


library(ppcor) #Partial correlation function
library(mgcv)
library(gplm)
library(FOCI) #Codec
library(rugarch) #To simulate garch
library(astsa) #To simulate sarima
library(dplyr)
library(tidyr)
library(ggplot2)
library(XICOR) #Xi correlation coefficient
library(forecast)
library(minerva)
library(tsDyn)
library(glmnet)

source("Functions.R") #Import some functions previously defined in the Functions.R file

# Define the folder path
folder_path <- "Results_sim/Results_simB" #For first 100 values 
#folder_path <- "Results_sim/Results_simC" #For last 100 values

# Check if the folder does NOT exist, then create it
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
}


### set seed for replicability -------------------

#For first 100 values (ResultsB)
set.seed(13072026) #Random seed, it was established due the first extension deadline
#For last 100 values (ResultsC)
#set.seed(02082026) #Random seed, it was established due the second extension deadline  


### Simulate time series -----------------------



### ------------------------- Simulations ----------------------------


#### Initial parameters -----------------------------------

#Fixed terms
n_sim=100 #Number os simulations
lag.max1=30 #Minimum of accepted las

# Association measures
pp=NULL
partp=NULL
psp=NULL
pk=NULL
pxi=NULL
pcodec<-NULL

# AIC/BIC Measures
paic  <- NULL   
pbic  <- NULL   
#LASSO regression measures
#plasso <- NULL

sizes<-(c(100,500,1000,2000,5000)) #Sample sizes to include in the simulations
n_size=0
for (size in sizes){
  n_size=n_size+1
  lag.max=max(lag.max1,round(12*(size/100)^{1/4})) #Schwert's Rule, to avoid a small sample of lags, we set a minimum of 30
  
  #### loop of Simulations---------------------------------
  for (sim in 1:n_sim){
    
    
    ##### Time series ----------------------------------------------------
    #Now, we are going to simulate the time series 
    
    #1. AR(5) 
    
    x1<-arima.sim(n=size,list(ar=c(0.5,-0.2,0.1,0.2,0.1)))
    
    #2. ARMA(1,3)
    x2<-arima.sim(n=size,list(ar=c(-0.75),ma=c(-0.5,0.5,.15)))
    
    #3.ARIMA(2,1,1)
    x3<-arima.sim(list(order=c(2,1,1),ar=c(0.7,-0.5),ma=0.4),n=size-1)
    
    #4. Non Linear  AR(3) with sine functions
    x4<-rnorm(3,0,20)
    for (i in 4:size){
      x4<-append(x4,3*sin(x4[i-1])+2*sin(x4[i-2])+sin(x4[i-3])+rnorm(1))
    }
    
    #5. Non linear AR(4) with sine and logistic function
    x5<-rnorm(4,0,1)
    for (i in 5:size){
      x5<-append(x5,3*sin(x5[i-1])+2*sin(x5[i-2]/3)+0.5*sin(x5[i-3]/2)-3*1/(1+exp(x5[i-4]))+rnorm(1))
    }
    
    
    #6. Non linear ARMA (2,2)
    x6<-rnorm(2,2,10)
    error<-rnorm(2)
    for (i in 3:size){
      error2<-rnorm(1,0,0.2)
      nuevox<-2*cos(x6[i-1])+0.5*sin(x6[i-2])+0.4*error[i-1]+0.8*1/(1+exp(error[i-2]))+error2
      error<-append(error,error2)
      x6<-append(x6,nuevox)
    }
    
    #7. SETAR(2;2,2)
    TvarMat <- c(2.9,-0.4,-0.1,-1.5, 0.2,0.3)
    x7<-setar.sim(n=size,B=TvarMat,lag=2, type="simul", nthresh=1, Thresh=2, starting=c(2.8,2.2))
    
    
    #8. ARIMA (2,1,1) with first order difference
    x8<-arima.sim(list(order=c(2,1,1),ar=c(0.7,-0.5),ma=0.4),n=size)
    x8<-diff(x8,lag=1)
    
    #9. SAR(4)_(2,0,0)_{12}
    x9<-sarima.sim(ar=c(0.6, -0.25, 0.1, -0.05),
                   S=12,
                   sar = c(.4,-.15),
                   burnin = 100,
                   n=6*size/5)
    
    #10. SAR additive decomposed
    x10<-decompose(x9,type="additive")
    x10<-x10$random[!is.na(x10$random)]
    x10<-tail(x10,size)
    
    #11. SAR multiplicative decomposed
    x11<- decompose(x9,type="multiplicative")
    x11<-x11$random[!is.na(x11$random)]
    
    x9<-tail(x9,size)
    x10<-tail(x10,size)
    x11<-tail(x11,size)
    
    
    #12. ARMA(1,0,1)-GARCH(1,1)
    spec <- ugarchspec(
      variance.model = list( model = "sGARCH", garchOrder = c(1, 1)),
      mean.model = list(armaOrder = c(1, 1),include.mean = FALSE),
      distribution.model = "norm",
      fixed.pars = list(
        ar1    = 0.75,
        ma1    = 0.5,
        omega  = 0.00001,
        alpha1 = 0.05,
        beta1  = 0.9)
    )
    
    sim_path <- ugarchpath(
      spec,
      n.sim   = size,
      n.start = 100,
      m.sim   = 1
    )
    
    x12<- fitted(sim_path)[, 1]             
    
    
    
    #13. ARIMA(3,1,0)
    x13<-arima.sim(list(order=c(3,1,0),ar=c(0.6,-0.25,.3)),n=size)
    
    #14. ARIMA(3,1,0) differentiated
    x14<-diff(x13,lag=1)
    x13<-tail(x13,size)
    
    
    series<-data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,
                       x8=x8,x9=x9,x10=x10,x11=x11,x12=x12,x13=x13,x14=x14)
    n_serie=0
    
    ###### Loop of series ------------------------
    for (serie in colnames(series)){
      x<-series[[serie]]
      n_serie<-n_serie+1
      
      corPearson<-c()
      corPartPearson<-c()
      corSpearman<-c()
      corKendall<-c()
      corXi<-c()
      
      for (lag in 1:lag.max){
        #This is just a function to create a progress bar
        update_progress((n_size-1)*(n_sim)*ncol(series)*lag.max+(sim-1)*lag.max*ncol(series) + (n_serie - 1)*lag.max + lag, length(sizes)*n_sim*ncol(series)*lag.max)
        data<-data.frame()
        
        #Construction of the matrix with lag columns
        for (i in 0:(length(x)-lag-1)){
          xt<-x[(length(x)-i)]
          xt1<-t(data.frame(rev(x[(length(x)-(lag)-i):(length(x)-(i+1))])))
          data<-rbind(data,cbind(xt,xt1))
        }
        #Modify the data frame created
        row.names(data)<-0:(nrow(data)-1)
        VarReg<-paste0("x",1:(lag))
        colnames(data)<-c("xt",VarReg)
        
        
        #calculate the partial correlation matrix for each correlation
        
        corP<-pcor(data,method="pearson")
        corSp<-pcor(data,method="spearman")
        if (size<=2000){
          corK<-pcor(data,method="kendall")}
        corX<-partcor(data,method = "xi")
        
        
        #Since we are only interested in the correlation with xt (the last time observed), we only consider the first row
        
        #The pearson and Spearman will be made by the pacf function because it's faster  
        corP<-corP$estimate[1,]
        corSp<-corSp$estimate[1,]
        if (size<=2000){
          corK<-corK$estimate[1,]}
        corX<-corX[1,]
        
        corPearson<-append(corPearson, corP[lag+1])
        corSpearman<-append(corSpearman,corSp[lag+1])
        if (size<=2000){corKendall<-append(corKendall,corK[lag+1])}
        corXi<-append(corXi,corX[lag+1])
      } #End of lag.max
      
      ##### AIC/BIC -------------------------------------------------------
      # AIC/BIC-based AR(p) order selection 
      
      #Get the p lag that minimizes the AIC and BIC
      aic_vals <- rep(NA_real_, lag.max)
      bic_vals <- rep(NA_real_, lag.max)
      for (p in 1:lag.max){
        #Adjust the serie to an AR(p) process
        fit <- tryCatch(
          arima(x, order=c(p,0,0), method="ML"),
          error = function(e) NULL
        )
        if (!is.null(fit)){
          aic_vals[p] <- AIC(fit)
          bic_vals[p] <- BIC(fit)
        } else {
          aic_vals[p] <- NA
          bic_vals[p] <- NA
        }
      }
      p_aic <- which.min(aic_vals)
      p_bic <- which.min(bic_vals)
      paic <- rbind(paic, data.frame(size=size, sim=sim, serie=serie, p1=p_aic))
      pbic <- rbind(pbic, data.frame(size=size, sim=sim, serie=serie, p1=p_bic))
      
      
      ##### Penalized Regression - LASSO ------------------------------------------------
      
      # Penalized regression (LASSO) 
      #fit_lasso <- cv.glmnet(
      #  x = as.matrix(data[, -1]),
      #  y = data$xt,
      #  alpha = 1
      #)
      
      #Identify significant coefficients
      #coefs_lasso <- coef(fit_lasso, s = "lambda.min")[-1]
      #lasso_lags <- which(coefs_lasso != 0)
      #nlasso <- length(lasso_lags)
      
      #Get three hights orders
      #plasso <- rbind(plasso, data.frame(
      #  size  = size,
      #  sim   = sim,
      #  serie = serie,
      #  p1 = ifelse(nlasso >= 1, max(lasso_lags), 0),
      #  p2 = ifelse(nlasso >= 2, sort(lasso_lags, TRUE)[2], 0),
      #  p3 = ifelse(nlasso >= 3, sort(lasso_lags, TRUE)[3], 0)
      #))
      
      
      
      
      
      
      ##### Confidence intervals ---------------------------------------------------
      #corSpearman<-pacf(rank(series[serie]),plot=F,lag.max = lag.max)
      #corSpearman<-corSpearman$acf
      #For Pearson the variance is 1/N
      Plags<-data.frame(lag=(index(corPearson)-1),value=corPearson)
      #plot(corPearson[1:lag],type="h")
      #abline(h=qnorm(0.975)*sqrt(1/size))
      #abline(h=-qnorm(0.975)*sqrt(1/size))
      Plags<-Plags[abs(Plags$value)>qnorm(0.975)*sqrt(1/size),]
      Plags<-Plags$lag
      npp<-length(Plags)
      pp<-rbind(pp,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(npp>=1, Plags[[length(Plags)]], 0),p2=ifelse(npp>=2, Plags[(length(Plags)-1)], 0),p3=ifelse(npp>=3, Plags[(length(Plags)-2)], 0)))
      
      
      # For PACF with pacf() function (its should be approx to Plags, but diference can exsts due the method used (pacf uses Durbin Levinson))
      corPartPearson<-pacf(x,plot=F,lag.max = lag.max)
      corPartPearson<-corPartPearson$acf
      corPartPearson_lags <- which(abs(corPartPearson) > qnorm(0.975)*sqrt(1/size))
      npartp<-length(corPartPearson_lags)
      partp<-rbind(partp,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(npartp>=1, corPartPearson_lags[[length(corPartPearson_lags)]], 0),p2=ifelse(npartp>=2, corPartPearson_lags[(length(corPartPearson_lags)-1)], 0),p3=ifelse(npartp>=3, corPartPearson_lags[(length(corPartPearson_lags)-2)], 0)))
      
      
      #For Spearman the variance is 1/(N-1)
      #plot(corSpearman,type="h")
      #abline(h=-qnorm(0.975)*sqrt(1/(size-1)))
      #abline(h=qnorm(0.975)*sqrt(1/(size-1)))
      
      SPlags<-data.frame(lag=(index(corSpearman)-1),value=corSpearman)
      SPlags<-SPlags[abs(SPlags$value)>qnorm(0.975)*sqrt(1/(size-1)),]
      SPlags<-SPlags$lag
      nsp<-length(SPlags)
      psp<-rbind(psp,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(nsp>=1, SPlags[[length(SPlags)]], 0),p2=ifelse(nsp>=2, SPlags[(length(SPlags)-1)], 0),p3=ifelse(nsp>=3, SPlags[(length(SPlags)-2)], 0)))
      
      
      #For Kendall the variance is 2(2n+5)/9n(n-1)
      #plot(corKendall[1:lag],type="h")
      #abline(h=qnorm(0.975)*sqrt(2*(2*size+5)/(9*size*(size-1))))
      #abline(h=-qnorm(0.975)*sqrt(2*(2*size+5)/(9*size*(size-1))))
      if (size<=2000){ #To avoid extensive computing
        Klags<-data.frame(lag=(index(corKendall)-1),value=corKendall)
        Klags<-Klags[abs(Klags$value)>(qnorm(0.975)*sqrt(2*(2*size+5)/(9*size*(size-1)))),]
        Klags<-Klags$lag
        nk<-length(Klags)
        pk<-rbind(pk,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(nk>=1, Klags[[length(Klags)]], 0),p2=ifelse(nk>=2, Klags[[(length(Klags)-1)]], 0),p3=ifelse(nk>=3, Klags[[(length(Klags)-2)]], 0)))
      }
      #For Xi the variance is 2/5N
      #plot(corXi[1:lag],type="h")
      #abline(h=qnorm(0.975)*sqrt(2/(5*size)),col="red",lw=1)
      #abline(h=-qnorm(0.975)*sqrt(2/(5*size)),col="red",lw=1)
      
      #Get the 3 maximum autocorrelations
      Xilags<-data.frame(lag=(index(corXi)-1),value=corXi)
      Xilags<-Xilags[abs(Xilags$value)>(qnorm(0.975)*sqrt(2/(5*size))),]
      Xilags<-Xilags$lag#[(nrow(Xilags)-2):nrow(Xilags)]
      nxi=length(Xilags)
      pxi<-rbind(pxi,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(nxi>=1, Xilags[[length(Xilags)]], 0),p2=ifelse(nxi>=2, Xilags[[(length(Xilags)-1)]], 0),p3=ifelse(nxi>=3, Xilags[[(length(Xilags)-2)]], 0)))
      
      #Now, with the coefficient of Azadka-Chatterjee
      codecM<-foci(data$xt,data[2:(lag+1)],numCores = 1)
      codecM<-codecM$selectedVar$index
      ncodec=length(codecM)
      pcodec<-rbind(pcodec,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(ncodec>=1, max(codecM), 0),p2=ifelse(ncodec>=2, sort(codecM, TRUE)[2], 0),p3=ifelse(ncodec>=3,sort(codecM, TRUE)[3], 0)))
    } #End series iterations
    if(sim%%15==0 | sim==n_sim){
      ## save results in different files
      write.csv(pp, file=paste0(folder_path,"/Results_Pearson.csv"),row.names = FALSE)
      write.csv(partp, file=paste0(folder_path,"/Results_Pearson_pacf.csv"),row.names = FALSE)
      write.csv(psp, file=paste0(folder_path,"/Results_Spearman.csv"),row.names = FALSE)
      write.csv(pk, file=paste0(folder_path,"/Results_Kendall.csv"),row.names = FALSE)
      write.csv(pxi, file=paste0(folder_path,"/Results_Chatterjee.csv"),row.names = FALSE)
      write.csv(paic, file=paste0(folder_path,"/Results_AIC.csv"),row.names = FALSE)
      write.csv(pbic, file=paste0(folder_path,"/Results_BIC.csv"),row.names = FALSE)
      write.csv(pcodec, file=paste0(folder_path,"/Results_FOCI.csv"),row.names = FALSE)
    }
  } #End simulations
} #End sizes 

