### Load required packages ------------------------

# List of packages to check
required_packages <- c(
  "ppcor", "mgcv", "gplm", "FOCI", "dplyr", "tidyr", 
  "ggplot2", "XICOR", "forecast", "minerva", "tsDyn", "glmnet"
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


library(ppcor)
library(mgcv)
library(gplm)
library(FOCI) #Codec
library(dplyr)
library(tidyr)
library(ggplot2)
library(XICOR)
library(forecast)
library(minerva)
library(tsDyn)
library(glmnet)

source("Functions.R")
### set seed for replicability -------------------
set.seed(131101) #Random seed, it was established due the hour the code was created 13:11:01 


### Simulate time series -----------------------





matrizsignos<-function(x){
  if (class(x)=="numeric"){
matrizX<-matrix(0,ncol=length(x),nrow = length(x))}else{if(class(x)=="data.frame"){
  matrizX<-matrix(0,ncol=nrow(x),nrow = nrow(x))
}}
  
n<-0
for (j in 1:nrow(matrizX)){
  for (i in 1:ncol(matrizX)){
    if(j>i){
      if (class(x)=="numeric"){matrizX[i,j]<-sign(x[i]-x[j])}else{matrizX[i,j]<-sign(x[i,]-x[j,])}
    
    n=n+1
    }
    }
}
return(list(matrizX,n))
}

kendallPartial<-function(x,y,z){
n<-matrizsignos(x)[[2]]
Y<-matrizsignos(y)[[1]]
X<-matrizsignos(x)[[1]]
Z<-matrizsignos(z)[[1]]

tau12<-sum(diag(t(X)%*%Y))/n
tau13<-sum(diag(t(X)%*%Z))/n
tau23<-sum(diag(t(Y)%*%Z))/n

taupartial<-(tau12-tau13*tau23)/(sqrt(1-tau13^2)*sqrt(1-tau23^2))
return(taupartial)
}











### ------------------------- Simulations ----------------------------


#### Initial parameters -----------------------------------

#Fixed terms
n_sim=200
lag.max=30

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
plasso <- NULL

#autop<-NULL
#autosp<-NULL
#autok<-NULL
#autoxi<-NULL
#autocodec<-NULL
sizes<-rev(c(100,500,1000,2000,5000))
sizes<-(c(100,500,1000,2000,5000))
n_sim<-100

sizes<-(c(100,500,1000))
n_sim<-10
n_size=0
for (size in sizes){
  n_size=n_size+1
#### loop of Simulations---------------------------------
for (sim in 1:n_sim){
  
  
  ##### Time series ----------------------------------------------------
  #Now, we are going to simulate the time series 
  
  #1. AR(5) 
  
  x1<-arima.sim(n=size,list(ar=c(0.5,-0.2,0.1,0.2,0.1)))
  
  #2. ARMA(1,3)
  x2<-arima.sim(n=size,list(ar=c(-0.75),ma=c(-0.5,1.5,1)))
  
  #3.ARIMA(2,1,1)
  x3<-arima.sim(list(order=c(2,1,1),ar=c(0.7,-0.5),ma=0.4),n=size-1)
  #4. Serie AR(3) con función seno
  x4<-rnorm(3,0,20)
  for (i in 4:size){
    x4<-append(x4,3*sin(x4[i-1])+2*sin(x4[i-2])+sin(x4[i-3])+rnorm(1))
  }
  
  #5. AR(4) con funciones sin(x)
  x5<-rnorm(4,0,1)
  for (i in 4:size){
    x5<-append(x5,3*sin(x5[i-1])+2*sin(x5[i-2]/3)+0.5*sin(x5[i-3]/2)-3*1/(1+exp(x5[i-4]))+rnorm(1))
  }

  
  #6. ARMA FUNCIONAL (2,2) 
  x6<-rnorm(2,2,10)
  error<-rnorm(2)
  for (i in 3:size){
    error2<-rnorm(1,0,0.2)
    nuevox<-2*cos(x6[i-1])+0.5*sin(x6[i-2])+0.4*error[i-1]+0.8*1/(1+exp(error[i-2]))+error2
    error<-append(error,error2)
    x6<-append(x6,nuevox)
  }
  
  #7. Simulamos un SETAR
  TvarMat <- c(2.9,-0.4,-0.1,-1.5, 0.2,0.3)
  x7<-setar.sim(n=size,B=TvarMat,lag=2, type="simul", nthresh=1, Thresh=2, starting=c(2.8,2.2))
  
  
  
  
  series<-data.frame(x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7)
  n_serie=0
  for (serie in colnames(series)){
    x<-series[[serie]]
    n_serie<-n_serie+1
    
    corPearson<-c()
    corPartPearson<-c()
    corSpearman<-c()
    corKendall<-c()
    corXi<-c()
    
    for (lag in 1:lag.max){
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
        
      #We found a faster method for Pearson and Spearman
      corP<-pcor(data,method="pearson")
      corSp<-pcor(data,method="spearman")
      if (size<=2000){
      corK<-pcor(data,method="kendall")}
      corX<-partcor(data,method = "xi")
  
      
      #Now, we calculate the autocorrelation  matrix for each correlation
      
      #autocorP<-cor(data,method="pearson")
      #autocorSp<-cor(data,method="spearman")
      #autocorK<-cor(data,method="kendall")
      #autocorX<-xicor(data)
      
      #Since we are only interested in the correlation with xt, we only consider the first row
      
      #The pearson and Spearman will be made by the pacf function because is faster  
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
  # AIC/BIC-based AR(p) order selection (only meaningful for linear/stationary series x1, x2 and x3)
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
    coefs_lasso <- coef(fit_lasso, s = "lambda.min")[-1]
    lasso_lags <- which(coefs_lasso != 0)
    nlasso <- length(lasso_lags)
    
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
      if (size<=2000){
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
} #End simulations
} #End sizes 


# Define the folder path
folder_path <- "Results"

# Check if the folder does NOT exist, then create it
if (!dir.exists(folder_path)) {
  dir.create(folder_path)
}

## save results in different files
write.csv(pp, file="Results/Results_Pearson.csv",col.names = TRUE,row.names = FALSE)
write.csv(ppart, file="Results/Results_Pearson_pacf.csv",col.names = TRUE,row.names = FALSE)
write.csv(psp, file="Results/Results_Spearman.csv",col.names = TRUE,row.names = FALSE)
write.csv(pk, file="Results/Results_Kendall.csv",col.names = TRUE,row.names = FALSE)
write.csv(pxi, file="Results/Results_Chatterjee.csv",col.names = TRUE,row.names = FALSE)
write.csv(paic, file="Results/Results_AIC.csv",col.names = TRUE,row.names = FALSE)
write.csv(pbic, file="Results/Results_BIC.csv",col.names = TRUE,row.names = FALSE)
write.csv(pcodec, file="Results/Results_FOCI.csv",col.names = TRUE,row.names = FALSE)
cat("\n")
#Iteraciones orden size>n_sim>series>lag



pp_long <- tidyr::pivot_longer(pp, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
psp_long <- tidyr::pivot_longer(psp, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
pk_long <- tidyr::pivot_longer(pk, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")

pxi_long <- tidyr::pivot_longer(pxi, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")

pcodec_long <- tidyr::pivot_longer(pcodec, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")


partp_long<- tidyr::pivot_longer(partp, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
paic_long <- tidyr::pivot_longer(paic, cols = c(p1), names_to = "Estimator", values_to = "Value")
pbic_long <- tidyr::pivot_longer(pbic, cols = c(p1), names_to = "Estimator", values_to = "Value")
plasso_long<- tidyr::pivot_longer(plasso, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")

### ------------------- Plots of the series -----------------------



#Titles for the plots
labs_est <- c(
  p1 = "Max~lag~p[1]",
  p2 = "Second~Max~lag~p[2]",
  p3 = "Third~Max~lag~p[3]"
)


comparison_plot<-function(series,estimators,inter, data,factor_plot="size"){
  plot<-ggplot(data%>%filter(serie==series),aes(x = factor(!!sym(factor_plot)), y = Value)) +
    geom_violin(trim = F,   
                fill = "grey60",
                color = "grey40",
                alpha = 0.25,
                lwd=.6) +
    geom_boxplot(width = 0.1, outlier.size = 0.9,
                 fill = "gray20",
                 color = "gray20",
                 lwd=.9,
                 median.colour = "white")+
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 25,     
      size = 2.5,
      fill = "white",  
      color = "black"
    )+
    geom_hline(yintercept = inter, colour = "black", lwd = 1, linetype=2, alpha=0.5) +
    labs(x="Size",y="Lag") +
    facet_wrap(
      vars(!!sym(estimators)),
      scales = "free_y",
      #labeller = labeller(!!estimators := as_labeller(labs_est, label_parsed))
    ) +
    theme_linedraw(paper = "white", ink = "gray20") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size=12)
    )
  return(plot)
}



#### ------------------------------- Serie 1 -------------------------

ggplot(pp_long[pp_long$serie=="x1",], aes(x = factor(size), y = Value)) +
  geom_violin(trim = F,   
              fill = "grey60",
              color = "grey40",
              alpha = 0.25,
              lwd=.6) +
  geom_boxplot(width = 0.1, outlier.size = 0.9,
               fill = "gray20",
               color = "gray20",
               lwd=.9,
               median.colour = "white") +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 25,       # triángulo hacia arriba
    size = 2.5,
    fill = "white",   # relleno del triángulo
    color = "black"
  )+
  geom_hline(yintercept = 5, colour = "black", lwd = 1, linetype=2, alpha=0.5) +
  labs(x="Size",y="Lag") +
  facet_wrap(
    ~Estimator,
    scales = "free_y",
    labeller = labeller(Estimator = as_labeller(labs_est, label_parsed))
  ) +
  theme_linedraw(paper = "white", ink = "gray20") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(size=12)
  )

ggplot(psp_long[psp_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 5,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 5,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ 
  geom_hline(yintercept = 5,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 5,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))


ggplot(partp_long[partp_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 5,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(paic_long[paic_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 5,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(pbic_long[pbic_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 5,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(plasso_long[plasso_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 5,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))



estimators<-"name"
data<-rbind(pp_long%>%mutate(name="Pearson"),pcodec_long%>%mutate(name="CODEC"), pk_long%>%mutate(name="Kendall"))
comparison_plot("x7","name",data=data, inter=2)
#### --------------------- Serie 2 ---------------------------

ggplot(pp_long[pp_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 1,colour="blue",lwd=1.25)+
  xlab("Size")+
  facet_wrap(~Estimator,scales = "free_y")+
  theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 1,colour="blue",lwd=1.25)+
  xlab("Size")+
  facet_wrap(~Estimator,scales = "free_y")+
  theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 1,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 1,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 1,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

#### ---------------------- Serie 3 --------------------------
ggplot(pp_long[pp_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))


#### ----------------------- Serie 4 --------------------

ggplot(pp_long[pp_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))



ggplot(partp_long[partp_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(paic_long[paic_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(pbic_long[pbic_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(plasso_long[plasso_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))


#### ----------------------- Serie 5 -------------------

ggplot(pp_long[pp_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))



ggplot(partp_long[partp_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(paic_long[paic_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(pbic_long[pbic_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(plasso_long[plasso_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

#### --------------- Serie 6 ------------------


ggplot(pp_long[pp_long$serie=="x6",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x6",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x6",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(pxi_long[pxi_long$serie=="x6",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x6",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))



ggplot(partp_long[partp_long$serie=="x6",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(paic_long[paic_long$serie=="x6",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(pbic_long[pbic_long$serie=="x6",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(plasso_long[plasso_long$serie=="x6",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

#### --------------- Serie 7 ------------------

ggplot(pp_long[pp_long$serie=="x7",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x7",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(pk_long[pk_long$serie=="x7",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(pxi_long[pxi_long$serie=="x7",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x7",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))



ggplot(partp_long[partp_long$serie=="x7",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(paic_long[paic_long$serie=="x7",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(pbic_long[pbic_long$serie=="x7",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

ggplot(plasso_long[plasso_long$serie=="x7",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

## -----Meausing the error ----------------
  
# Real value of the parameters
valores_reales <- c(5,1,2,3, 4, 2, 2)
names(valores_reales) <- c("x1","x2","x3","x4","x5","x6","x7")




# Add the real parameter values to the table
pp <- pp %>%
  mutate(valor_real = valores_reales[as.character(serie)])
psp <- psp %>%
  mutate(valor_real = valores_reales[as.character(serie)])
pk <- pk %>%
  mutate(valor_real = valores_reales[as.character(serie)])
pxi <- pxi %>%
  mutate(valor_real = valores_reales[as.character(serie)])
pcodec <- pcodec %>%
  mutate(valor_real = valores_reales[as.character(serie)])


#### ----------------- Root Mean Squared Error (RMSE) ----------------------- 
results_rmsePearson <- pp %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()

results_rmseSpearman <- psp %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()

results_rmseKendall <- pk %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()

results_rmseXi <- pxi %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()

results_rmseCodec <- pcodec %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()


meanCodec <- pcodec %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = mean(p1),
    RMSE_p2 = mean(p2),
    RMSE_p3 = mean(p3)
  ) %>%
  ungroup()

#### -------- Mean Absolute Error (MAE) -----------

results_maePearson <- pp %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, valor_real),
    MAE_p2 = mae(p2, valor_real),
    MAE_p3 = mae(p3, valor_real)
  ) %>%
  ungroup()

results_maeSpearman <- psp %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, valor_real),
    MAE_p2 = mae(p2, valor_real),
    MAE_p3 = mae(p3, valor_real)
  ) %>%
  ungroup()

results_maeKendall <- pk %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, valor_real),
    MAE_p2 = mae(p2, valor_real),
    MAE_p3 = mae(p3, valor_real)
  ) %>%
  ungroup()

results_maeXi <- pxi %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, valor_real),
    MAE_p2 = mae(p2, valor_real),
    MAE_p3 = mae(p3, valor_real)
  ) %>%
  ungroup()

results_maeCodec <- pcodec %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, valor_real),
    MAE_p2 = mae(p2, valor_real),
    MAE_p3 = mae(p3, valor_real)
  ) %>%
  ungroup()



#### -------- Median Absolute Deviation (MAD) -------
results_madPearson <- pp %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, valor_real),
    MAD_p2 = mad(p2, valor_real),
    MAD_p3 = mad(p3, valor_real)
  ) %>%
  ungroup()


results_madSpearman <- psp %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, valor_real),
    MAD_p2 = mad(p2, valor_real),
    MAD_p3 = mad(p3, valor_real)
  ) %>%
  ungroup()

results_madKendall <- pk %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, valor_real),
    MAD_p2 = mad(p2, valor_real),
    MAD_p3 = mad(p3, valor_real)
  ) %>%
  ungroup()

results_madXi <- pxi %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, valor_real),
    MAD_p2 = mad(p2, valor_real),
    MAD_p3 = mad(p3, valor_real)
  ) %>%
  ungroup()

results_madCodec <- pcodec %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, valor_real),
    MAD_p2 = mad(p2, valor_real),
    MAD_p3 = mad(p3, valor_real)
  ) %>%
  ungroup()


#### --------- plots of errors -------------------------

## Unifying the rmse error for all coefficients
results_rmseCodec$Coefficient<-"Codec"
results_rmseKendall$Coefficient<-"Kendall"
results_rmsePearson$Coefficient<-"Pearson"
results_rmseSpearman$Coefficient<-"Spearman"
results_rmseXi$Coefficient<-"Xi"
resultsRMSE<-rbind(results_rmsePearson,results_rmseSpearman,results_rmseKendall,results_rmseXi,results_rmseCodec)

## Now, unifying the MAE
results_maeCodec$Coefficient<-"Codec"
results_maeKendall$Coefficient<-"Kendall"
results_maePearson$Coefficient<-"Pearson"
results_maeSpearman$Coefficient<-"Spearman"
results_maeXi$Coefficient<-"Xi"
resultsMAE<-rbind(results_maePearson,results_maeSpearman,results_maeXi,results_maeCodec)

## Finally, unifying the MAD
results_madCodec$Coefficient<-"Codec"
results_madKendall$Coefficient<-"Kendall"
results_madPearson$Coefficient<-"Pearson"
results_madSpearman$Coefficient<-"Spearman"
results_madXi$Coefficient<-"Xi"
resultsMAD<-rbind(results_madPearson,results_madSpearman,results_madXi,results_madCodec)

##Now, unifying the errors
resultsRMSE$Error<-"RMSE"
resultsMAE$Error<-"MAE"
resultsMAD$Error<-"MAD"
colnames(resultsRMSE)<-c("Size","Serie","p1","p2","p3","Coefficient","Error")
colnames(resultsMAE)<-colnames(resultsMAD)<-colnames(resultsRMSE)
Errors<-rbind(resultsRMSE,resultsMAE,resultsMAD)

Errors <- tidyr::pivot_longer(Errors, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")



ggplot(data=Errors[Errors$Error=="RMSE",], aes(x=as.numeric(Size), y=Value,colour = Serie)) +
  #theme_bw()+
  geom_line()+
  facet_grid(Coefficient~Estimator,scales="free_y")+
  xlab("Size")+theme_bw()
theme(legend.position = "none")

