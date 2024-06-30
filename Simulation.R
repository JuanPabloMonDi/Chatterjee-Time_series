### Load required packages ------------------------

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
source("Functions.R")
### set seed for replicability -------------------
set.seed(131101)


### Simulate time series -----------------------

x1<-arima.sim(list(2,1,1),n=1000)
x2<-arima.sim(list(3,0,0),n=1000)



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




#Fixed terms
n_sim=200
lag.max=30

pp=NULL
psp=NULL
pk=NULL
pxi=NULL
pcodec<-NULL

#autop<-NULL
#autosp<-NULL
#autok<-NULL
#autoxi<-NULL
#autocodec<-NULL
sizes<-rev(c(50,100,500,1000,2000,5000,10000))
n_sim<-250
n_size=0
for (size in sizes){
  n_size=n_size+1
#Simulations
for (sim in 1:n_sim){
  
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
      corP<-pcor(data,method="pearson")
      corSp<-pcor(data,method="spearman")
      corK<-pcor(data,method="kendall")
      corX<-partcor(data,method = "xi")
  
      
      #Now, we calculate the autocorrelation  matrix for each correlation
      
      #autocorP<-cor(data,method="pearson")
      #autocorSp<-cor(data,method="spearman")
      #autocorK<-cor(data,method="kendall")
      #autocorX<-xicor(data)
      
      #Since we are only interested in the correlation with xt, we only consider the first row
      corP<-corP$estimate[1,]
      corSp<-corSp$estimate[1,]
      corK<-corK$estimate[1,]
      corX<-corX[1,]
      #autocorP<-autocorP$estimate[1,]
      #autocorSp<-autocorSp$estimate[1,]
      #autocorK<-autocorK$estimate[1,]
      #autocorX<-aut
      
      
      corPearson<-append(corPearson, corP[lag+1])
      corSpearman<-append(corSpearman,corSp[lag+1])
      corKendall<-append(corKendall,corK[lag+1])
      corXi<-append(corXi,corX[lag+1])
      #print(paste("|--lag ",lag,"/",lag.max ,"|--serie",serie, "|--simulation ",sim,"/",n_sim," |-- size",size))
      }
      #Confidence intervals
      #For Pearson the variance is 1/N
      Plags<-data.frame(lag=(index(corPearson)-1),value=corPearson)
      #plot(corPearson[1:lag],type="h")
      #abline(h=qnorm(0.975)*sqrt(1/size))
      #abline(h=-qnorm(0.975)*sqrt(1/size))
      Plags<-Plags[abs(Plags$value)>qnorm(0.975)*sqrt(1/size),]
      Plags<-Plags$lag
      npp<-length(Plags)
      pp<-rbind(pp,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(npp>=1, Plags[[length(Plags)]], 0),p2=ifelse(npp>=2, Plags[(length(Plags)-1)], 0),p3=ifelse(npp>=3, Plags[(length(Plags)-2)], 0)))

  
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
      Klags<-data.frame(lag=(index(corKendall)-1),value=corKendall)
      Klags<-Klags[abs(Klags$value)>(qnorm(0.975)*sqrt(2*(2*size+5)/(9*size*(size-1)))),]
      Klags<-Klags$lag
      nk<-length(Klags)
      pk<-rbind(pk,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(nk>=1, Klags[[length(Klags)]], 0),p2=ifelse(nk>=2, Klags[[(length(Klags)-1)]], 0),p3=ifelse(nk>=3, Klags[[(length(Klags)-2)]], 0)))
    
      #For Xi the variance is 2/5N
      plot(corXi[1:lag],type="h")
      abline(h=qnorm(0.975)*sqrt(2/(5*size)),col="red",lw=1)
      abline(h=-qnorm(0.975)*sqrt(2/(5*size)),col="red",lw=1)
      #Get the 3 maximum autocorrelations
      Xilags<-data.frame(lag=(index(corXi)-1),value=corXi)
      Xilags<-Xilags[abs(Xilags$value)>(qnorm(0.975)*sqrt(2/(5*size))),]
      Xilags<-Xilags$lag#[(nrow(Xilags)-2):nrow(Xilags)]
      nxi=length(Xilags)
      pxi<-rbind(pxi,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(nxi>=1, Xilags[[length(Xilags)]], 0),p2=ifelse(nxi>=2, Xilags[[(length(Xilags)-1)]], 0),p3=ifelse(nxi>=3, Xilags[[(length(Xilags)-2)]], 0)))

      #codecM<-append(codecM,codec(data[,1],data[,ncol(data)],data[,2:(ncol(data)-1)]))
        
      #Now, with the coefficient of Azadka-Chatterjee
      codecM<-foci(data$xt,data[2:(lag+1)],numCores = 1)
      codecM<-codecM$selectedVar$index
      ncodec=length(codecM)
      pcodec<-rbind(pcodec,data.frame(size=size,sim=sim, serie=serie,p1=ifelse(ncodec>=1, max(codecM), 0),p2=ifelse(ncodec>=2, sort(codecM, TRUE)[2], 0),p3=ifelse(ncodec>=3,sort(codecM, TRUE)[3], 0)))
  } #End series iterations
} #End simulations
}#End sizes 

cat("\n")
#Iteraciones orden size>n_sim>series>lag



pp_long <- tidyr::pivot_longer(pp, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
psp_long <- tidyr::pivot_longer(psp, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
pk_long <- tidyr::pivot_longer(pk, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")

pxi_long <- tidyr::pivot_longer(pxi, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")

pcodec_long <- tidyr::pivot_longer(pcodec, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")

#Plots of the series


#Serie 1

ggplot(pp_long[pp_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ 
  geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x1",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 3,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

#Serie 2

ggplot(pp_long[pp_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+
  facet_wrap(~Estimator,scales = "free_y")+
  theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+
  facet_wrap(~Estimator,scales = "free_y")+
  theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x2",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 4,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))

#Serie 3
ggplot(pp_long[pp_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x3",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  #geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))


#Serie 4

ggplot(pp_long[pp_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x4",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))





#Serie 5

ggplot(pp_long[pp_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")


ggplot(psp_long[psp_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pk_long[pk_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")



ggplot(pxi_long[pxi_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.5)+ geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+facet_wrap(~Estimator,scales = "free_y")+theme_bw()+theme(legend.position = "none")

ggplot(pcodec_long[pcodec_long$serie=="x5",], aes(x=factor(size), y=Value,fill = factor(size))) + 
  geom_violin(trim = F)+
  geom_boxplot(width=0.3)+
  geom_hline(yintercept = 2,colour="blue",lwd=1.25)+
  xlab("Size")+
  theme_bw()+theme(legend.position = "none")+
  facet_wrap(~Estimator,scales = "free_y",labeller = labeller(variable=c(p1="1st Max",p2="2nd Max",p3="3rd Max")))



  
  
  #rmse
# Valores reales para las series
valores_reales <- c(3, 4, 2, 2, 2)
names(valores_reales) <- c("x1","x2","x3","x4","x5")


# Función para calcular el RMSE
rmse <- function(estimados, reales) {
  sqrt(mean((estimados - reales) ^ 2,na.rm=T))
}

# Añadir los valores reales a los datos
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


# Calcular RMSE para cada combinación de tamaño y serie
resultados_rmsePearson <- pp %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()

resultados_rmseSpearman <- psp %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()

resultados_rmseKendall <- pk %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()

resultados_rmseXi <- pxi %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()

resultados_rmseCodec <- pcodec %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, valor_real),
    RMSE_p2 = rmse(p2, valor_real),
    RMSE_p3 = rmse(p3, valor_real)
  ) %>%
  ungroup()

