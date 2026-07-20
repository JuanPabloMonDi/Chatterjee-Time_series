#Lets import the necessary packages

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


### Create symmetric XI coefficient -------------------- 

symxi<-function(X,Y){max(xicor(X,Y),xicor(Y,X))}

### Creating the function of autocorrelation ------------------

autocor<-function(SerieX, lag.max=30, method="pearson",plot=TRUE){
  if (method=="codec"){
    autocorrelacionesmanual<-c()
    for (i in 0:lag.max){
      rho<-codec(SerieX[1:(length(SerieX)-i)],SerieX[(i+1):length(SerieX)])
      autocorrelacionesmanual<-append(autocorrelacionesmanual,rho)
    }
    autocorrelacionesmanual<-data.frame(autocorrelacionesmanual)
  }
  if(method=="xi"){
    autocorrelacionesmanual<-c()
    for (i in 0:lag.max){
      rho<-xicor(SerieX[1:(length(SerieX)-i)],SerieX[(i+1):length(SerieX)])
      autocorrelacionesmanual<-append(autocorrelacionesmanual,rho)
    }
    autocorrelacionesmanual<-data.frame(autocorrelacionesmanual)
  }
  if(method=="symxi"){
    autocorrelacionesmanual<-c()
    for (i in 0:lag.max){
      rho<-symxi(SerieX[1:(length(SerieX)-i)],SerieX[(i+1):length(SerieX)])
      autocorrelacionesmanual<-append(autocorrelacionesmanual,rho)
    }
    autocorrelacionesmanual<-data.frame(autocorrelacionesmanual)
  }
  if(method=="kendall"){
    autocorrelacionesmanual<-c()
    for (i in 0:lag.max){
      rho<-cor(SerieX[1:(length(SerieX)-i)],SerieX[(i+1):length(SerieX)],method="kendall")
      autocorrelacionesmanual<-append(autocorrelacionesmanual,rho)
    }
    autocorrelacionesmanual<-data.frame(autocorrelacionesmanual)
  }else{
    if (method=="pearson"| method=="spearman"){
      if(method=="spearman"){SerieX<-rank(SerieX)}
      media<-mean(SerieX)
      autocorrelacionesmanual<-c()
      
      for (i in 0:lag.max){
        suma=0
        for (j in (i+1):length(SerieX)){
          suma<-suma+(SerieX[j]-media)*(SerieX[j-i]-media)}
        rho<-suma/sum((SerieX-mean(SerieX))^2)
        autocorrelacionesmanual<-append(autocorrelacionesmanual,rho)}
      autocorrelacionesmanual<-data.frame(autocorrelacionesmanual)
    }}
  if (plot==TRUE){
    plot(0:lag.max,autocorrelacionesmanual$autocorrelacionesmanual,
         main=paste("ACF Time Series -",method),
         ylab="ACF",
         xlab="lag",type="n")+
      abline(h=0,col="black")+
      segments(0:(nrow(autocorrelacionesmanual)-1),
               rep(0,nrow(autocorrelacionesmanual)),
               0:(nrow(autocorrelacionesmanual)-1), 
               autocorrelacionesmanual$autocorrelacionesmanual, col="black")
  }
  return(autocorrelacionesmanual)
}


### Creating a function for partial correlation


partcor<-function(datos, method="kendall"){
  if (method=="xi"){a<-matrix(0,ncol=ncol(datos),nrow=ncol(datos))
  for (i in 1:ncol(datos)){
    for (j in 1:ncol(datos)){
      a[i,j]<-XICOR::xicor(datos[i],datos[j])
    }
  }}else{a<-cor(datos,method=method)}
  
  b<-solve(a)
  c<-matrix(0,ncol = ncol(a), nrow = ncol(a))
  for (i in 1:ncol(a)){
    for (j in 1:ncol(a)){
      c[i,j]=-b[i,j]/(sqrt(b[i,i])*sqrt(b[j,j]))
      if (i==j){c[i,j]=-c[i,j]}
    }
  }
  return(c)}



### Creating a function for partial autocorrelation --------------
pautocor<- function(Serie,lag.max=30,method="pearson", plot=TRUE){
  if (method=="codec"){
    codecM<-c()
    for (lag in 1:lag.max){
      data<-data.frame()
      for (i in 0:(length(Serie)-lag-1)){
        xt<-Serie[(length(Serie)-i)]
        xt1<-t(data.frame(rev(Serie[(length(Serie)-(lag)-i):(length(Serie)-(i+1))])))
        data<-rbind(data,cbind(xt,xt1))
      }
      row.names(data)<-0:(nrow(data)-1)
      VarReg<-paste0("x",1:(lag))
      colnames(data)<-c("xt",VarReg)
      
      codecM<-append(codecM,codec(data[,1],data[,ncol(data)],data[,2:(ncol(data)-1)]))
      
    }
    codecM<-c(1,codecM)
    if (plot==T){plot(0:(length(codecM)-1),codecM,type="h")}
    return(codecM)
  }else{
    matrixseries<-drop_na(data.frame(h0=rev(Serie[lag.max+1:length(Serie)])))
    autocorrelaciones<-autocor(Serie,lag.max+1,method=method,plot=F)[["autocorrelacionesmanual"]]
    columna0<-autocorrelaciones[1:length(autocorrelaciones)-1]
    autocorrelaciones<-autocorrelaciones[2:length(autocorrelaciones)]
    for (i in 1:lag.max){
      columna<-drop_na(data.frame(rev(Serie[(lag.max+1-i):(length(Serie)-i)])))
      colnames(columna)<-as.character(i)
      matrixseries<-cbind(matrixseries,columna)
      columna2<-c(rev(columna0[1:i+1]), columna0[1:(length(columna0) - (i))])
      
      columna0<-cbind(columna0,columna2)
    }
    partial<-(solve(columna0)%*%autocorrelaciones)
    rownames(partial)<-NULL
    return(partial)}
}



# Function to show the progress bar of the simulations
update_progress <- function(current, total) {
  message <- paste("|--", current, "/", total, "--|")
  cat("\r", message)
  flush.console()
}

# Function to calculate the RMSE
rmse <- function(estimados, reales) {
  sqrt(mean((estimados - reales) ^ 2,na.rm=T))
}

#Function to calculate the MAE
mae<- function(estimate,real_values){
  mean(abs(estimate-real_values),na.rm=T)
}

#Function to calculate the MAD (Median Absolute Deviation)

mad<-function(estimate,real_values){
  median(abs(estimate-real_values))
}

