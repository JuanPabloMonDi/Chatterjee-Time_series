### Load required packages ------------------------
library(tidyr)
library(tidyverse)
library(Metrics)

### ------------------------------------------------------
### Errors Code ------------------------------------------
### ------------------------------------------------------

# First, we get and unify the simulation results stored in the Results folder

## Association measures ---------------------------------

#### Pearson coefficient ----------------------------------------------------------------
results_pearsonB<- read_csv("Results_sim/Results_simB/Results_Pearson.csv")
results_pearsonC<- read_csv("Results_sim/Results_simC/Results_Pearson.csv")
results_pearsonC<-results_pearsonC%>%mutate(sim=sim+100)
results_pearson<-rbind(results_pearsonB,results_pearsonC)
rm(results_pearsonB,results_pearsonC)


#### Spearman coefficient  ---------------------------------------------------------
results_spearmanB<- read_csv("Results_sim/Results_simB/Results_Spearman.csv")
results_spearmanC<- read_csv("Results_sim/Results_simC/Results_Spearman.csv")
results_spearmanC<-results_spearmanC%>%mutate(sim=sim+100) #this is just to enumerate the simulations
results_spearman<-rbind(results_spearmanB,results_spearmanC)
rm(results_spearmanB,results_spearmanC)

#### Kendall Rank coefficient ------------------------------------------------
results_kendallB<- read_csv("Results_sim/Results_simB/Results_Kendall.csv")
results_kendallC<- read_csv("Results_sim/Results_simC/Results_Kendall.csv")
results_kendallC<-results_kendallC%>%mutate(sim=sim+100)
results_kendall<-rbind(results_kendallB,results_kendallC)
rm(results_kendallB,results_kendallC)


#### Chatterjee's Xi ----------------------------------------------------
results_xiB<- read_csv("Results_sim/Results_simB/Results_Chatterjee.csv")
results_xiC<- read_csv("Results_sim/Results_simC/Results_Chatterjee.csv")
results_xiC<-results_xiC%>%mutate(sim=sim+100)
results_xi<-rbind(results_xiB,results_xiC)
rm(results_xiB,results_xiC)


#### Azadkia-Chatterjee FOCI-CODEC ----------------------------------------------
results_fociB<- read_csv("Results_sim/Results_simB/Results_FOCI.csv")
results_fociC<- read_csv("Results_sim/Results_simC/Results_FOCI.csv")
results_fociC<-results_fociC%>%mutate(sim=sim+100)
results_foci<-rbind(results_fociB,results_fociC)
rm(results_fociB,results_fociC)


#### AIC & BIC ----------------------------------------------------------
results_aicB<- read_csv("Results_sim/Results_simB/Results_AIC.csv")
results_aicB<-results_aicB%>%mutate(metric="AIC")
results_aicC<- read_csv("Results_sim/Results_simC/Results_AIC.csv")
results_aicC<-results_aicC%>%mutate(sim=sim+100, metric="AIC")
results_bicB<- read_csv("Results_sim/Results_simB/Results_BIC.csv")
results_bicB<-results_bicB%>%mutate(metric="BIC")
results_bicC<- read_csv("Results_sim/Results_simC/Results_BIC.csv")
results_bicC<-results_bicC%>%mutate(sim=sim+100, metric="BIC")


results_aic<-rbind(results_aicB,results_aicC,
                   results_bicB,results_bicC)
rm(results_aicB,results_aicC,results_bicB,results_bicC)



## -----Meausing the error ----------------

# Real value of the parameters
real_values <- c(5,1,2,3, 4, 2, 2,2,24,24,24,1,3,3)
names(real_values) <- c("x1", "x2",  "x3", "x4",  "x5", "x6",
                        "x7", "x8", "x9",  "x10", "x11", "x12",
                        "x13", "x14")




# Add the real parameter values to the table
results_pearson <- results_pearson %>%
  mutate(sim_value = real_values[as.character(serie)])
results_spearman <- results_spearman %>%
  mutate(sim_value = real_values[as.character(serie)])
results_kendall <- results_kendall %>%
  mutate(sim_value = real_values[as.character(serie)])
results_xi <- results_xi %>%
  mutate(sim_value = real_values[as.character(serie)])
results_foci <- results_foci %>%
  mutate(sim_value = real_values[as.character(serie)])
results_aic <- results_aic %>%
  mutate(sim_value = real_values[as.character(serie)])


#### ----------------- Root Mean Squared Error (RMSE) ----------------------- 
results_rmsePearson <- results_pearson %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, sim_value),
    RMSE_p2 = rmse(p2, sim_value),
    RMSE_p3 = rmse(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_rmseSpearman <- results_spearman %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, sim_value),
    RMSE_p2 = rmse(p2, sim_value),
    RMSE_p3 = rmse(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_rmseKendall <- results_kendall %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, sim_value),
    RMSE_p2 = rmse(p2, sim_value),
    RMSE_p3 = rmse(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_rmseXi <- results_xi %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, sim_value),
    RMSE_p2 = rmse(p2, sim_value),
    RMSE_p3 = rmse(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_rmseCodec <- results_foci %>%
  group_by(size, serie) %>%
  summarise(
    RMSE_p1 = rmse(p1, sim_value),
    RMSE_p2 = rmse(p2, sim_value),
    RMSE_p3 = rmse(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_rmseAIC<- results_aic %>%
  group_by(size, serie,metric) %>%
  summarise(
    RMSE_p1 = rmse(p1, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

#### -------- Mean Absolute Error (MAE) -----------

results_maePearson <- results_pearson %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, sim_value),
    MAE_p2 = mae(p2, sim_value),
    MAE_p3 = mae(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_maeSpearman <- results_spearman %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, sim_value),
    MAE_p2 = mae(p2, sim_value),
    MAE_p3 = mae(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_maeKendall <- results_kendall %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, sim_value),
    MAE_p2 = mae(p2, sim_value),
    MAE_p3 = mae(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_maeXi <- results_xi %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, sim_value),
    MAE_p2 = mae(p2, sim_value),
    MAE_p3 = mae(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_maeCodec <- results_foci %>%
  group_by(size, serie) %>%
  summarise(
    MAE_p1 = mae(p1, sim_value),
    MAE_p2 = mae(p2, sim_value),
    MAE_p3 = mae(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_maeAIC <- results_aic %>%
  group_by(size, serie,metric) %>%
  summarise(
    MAE_p1 = mae(p1, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()


#### -------- Median Absolute Deviation (MAD) -------
results_madPearson <- results_pearson %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, sim_value),
    MAD_p2 = mad(p2, sim_value),
    MAD_p3 = mad(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()


results_madSpearman <- results_spearman %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, sim_value),
    MAD_p2 = mad(p2, sim_value),
    MAD_p3 = mad(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_madKendall <- results_kendall %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, sim_value),
    MAD_p2 = mad(p2, sim_value),
    MAD_p3 = mad(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_madXi <- results_xi %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, sim_value),
    MAD_p2 = mad(p2, sim_value),
    MAD_p3 = mad(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

results_madCodec <- results_foci %>%
  group_by(size, serie) %>%
  summarise(
    MAD_p1 = mad(p1, sim_value),
    MAD_p2 = mad(p2, sim_value),
    MAD_p3 = mad(p3, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()


results_madAIC <- results_aic %>%
  group_by(size, serie,metric) %>%
  summarise(
    MAD_p1 = mad(p1, sim_value),
    .groups = "drop"
  ) %>%
  ungroup()

### --------- Plots of errors -------------------------

## Unifying the rmse error for all coefficients
results_rmseCodec$Coefficient<-"Codec"
results_rmseKendall$Coefficient<-"Kendall"
results_rmsePearson$Coefficient<-"Pearson"
results_rmseSpearman$Coefficient<-"Spearman"
results_rmseXi$Coefficient<-"Xi"
resultsRMSE<-rbind(results_rmsePearson,results_rmseSpearman,
                   results_rmseKendall,results_rmseXi,
                   results_rmseCodec)

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

# Now the AIC
results_rmseAIC<-results_rmseAIC%>%rename(Coefficient=metric,Value=RMSE_p1)%>%mutate(Error="RMSE",Estimator="p1")
results_maeAIC<-results_maeAIC%>%rename(Coefficient=metric,Value=MAE_p1)%>%mutate(Error="MAE", Estimator="p1")
results_madAIC<-results_madAIC%>%rename(Coefficient=metric,Value=MAD_p1)%>%mutate(Error="MAD", Estimator="p1")
errors_AIC<-rbind(results_rmseAIC,results_maeAIC,results_madAIC)
colnames(errors_AIC)[c(1,2)]<-c("Size","Serie")
Errors<-rbind(Errors,errors_AIC)



Errors<-Errors%>%mutate(name_serie=case_when(
  Serie=="x1"~"AR(5)",
  Serie=="x2"~"ARMA(1,3)",
  Serie=="x3"~"ARIMA(2,1,1)",
  Serie=="x4"~"NLAR(3)",
  Serie=="x5"~"NLAR(4)",
  Serie=="x6"~"NLARMA(2,2)",
  Serie=="x7"~"SETAR(2,2,[2])",
  Serie=="x8"~"ARIMA(2,1,1) \nDifferenciated",
  Serie=="x9"~"SARIMA(4,0,0) × (2,0,0)[12]",
  Serie=="x10"~"SARIMA(4,0,0) × (2,0,0)[12] \nAdditive decomposition",
  Serie=="x11"~"SARIMA(4,0,0) × (2,0,0)[12] \n Multiplicative decomposition",
  Serie=="x12"~"ARMA(1,0,1)-GARCH(1,1)",
  Serie=="x13"~"ARIMA(3,1,0)",
  Serie=="x14"~"ARIMA(3,1,0) \nDifferenciated",
))


error_plot<-function(error="RMSE",estimator="p1",measures=c("AIC","BIC","Codec","Pearson"),data=Errors,geometry="line", save=TRUE){
  
  k <- as.numeric(sub("^p", "", estimator))
  plot<-ggplot(data=data%>%filter(Error==error,
                            Estimator==estimator,
                            Coefficient%in%measures
                            ),
               aes(
                 x = as.numeric(Size),
                 y = Value,
                 colour = Coefficient,
                 linetype = Coefficient,
                 linewidth = Coefficient
                 )) +
    geom_point(size=1.2)
  if(geometry=="line"){plot<-plot+geom_line()}
  if(geometry=="smooth"){plot<-plot+geom_smooth(se=FALSE)}
  plot<-plot+scale_colour_manual(
    values = c(
      "Codec" = "black",
      "AIC" = "gray30",
      "BIC" = "gray50",
      "Pearson" = "gray20"
      )
    )+
    scale_linetype_manual(
      values = c(
        "Codec" = "solid",
        "AIC" = "dashed",
        "BIC" = "dotdash",
        "Pearson" = "dotted"
        )
      ) +
  
  scale_linewidth_manual(
    values = c(
      "Codec" = 1.2,
      "AIC" = 0.8,
      "BIC" = 0.8,
      "Pearson" = 0.8
    ),
    guide = "none"
  ) +
  facet_wrap(~name_serie, 
             nrow=2)+
  labs(x="Size",
       y=bquote(hat(italic(p))[.(k)]),
       colour="Measure",
       linetype="Measure")+
  theme_linedraw(paper = "white", ink = "gray20") +
  theme(
    panel.grid.major = element_blank(), #Panel of the background
    panel.grid.minor = element_blank(), #Panel inside each serie plot
    strip.background = element_rect(fill = "white", color = "black"), #Back ground of each serie title
    strip.text = element_text(color = "black", size = 8.5, face = "bold"), #Title of each serie
    axis.title.x = element_text(size=11, face = "bold"), # hat p axis
    axis.title.y = element_text(size=11, face = "bold"), #"Size" Axis
    legend.position = "bottom",
    legend.text = element_text(size = 11), #AIC, BIC, Codec and Pearson text
    legend.title = element_text(size = 12, face = "bold"), #"Measure" text
    legend.key.width = unit(1.5, "cm") #Width of legend icon
  )
  if(save==TRUE){
    ggsave(filename=paste0("plots/errors/",error,"_",estimator,"_",geometry,".pdf"),
           plot=plot,units = "cm",
           width = 42,height = 18, dpi=600)
  }
  
  return(plot)
}


#Now, for each error and estimator, we save the plot of errors
for(er in unique(Errors$Error)){
  for(param in unique(Errors$Estimator)){
    if(param!="p1"){assoc<-c("Codec","Pearson")}else{
      assoc<-c("AIC","BIC","Codec","Pearson")
    }
    for(geom in c("line","smooth")){
      error_plot(error=er,estimator = param,measures = assoc,geometry = geom)
    }
  }
}

