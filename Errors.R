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
names(real_values) <- c("x1",
                           "x2",
                           "x3",
                           "x4",
                           "x5",
                           "x6",
                           "x7",
                           "x8",
                           "x9",
                           "x10",
                           "x11",
                           "x12",
                           "x13",
                           "x14")




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
    MAE_p = mae(p1, sim_value),
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
  facet_grid(Coefficient~Estimator,scales="fixed")+
  xlab("Size")+theme_bw()
theme(legend.position = "none")

