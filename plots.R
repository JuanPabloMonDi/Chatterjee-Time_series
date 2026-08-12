#This file uses the results files save in the Results folder 
library(tidyverse)
library(Metrics)

# First, we get and unify the simulation results stored in the Results folder

## Pearson coefficient ----------------------------------------------------------------
results_pearsonB<- read_csv("Results_sim/Results_simB/Results_Pearson.csv")
results_pearsonC<- read_csv("Results_sim/Results_simC/Results_Pearson.csv")
results_pearsonC<-results_pearsonC%>%mutate(sim=sim+100)
results_pearson<-rbind(results_pearsonB,results_pearsonC)
rm(results_pearsonB,results_pearsonC)


## Spearman coefficient  ---------------------------------------------------------
results_spearmanB<- read_csv("Results_sim/Results_simB/Results_Spearman.csv")
results_spearmanC<- read_csv("Results_sim/Results_simC/Results_Spearman.csv")
results_spearmanC<-results_spearmanC%>%mutate(sim=sim+100) #this is just to enumerate the simulations
results_spearman<-rbind(results_spearmanB,results_spearmanC)
rm(results_spearmanB,results_spearmanC)

## Kendall Rank coefficient ------------------------------------------------
results_kendallB<- read_csv("Results_sim/Results_simB/Results_Kendall.csv")
results_kendallC<- read_csv("Results_sim/Results_simC/Results_Kendall.csv")
results_kendallC<-results_kendallC%>%mutate(sim=sim+100)
results_kendall<-rbind(results_kendallB,results_kendallC)
rm(results_kendallB,results_kendallC)


## Chatterjee's Xi ----------------------------------------------------
results_xiB<- read_csv("Results_sim/Results_simB/Results_Chatterjee.csv")
results_xiC<- read_csv("Results_sim/Results_simC/Results_Chatterjee.csv")
results_xiC<-results_xiC%>%mutate(sim=sim+100)
results_xi<-rbind(results_xiB,results_xiC)
rm(results_xiB,results_xiC)


## Azadkia-Chatterjee FOCI-CODEC ----------------------------------------------
results_fociB<- read_csv("Results_sim/Results_simB/Results_FOCI.csv")
results_fociC<- read_csv("Results_sim/Results_simC/Results_FOCI.csv")
results_fociC<-results_fociC%>%mutate(sim=sim+100)
results_foci<-rbind(results_fociB,results_fociC)
rm(results_fociB,results_fociC)


## AIC & BIC ----------------------------------------------------------
results_aicB<- read_csv("Results_sim/Results_simB/Results_AIC.csv")
results_aicC<- read_csv("Results_sim/Results_simC/Results_AIC.csv")
results_aicC<-results_aicC%>%mutate(sim=sim+100, metric="AIC")
results_bicB<- read_csv("Results_sim/Results_simB/Results_BIC.csv")
results_bicC<- read_csv("Results_sim/Results_simC/Results_BIC.csv")
results_bicC<-results_bicC%>%mutate(sim=sim+100, metric="BIC")


results_aic<-rbind(results_aicB,results_aicC,
                   results_bicB,results_bicC)
rm(results_aicB,results_aicC,results_bicB,results_bicC)


pp_long     <- tidyr::pivot_longer(results_pearson, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
psp_long    <- tidyr::pivot_longer(results_spearman, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
pk_long     <- tidyr::pivot_longer(results_kendall, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
pxi_long    <- tidyr::pivot_longer(results_xi, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
pcodec_long <- tidyr::pivot_longer(results_foci, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
paic_long   <- tidyr::pivot_longer(results_aic, cols = c(p1), names_to = "Estimator", values_to = "Value")




### ------------------- Plots of the series -----------------------



#Titles for the plots
labs_est <- c(
  p1 = "Max~lag~p[1]",
  p2 = "Second~Max~lag~p[2]",
  p3 = "Third~Max~lag~p[3]"
)

#To avoid repetition of code, we create a function to make the plot of each serie-autocorrelation measure combination
comparison_plot<-function(series,inter,data,estimators="Estimator", factor_plot="size",label_est=NULL){
  #labels for facet wrap boxes
  if (is.null(label_est)) {my_labeller <- "label_value"}else{
    my_labeller <- labeller(!!sym(estimators) := as_labeller(label_est, label_parsed))}
  
  
  plot<-ggplot(data%>%filter(serie==series),aes(x = factor(!!sym(factor_plot)), y = Value)) +
    #Violin plot
    geom_violin(trim = F,   
                fill = "grey60",
                color = "grey40",
                alpha = 0.25,
                lwd=.6) +
    #Boxplot
    geom_boxplot(width = 0.1, outlier.size = 0.9,
                 fill = "gray20",
                 color = "gray20",
                 lwd=.9,
                 median.colour = "white")+
    #Add triangle to indicate the mean of the distribution
    stat_summary(
      fun = mean,
      geom = "point",
      shape = 25,     
      size = 2.5,
      fill = "white",  
      color = "black"
    )+
    #Add line ot indicate the real value of the parameter
    geom_hline(yintercept = inter, colour = "black", lwd = 1, linetype=2, alpha=0.5) +
    labs(x="Size",y="Lag") +
    facet_wrap(
      vars(!!sym(estimators)),
      scales = "fixed",
      labeller = my_labeller  ) +
    theme_linedraw(paper = "white", ink = "gray20") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "white", color = "black"),
      strip.text = element_text(color = "black", size = 10, face = "bold"),
      axis.title.x = element_text(size=12)
    )
  return(plot)
}

# Pearson plots -----------------------------------------------------------------
pearson_plot1<-comparison_plot(serie="x1",inter=5, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot2<-comparison_plot(serie="x2",inter=1, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot3<-comparison_plot(serie="x3",inter=2, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot4<-comparison_plot(serie="x4",inter=3, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot5<-comparison_plot(serie="x5",inter=4, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot6<-comparison_plot(serie="x6",inter=2, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot7<-comparison_plot(serie="x7",inter=2, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot8<-comparison_plot(serie="x8",inter=2, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot9<-comparison_plot(serie="x9",inter=24, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot10<-comparison_plot(serie="x10",inter=24, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot11<-comparison_plot(serie="x11",inter=24, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot12<-comparison_plot(serie="x12",inter=1, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot13<-comparison_plot(serie="x13",inter=3, data=pp_long, estimators = "Estimator",label_est = labs_est)
pearson_plot14<-comparison_plot(serie="x14",inter=3, data=pp_long, estimators = "Estimator",label_est = labs_est)



# Spearman plots ----------------------------------------------------------------
spearman_plot1<-comparison_plot(serie="x1",inter=5, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot2<-comparison_plot(serie="x2",inter=1, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot3<-comparison_plot(serie="x3",inter=2, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot4<-comparison_plot(serie="x4",inter=3, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot5<-comparison_plot(serie="x5",inter=4, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot6<-comparison_plot(serie="x6",inter=2, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot7<-comparison_plot(serie="x7",inter=2, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot8<-comparison_plot(serie="x8",inter=2, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot9<-comparison_plot(serie="x9",inter=24, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot10<-comparison_plot(serie="x10",inter=24, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot11<-comparison_plot(serie="x11",inter=24, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot12<-comparison_plot(serie="x12",inter=1, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot13<-comparison_plot(serie="x13",inter=3, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot14<-comparison_plot(serie="x14",inter=3, data=psp_long, estimators = "Estimator",label_est = labs_est)


# Kendall plots ------------------------------------------------------
kendall_plot1<-comparison_plot(serie="x1",inter=5, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot2<-comparison_plot(serie="x2",inter=1, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot3<-comparison_plot(serie="x3",inter=2, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot4<-comparison_plot(serie="x4",inter=3, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot5<-comparison_plot(serie="x5",inter=4, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot6<-comparison_plot(serie="x6",inter=2, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot7<-comparison_plot(serie="x7",inter=2, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot8<-comparison_plot(serie="x8",inter=2, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot9<-comparison_plot(serie="x9",inter=24,data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot10<-comparison_plot(serie="x10",inter=24, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot11<-comparison_plot(serie="x11",inter=24, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot12<-comparison_plot(serie="x12",inter=1, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot13<-comparison_plot(serie="x13",inter=3, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot14<-comparison_plot(serie="x14",inter=3, data=pk_long, estimators = "Estimator",label_est = labs_est)

# Chatterjee's Xi plots ------------------------------------------------------
xi_plot1<-comparison_plot(serie="x1",inter=5, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot2<-comparison_plot(serie="x2",inter=1, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot3<-comparison_plot(serie="x3",inter=2, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot4<-comparison_plot(serie="x4",inter=3, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot5<-comparison_plot(serie="x5",inter=4, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot6<-comparison_plot(serie="x6",inter=2, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot7<-comparison_plot(serie="x7",inter=2, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot8<-comparison_plot(serie="x8",inter=2, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot9<-comparison_plot(serie="x9",inter=24,data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot10<-comparison_plot(serie="x10",inter=24, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot11<-comparison_plot(serie="x11",inter=24, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot12<-comparison_plot(serie="x12",inter=1, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot13<-comparison_plot(serie="x13",inter=3, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot14<-comparison_plot(serie="x14",inter=3, data=pxi_long, estimators = "Estimator",label_est = labs_est)


# Azadkia-Chatterjee plots ------------------------------------------------------
codec_plot1<-comparison_plot(serie="x1",inter=5, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot2<-comparison_plot(serie="x2",inter=1, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot3<-comparison_plot(serie="x3",inter=2, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot4<-comparison_plot(serie="x4",inter=3, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot5<-comparison_plot(serie="x5",inter=4, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot6<-comparison_plot(serie="x6",inter=2, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot7<-comparison_plot(serie="x7",inter=2, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot8<-comparison_plot(serie="x8",inter=2, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot9<-comparison_plot(serie="x9",inter=24, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot10<-comparison_plot(serie="x10",inter=24, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot11<-comparison_plot(serie="x11",inter=24, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot12<-comparison_plot(serie="x12",inter=1, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot13<-comparison_plot(serie="x13",inter=3, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot14<-comparison_plot(serie="x14",inter=3, data=pcodec_long, estimators = "Estimator",label_est = labs_est)

# AIC-BIC plots -------------------------------------------------------------------
aic_plot1<-comparison_plot(serie="x1",inter=5, data=paic_long, estimators = "metric")
aic_plot2<-comparison_plot(serie="x2",inter=1, data=paic_long, estimators = "metric")
aic_plot3<-comparison_plot(serie="x3",inter=2, data=paic_long, estimators = "metric")
aic_plot4<-comparison_plot(serie="x4",inter=3, data=paic_long, estimators = "metric")
aic_plot5<-comparison_plot(serie="x5",inter=4, data=paic_long, estimators = "metric")
aic_plot6<-comparison_plot(serie="x6",inter=2, data=paic_long, estimators = "metric")
aic_plot7<-comparison_plot(serie="x7",inter=2, data=paic_long, estimators = "metric")
aic_plot8<-comparison_plot(serie="x8",inter=2, data=paic_long, estimators = "metric")
aic_plot9<-comparison_plot(serie="x9",inter=4, data=paic_long, estimators = "metric")
aic_plot10<-comparison_plot(serie="x10",inter=4, data=paic_long, estimators = "metric")
aic_plot11<-comparison_plot(serie="x11",inter=4, data=paic_long, estimators = "metric")
aic_plot12<-comparison_plot(serie="x12",inter=1, data=paic_long, estimators = "metric")
aic_plot13<-comparison_plot(serie="x13",inter=3, data=paic_long, estimators = "metric")
aic_plot14<-comparison_plot(serie="x14",inter=3, data=paic_long, estimators = "metric")





#The function also works to compare measures on the same plot, as example
estimators<-"name"
data<-rbind(pp_long%>%mutate(name="Pearson"),
            pcodec_long%>%mutate(name="CODEC") 
            #pk_long%>%mutate(name="Kendall"),
            #psp_long%>%mutate(name="Spearman")
            )
comparison_plot("x7","name",data=data, inter=2)



# Saving the plots

#to avoid the repetition of lines, we are going to create a function to save all the plots by its respective measure

# Put all plots into named lists
pearson_plots <- mget(ls(pattern = "^pearson_plot[0-9]+$"))
spearman_plots <- mget(ls(pattern = "^spearman_plot[0-9]+$"))
kendall_plots  <- mget(ls(pattern = "^kendall_plot[0-9]+$"))
xi_plots       <- mget(ls(pattern = "^xi_plot[0-9]+$"))
codec_plots    <- mget(ls(pattern = "^codec_plot[0-9]+$"))
aic_plots      <- mget(ls(pattern = "^aic_plot[0-9]+$"))


# Function to save a list of plots
save_plots <- function(plot_list, folder, width = 10, height = 6, dpi = 600,units="cm") {
  
  # Create folder if it doesn't exist
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE)
  }
  
  # Save each plot in the respective folder 
  for (plot_name in names(plot_list)) {
    
    ggsave(
      filename = file.path(folder, paste0(plot_name, ".pdf")),
      plot = plot_list[[plot_name]],
      width = width,
      height = height,
      dpi = dpi,
      units = units
    )
  }
}


# Save each group in its own folder
save_plots(pearson_plots,  "plots/Pearson", units = "in")
save_plots(spearman_plots, "plots/Spearman", units = "in")
save_plots(kendall_plots,  "plots/Kendall", units = "in")
save_plots(xi_plots,       "plots/Xi", units = "in")
save_plots(codec_plots,    "plots/Azadkia-Chatterjee", units = "in")
save_plots(aic_plots,      "plots/AIC-BIC", units = "in")

