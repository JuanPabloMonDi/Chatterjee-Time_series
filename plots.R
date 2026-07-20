#This file uses the results files save in the Chatterjee_sim/Results folder 
library(tidyverse)


#Pearson
pp<-read.csv()
#Spearman
ppart<-read.csv()
# Kendall
psp<-read.csv()
#Chatterjee's Xi
pk
#Azadkia-Chatterjee CODEC
pcodec
#AIC
paic
#BIC
pib
#LASSO
#plasso


pp_long <- tidyr::pivot_longer(pp, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
psp_long <- tidyr::pivot_longer(psp, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
pk_long <- tidyr::pivot_longer(pk, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")

pxi_long <- tidyr::pivot_longer(pxi, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")

pcodec_long <- tidyr::pivot_longer(pcodec, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")


partp_long<- tidyr::pivot_longer(partp, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")
paic_long <- tidyr::pivot_longer(paic, cols = c(p1), names_to = "Estimator", values_to = "Value")
pbic_long <- tidyr::pivot_longer(pbic, cols = c(p1), names_to = "Estimator", values_to = "Value")

paic_long<-rbind(paic_long%>%mutate(metric="AIC"),pbic_long%>%mutate(metric="BIC"))

#plasso_long<- tidyr::pivot_longer(plasso, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")




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
      scales = "free_y",
      labeller = my_labeller  ) +
    theme_linedraw(paper = "white", ink = "gray20") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
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

# Spearman plots ----------------------------------------------------------------
spearman_plot1<-comparison_plot(serie="x1",inter=5, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot2<-comparison_plot(serie="x2",inter=1, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot3<-comparison_plot(serie="x3",inter=2, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot4<-comparison_plot(serie="x4",inter=3, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot5<-comparison_plot(serie="x5",inter=4, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot6<-comparison_plot(serie="x6",inter=2, data=psp_long, estimators = "Estimator",label_est = labs_est)
spearman_plot7<-comparison_plot(serie="x7",inter=2, data=psp_long, estimators = "Estimator",label_est = labs_est)


# Kendall plots ------------------------------------------------------
kendall_plot1<-comparison_plot(serie="x1",inter=5, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot2<-comparison_plot(serie="x2",inter=1, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot3<-comparison_plot(serie="x3",inter=2, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot4<-comparison_plot(serie="x4",inter=3, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot5<-comparison_plot(serie="x5",inter=4, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot6<-comparison_plot(serie="x6",inter=2, data=pk_long, estimators = "Estimator",label_est = labs_est)
kendall_plot7<-comparison_plot(serie="x7",inter=2, data=pk_long, estimators = "Estimator",label_est = labs_est)

# Chatterjee's Xi plots ------------------------------------------------------
xi_plot1<-comparison_plot(serie="x1",inter=5, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot2<-comparison_plot(serie="x2",inter=1, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot3<-comparison_plot(serie="x3",inter=2, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot4<-comparison_plot(serie="x4",inter=3, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot5<-comparison_plot(serie="x5",inter=4, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot6<-comparison_plot(serie="x6",inter=2, data=pxi_long, estimators = "Estimator",label_est = labs_est)
xi_plot7<-comparison_plot(serie="x7",inter=2, data=pxi_long, estimators = "Estimator",label_est = labs_est)

# Azadkia-Chatterjee plots ------------------------------------------------------
codec_plot1<-comparison_plot(serie="x1",inter=5, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot2<-comparison_plot(serie="x2",inter=1, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot3<-comparison_plot(serie="x3",inter=2, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot4<-comparison_plot(serie="x4",inter=3, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot5<-comparison_plot(serie="x5",inter=4, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot6<-comparison_plot(serie="x6",inter=2, data=pcodec_long, estimators = "Estimator",label_est = labs_est)
codec_plot7<-comparison_plot(serie="x7",inter=2, data=pcodec_long, estimators = "Estimator",label_est = labs_est)

# AIC-BIC plots -------------------------------------------------------------------
aic_plot1<-comparison_plot(serie="x1",inter=5, data=paic_long, estimators = "metric")
aic_plot2<-comparison_plot(serie="x2",inter=1, data=paic_long, estimators = "metric")
aic_plot3<-comparison_plot(serie="x3",inter=2, data=paic_long, estimators = "metric")
aic_plot4<-comparison_plot(serie="x4",inter=3, data=paic_long, estimators = "metric")
aic_plot5<-comparison_plot(serie="x5",inter=4, data=paic_long, estimators = "metric")
aic_plot6<-comparison_plot(serie="x6",inter=2, data=paic_long, estimators = "metric")
aic_plot7<-comparison_plot(serie="x7",inter=2, data=paic_long, estimators = "metric")





#The function also works to compare measures on the same plot, as example
estimators<-"name"
data<-rbind(pp_long%>%mutate(name="Pearson"),
            pcodec_long%>%mutate(name="CODEC") 
            #pk_long%>%mutate(name="Kendall"),
            #psp_long%>%mutate(name="Spearman")
            )
comparison_plot("x7","name",data=data, inter=2)
