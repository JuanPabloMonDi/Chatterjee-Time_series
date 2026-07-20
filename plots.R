#This file uses the results files save in the Chatterjee_sim/Results folder 

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

#plasso_long<- tidyr::pivot_longer(plasso, cols = c(p1,p2,p3), names_to = "Estimator", values_to = "Value")




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


