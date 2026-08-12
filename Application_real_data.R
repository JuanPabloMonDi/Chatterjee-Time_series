library(ppcor)
library(mgcv)
library(gplm)
library(FOCI) #Codec
library(dplyr)
library(tidyr)
library(ggplot2)
library(XICOR)
library(forecast)
library(astsa)
library(rugarch)
library(minerva)
library(tsDyn)
library(TSA)
library(datasets)

library(patchwork)

set.seed(1234)

#Functions to apply FOCI
foci_path <- function(y, X, stop_when_negative = FALSE) {
  if (!is.matrix(X)) X <- as.matrix(X)
  
  p <- ncol(X)
  selected_vars <- c()
  remaining_vars <- 1:p
  codec_path <- data.frame(step = integer(), variable = integer(), codec = numeric())
  
  step <- 1
  repeat {
    if (length(remaining_vars) == 0) break
    
    # Compute marginal codec gain for each remaining variable
    codec_vals <- sapply(remaining_vars, function(j) {
      ifelse(length(selected_vars)==0,codec(Y=y,Z=X[, j, drop = FALSE]),codec(Y=y,X=X[, selected_vars, drop = FALSE],Z=X[, j, drop = FALSE]))
    })
    
    # Find the best variable
    best_idx <- which.max(codec_vals)
    best_var <- remaining_vars[best_idx]
    best_val <- codec_vals[best_idx]
    
    # Optionally stop if CODEC gain is not positive
    if (stop_when_negative && best_val < 0) break
    
    # Store the step result
    codec_path <- rbind(codec_path, data.frame(
      step = step,
      variable = best_var,
      codec = best_val
    ))
    
    # Update selected and remaining variables
    selected_vars <- c(selected_vars, best_var)
    remaining_vars <- setdiff(remaining_vars, best_var)
    step <- step + 1
  }
  
  return(codec_path)
}

chatterjee<-function(serie, method="foci"){
  size<-max(length(serie),nrow(serie))
  lag.max=floor(12*(size/100)^(1/4))
  x<-as.vector(serie)
  corPearson<-c()
  corSpearman<-c()
  corCodec<-c()
  for (lag in 1:lag.max){
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
  }
  
  #Now, with the coefficient of Azadka-Chatterjee
  if(method=="foci"){
    codecM<-foci_path(data$xt,data[2:(lag+1)])
    codecM2<-foci(data$xt,data[2:(lag+1)],numCores = 1)
    codecM2<-codecM2$selectedVar$index
    ncodec=length(codecM2)
    p1=ifelse(ncodec>=1, max(codecM2), 0)
    p2=ifelse(ncodec>=2, sort(codecM2, TRUE)[2], 0)
    p3=ifelse(ncodec>=3,sort(codecM2, TRUE)[3], 0)
    correlation=corCodec
  }
  return(list("p1"=p1,"p2"=p2,"p3"=p3,"foci"=codecM2,"foci_path"=codecM))
}

#Implement the function in the datasets

passengers<-datasets::AirPassengers
lynx<-datasets::lynx
sunspot<-datasets::sunspot.year

url <- "https://raw.githubusercontent.com/jbrownlee/Datasets/master/daily-min-temperatures.csv"
melbourne_temps <- read.csv(url, stringsAsFactors = FALSE)
head(melbourne_temps)
melbourne_clean <- melbourne_temps %>%
  mutate(
    Date = as.Date(Date, format="%Y-%m-%d"),
    Temp = as.numeric(Temp)
  )
melbourne_ts <- ts(melbourne_clean$Temp, start = c(1981, 1), frequency = 365)


#####ANALYSIS FOR THE SUNSPOT COUNT#####

series<-sunspot
chatPasse<-chatterjee(series)
obj<-chatPasse$foci_path
df <- data.frame(lag = as.factor(obj$variable),correlation =obj$codec,step=obj$step)
df$lag <- factor(df$lag, levels = df$lag[order(df$step)])
df$step <- factor(df$step)
df$color <- ifelse(df$correlation < 0, "red", "steelblue")

codec_plot <- ggplot(df, aes(x = step, y = correlation)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_segment(aes(xend = step, yend = 0, color = color), linewidth = 1.2) +
  scale_color_identity() +  # use actual color values
  geom_text(aes(label = lag), vjust = ifelse(df$correlation >= 0, -0.5, 1.5), size = 4, color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal() +
  labs(
    title = "Selected lag in each iteration \n(CODEC value prop.to bar length)",
    x = "Iteration",
    y = ""
  )

df_ts <- data.frame(
  time = as.Date(time(series)),
  values = as.numeric(series)
)

# Time series plot
ts_plot <- ggplot(df_ts, aes(x = time, y = values)) +
  geom_area(fill = "#FFDBBB", alpha = 0.5) +  # filled area under the line
  geom_line(color = "orange", linewidth = 1.2) +  # smooth line
  geom_point(color = "black", size = 1) +  # dot at each observation
  theme(plot.title = element_text(hjust = 1)) +
  theme_minimal() +
  labs(
    title = "Original time series",
    x = "",
    y = "Sunspot count"
  )


combined_plot <- ts_plot + codec_plot +
  plot_layout(ncol = 2, widths = c(1.3, 1)) +
  plot_annotation(
    title = "Analysis of Annual Count of Sunspots",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
x11()
combined_plot


#####ANALYSIS FOR THE LYNX COUNT#####

series<-lynx
chatPasse<-chatterjee(series)
obj<-chatPasse$foci_path
df <- data.frame(lag = as.factor(obj$variable),correlation =obj$codec,step=obj$step)
df$lag <- factor(df$lag, levels = df$lag[order(df$step)])
df$step <- factor(df$step)
df$color <- ifelse(df$correlation < 0, "red", "steelblue")

codec_plot <- ggplot(df, aes(x = step, y = correlation)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_segment(aes(xend = step, yend = 0, color = color), linewidth = 1.2) +
  scale_color_identity() +  # use actual color values
  geom_text(aes(label = lag), vjust = ifelse(df$correlation >= 0, -0.5, 1.5), size = 4, color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal() +
  labs(
    title = "Selected lag in each iteration \n(CODEC value prop.to bar length)",
    x = "Iteration",
    y = ""
  )

df_ts <- data.frame(
  time = as.Date(time(series)),
  values = as.numeric(series)
)

# Time series plot
ts_plot <- ggplot(df_ts, aes(x = time, y = values)) +
  geom_area(fill = "#FFDBBB", alpha = 0.5) +  # filled area under the line
  geom_line(color = "orange", linewidth = 1.2) +  # smooth line
  geom_point(color = "black", size = 1) +  # dot at each observation
  theme(plot.title = element_text(hjust = 1)) +
  theme_minimal() +
  labs(
    title = "Original time series",
    x = "",
    y = "Lynx count"
  )


combined_plot <- ts_plot + codec_plot +
  plot_layout(ncol = 2, widths = c(1.3, 1)) +
  plot_annotation(
    title = "Analysis of Annual Count of Lynx",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
x11()
combined_plot


#####ANALYSIS FOR THE PASSENGER SERIES#####

series<-passengers
chatPasse<-chatterjee(series)
obj<-chatPasse$foci_path
df <- data.frame(lag = as.factor(obj$variable),correlation =obj$codec,step=obj$step)
df$lag <- factor(df$lag, levels = df$lag[order(df$step)])
df$step <- factor(df$step)
df$color <- ifelse(df$correlation < 0, "red", "steelblue")

codec_plot <- ggplot(df, aes(x = step, y = correlation)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_segment(aes(xend = step, yend = 0, color = color), linewidth = 1.2) +
  scale_color_identity() +  # use actual color values
  geom_text(aes(label = lag), vjust = ifelse(df$correlation >= 0, -0.5, 1.5), size = 4, color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal() +
  labs(
    title = "Selected lag in each iteration \n(CODEC value prop.to bar length)",
    x = "Iteration",
    y = ""
  )

df_ts <- data.frame(
  time = as.Date(time(series)),
  values = as.numeric(series)
)

# Time series plot
ts_plot <- ggplot(df_ts, aes(x = time, y = values)) +
  geom_area(fill = "#FFDBBB", alpha = 0.5) +  # filled area under the line
  geom_line(color = "orange", linewidth = 1.2) +  # smooth line
  geom_point(color = "black", size = 1) +  # dot at each observation
  theme(plot.title = element_text(hjust = 1)) +
  theme_minimal() +
  labs(
    title = "Original time series",
    x = "",
    y = "Passenger count"
  )


combined_plot <- ts_plot + codec_plot +
  plot_layout(ncol = 2, widths = c(1.3, 1)) +
  plot_annotation(
    title = "Analysis of Monthly Air Passengers",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
x11()
combined_plot
#####ANALYSIS FOR THE TEMPERATURE SERIES#####

series<-melbourne_ts
chatPasse<-chatterjee(series)
obj<-chatPasse$foci_path
df <- data.frame(lag = as.factor(obj$variable),correlation =obj$codec,step=obj$step)
df$lag <- factor(df$lag, levels = df$lag[order(df$step)])
df$step <- factor(df$step)
df$color <- ifelse(df$correlation < 0, "red", "steelblue")

codec_plot <- ggplot(df, aes(x = step, y = correlation)) +
  geom_hline(yintercept = 0, color = "gray50") +
  geom_segment(aes(xend = step, yend = 0, color = color), linewidth = 1.2) +
  scale_color_identity() +  # use actual color values
  geom_text(aes(label = lag), vjust = ifelse(df$correlation >= 0, -0.5, 1.5), size = 4, color = "black") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_minimal() +
  labs(
    title = "Selected lag in each iteration \n(CODEC value prop.to bar length)",
    x = "Iteration",
    y = ""
  )


df_ts <- data.frame(
  time = time(series), #only for Melbourne
  values = as.numeric(series)
)

# Time series plot
ts_plot <- ggplot(df_ts, aes(x = time, y = values)) +
  geom_area(fill = "#FFDBBB", alpha = 0.5) +  # filled area under the line
  geom_line(color = "orange", linewidth = 1.2) +  # smooth line
  geom_point(color = "black", size = 1) +  # dot at each observation
  theme(plot.title = element_text(hjust = 1)) +
  theme_minimal() +
  labs(
    title = "Original time series",
    x = "",
    y = "Min. Temperature (°C)"
  )


combined_plot <- ts_plot + codec_plot +
  plot_layout(ncol = 2, widths = c(1.3, 1)) +
  plot_annotation(
    title = "Analysis of Daily min. Temperatures in Melbourne",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
  )
x11()
combined_plot
