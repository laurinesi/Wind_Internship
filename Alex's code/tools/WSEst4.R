# take dependence
#all analytic

WSEst_model_to_measure_BM_4 <- function(data_model,data_measure, col_model='PF010', col_measure='F010', tag=NULL,length.out=30,cov_type="num",fixtailto0=F) {
  
  
  if (fixtailto0){
    tail_model <- list(tail = 0,tailStd = 0)
  }else{
    bm_model <- BM_select(data_model,height = col_model)

    tail_model <- BM_fit_cov(bm_model,plot = F,cov_type=cov_type)

    tail_model <- list(tail = tail_model$trueTail['tail'],tailStd = tail_model$cov['tail','tail'])
    }
  
  # Get the block maxima from measurement dataset
  bm <- BM_select(data_measure, height = col_measure)
  
  # Function to constrain the tail index parameter in the estimation
  fpar <- function(p, xpar) {
    loc <- matrix(p[1], xpar[1], xpar[2])
    scale <- matrix(p[2], xpar[1], xpar[2])
    shape <- matrix(tail_model$tail, xpar[1], xpar[2])
    list(loc = loc, scale = scale, shape = shape)
  }
  
  
  kwargs <- list(start = c(mean(bm$max), sd(bm$max)), fpar = fpar, xpar = c(length(bm$max), 1), tailStd = tail_model$tailStd)
 
  
  # Compute scale and loc for original and tail_index model estimated
  
  
  # Quantile for confidence interval (CI)
  z <- qnorm(1 - 0.05 / 2)
  
  fit_measure <- BM_fit_cov(bm,plot = F,cov_type=cov_type)

  fit_measure_model <- BM_fit_cov(bm,plot = F,kwargs = kwargs,cov_type=cov_type)
  print("test")
  
  # Return periods
  X <- 10^seq(log10(5), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  print(paste('Model estimation : tail =',tail_model$tail,'sd',tail_model$tailStd))
  # measurement original estimation
  parameters = list(mean = fit_measure$trueTail,cov = fit_measure$cov)
  measure_est = qgev_distrib_2(exp(-1/X),parameters)
  print('Original estimation : ')
  print(parameters)
  
  # measurement original estimation
  parameters = list(mean = fit_measure_model$trueTail,cov = fit_measure_model$cov)
  model_est = qgev_distrib_2(exp(-1/X),parameters)
  print('Measurements + model tail estimation : ')
  print(parameters)
  
  print(model_est)
  
  # Create a data frame for plotting
  df <- data.frame(
    return = X,
    original = measure_est$mean,
    lb_o = measure_est$lb,
    ub_o = measure_est$ub,
    sd_o = measure_est$sd,
    model_est = model_est$mean,
    lb = model_est$lb,
    ub = model_est$ub,
    sd = model_est$sd,
    proba = 1/X
  )
  
  
  
  # Create a data frame for observed values
  points <- data.frame(speed = sort(bm$max), return = 1 / (1 - seq(1 / length(bm$max), 1 - 1 / length(bm$max), length.out = length(bm$max))))
  
  
  
  {
  # Create the ggplot
  custom_color <- "blue"
  
  gg <- ggplot(df, aes(x = return)) +
    geom_line(aes(y = model_est, color = 'shape estimated with model')) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval, 95%'), alpha = 0.2) +
    geom_line(aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
    geom_line(aes(y = lb_o), lty = 'dashed', color = "red") +
    geom_line(aes(y = ub_o), lty = 'dashed', color = "red") +
    geom_point(data = points, aes(y = speed,shape = 'Measurements observed')) +
    scale_x_log10(labels = scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\n1y BM (GEV), on measurements')) +
    scale_color_manual(values = c('shape estimated with model' = custom_color, 'all parameters estimated\nfrom measurements' = 'red','TRUE'='black','FALSE'='lightgrey')) +
    scale_fill_manual(values = c('Confidence Interval, 95%' = custom_color)) +
    theme(legend.position = "right")
  
  last_values <- tail(df, 1)  # Get the last row of df_measure
  print(last_values)
  # Annotate the last values
  gg <- gg +
    geom_text(data = last_values, aes(label = sprintf("%.2f", original), 
                                      x = max(df$return), y = original ), 
              color = 'red', size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", model_est), 
                                      x = max(df$return), y = model_est), 
              color = custom_color, size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", lb), 
                                      x = max(df$return), y = lb), 
              color = custom_color, size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", ub), 
                                      x = max(df$return), y = ub), 
              color = custom_color, size = 3, hjust = 0)
  
  closest_to_50 <- df %>% 
    slice(which.min(abs(return - 50)))
  
  gg <- gg +  
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "red", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = model_est), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "red", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", model_est), 
                                        x = return, y = model_est + 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", lb), 
                                        x = return, y = lb - 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", ub), 
                                        x = return, y = ub + 1), 
              color = "blue", size = 3, hjust = 0)
  }#plot
  
  # Return a list containing the ggplot object, data frame, and observed points
  parameter_distributions = list(
    model = tail_model,
    measure = fit_measure,
    measure_model = fit_measure_model
  )
  return(list(gg = gg, df = df,timestep = 1,points=points, parameter_distributions = parameter_distributions))
}



WSEst_GP_4 <- function(data,timestep,threshold,col="F010",length.out = 30,fixtailto0=F,winter=F){
  
  ws <- na.approx(data[[col]])
  if (winter){
    ws <- na.approx(data[[col]][month(data$DateTime) %in% c(10:12,1:3)])
  }
  pot <- PoTselect_2(ws,0.3,12/timestep)
  print(length(pot$pot))
  l0 <- round(threshold*length(pot$pot))-1
  
  if (fixtailto0){fixedpar <- list(gamma0 = 0,gamma0Std = 0)}
  else{fixedpar <- NULL}
  

  # plot_extremes(pot,data)
  # plot(ecdf(pot$pot))
  # lines(seq(min(pot$pot),max(pot$pot),length.out=100),cgp((seq(min(pot$pot),max(pot$pot),length.out=100)-pot_fit$location)/pot_fit$scale,gamma = pot_fit$tailindex),col='red')
  
  
  timestep_selected_peaks = mean(diff(pot$ind)*timestep/(24*365.2425)) # mean of the timediff between each peak
  
  
  points <- data.frame(speed = sort(pot$pot))
  points$p <- seq(0,1-1/length(points$speed),length.out=length(points$speed))     # observed Values
  points$used <- 1-points$p <= l0/length(pot$pot)
  points$return <- timestep_selected_peaks / (1- points$p)
  
  
  # print(ggplot(points,aes(x=speed,y=p,color=used))+geom_line()+geom_point()+scale_color_manual(name="",values = c('TRUE'='black','FALSE'='lightgrey')))
  
  
  
  X <- 10^seq(log10(timestep_selected_peaks*length(pot$pot)/l0),log10(10000),length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  pot_fit <- FitGP_MLE2(pot$pot,timestep_selected_peaks/X, N= 0, r11= 1, fixedpar= fixedpar, l0= l0, sigma= Inf)
  
  distrib = data.frame(mean = t(pot_fit$quantile),lb = t(pot_fit$quantile - qnorm(0.975)*pot_fit$quantileStd), ub = t(pot_fit$quantile + qnorm(0.975)*pot_fit$quantileStd), sd = t(pot_fit$quantileStd))
  
  
  df <- data.frame(
    return = X,
    original = distrib$mean,
    lb = distrib$lb,
    ub = distrib$ub,
    sd = distrib$sd
  )
  
  
  gg = ggplot(df, aes(x = return)) +
    geom_line(aes(x=return,y=original),col='red') +
    geom_ribbon(aes(ymin=lb,ymax=ub),col='red',alpha=0.2)+
    geom_point(data=points,aes(x = return,y = speed,color=used)) +
    geom_vline(xintercept = timestep_selected_peaks*length(pot$pot)/l0,lty='dashed',alpha=0.3)+
    scale_x_log10(labels= scales::comma) +
    labs(x = 'Return Period (years)',
         y = 'Wind Speed (m/s)',
         title = paste0('Estimated wind speed for return period until 10000 years,\nPoT: ',threshold*100,'% of the peaks, GP')) +
    theme(legend.position = "right")+
    scale_color_manual(name="",values = c('TRUE'='black','FALSE'='lightgrey'))
  
  last_values <- tail(df, 1)  # Get the last row of df_measure
  print(last_values)
  # Annotate the last values
  gg <- gg +
    geom_text(data = last_values, aes(label = sprintf("%.2f", original), 
                                      x = max(df$return), y = original ), 
              color = 'red', size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", lb), 
                                      x = max(df$return), y = lb), 
              color = 'red', size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", ub), 
                                      x = max(df$return), y = ub), 
              color = 'red', size = 3, hjust = 0)
  
  closest_to_50 <- df %>% 
    slice(which.min(abs(return - 50)))
  
  gg <- gg +  
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", lb), 
                                        x = return, y = lb - 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", ub), 
                                        x = return, y = ub + 1), 
              color = "blue", size = 3, hjust = 0)
  return(list(gg=gg,df=df,points=points))
}


WSEst_model_to_measure_GP_4 <- function(data_model, data_measure, col_model='PF010',col_measure='F010',timestep_model=1,timestep_measure=1, th_model=NULL, th_measure=NULL, length.out = 100,peak_frac=0.3,winter=F) {
  
  ws_model <- na.approx(data_model[[col_model]])
  ws_measure <- na.approx(data_measure[[col_measure]])
  if (winter){
    ws_model <- na.approx(data_model[[col_model]][month(data_model$DateTime) %in% c(10:12,1:3)])
    ws_measure <- na.approx(data_measure[[col_measure]][month(data_measure$DateTime) %in% c(10:12,1:3)])
  }
  
  
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Estimation on model\n")
  flush.console()
  pot_model <- PoTselect_2(ws_model,peak_frac,12/timestep_model)
  l0_model <- round(th_model*length(pot_model$pot))-1
  fit_model <- FitGP_MLE2(pot_model$pot, 1, N= nrow(data_model), r11= 1, l0= l0_model)
  
  
  
  
  # Function to constrain the tail index parameter in the estimation
  fixedpar = list(gamma0 = fit_model$tailindex,gamma0Std = fit_model$tailindexStd)
  
  
  # POT and l0
  pot_measure <- PoTselect_2(ws_measure,peak_frac,12/timestep_measure)
  print(sum(is.na(pot_measure$pot)))
  l0_measure <- round(th_measure*length(pot_measure$pot))-1
  
  # mean timestep for POT select on measurements
  timestep_POT_measure <- (length(data_measure[[col_measure]])/length(pot_measure$pot))*timestep_measure/(24*365.2425)
  
  # Return periods
  X <- 10^seq(log10(timestep_POT_measure*length(pot_measure$pot)/l0_measure), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  
  # Compute scale and loc for original and tail_index model estimated
  
  fit_measure <- FitGP_MLE2(pot_measure$pot, p =timestep_POT_measure/X, r11= 1, l0= l0_measure)

  fit_measure_model <- FitGP_MLE2(pot_measure$pot, 1, N= nrow(data_measure),fixedpar = fixedpar, r11= 1, l0= l0_measure)
  

  
  

  

  print(paste('Model estimation : tail =',fit_model$tailindex,'sd',fit_model$tailindexStd,'scale =',fit_model$scale,' loc =',fit_model$location))
  
  
  # measurement original estimation (there is dependence)
  loc = list(mean = fit_measure$location, sd = fit_measure$locationStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$logdispStd)
  tail = list(mean = fit_measure$tailindex, sd = fit_measure$tailindexStd)
  
  print(paste('Original estimation (but there is dependence between) : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  measure_est = data.frame(mean = t(fit_measure$quantile),lb = t(fit_measure$quantile - qnorm(0.975)*fit_measure$quantileStd), ub = t(fit_measure$quantile + qnorm(0.975)*fit_measure$quantileStd), sd = t(fit_measure$quantileStd))
  
  # Estimate return values using the model-estimated tail index (all is independent)
  loc = list(mean = fit_measure_model$location, sd = fit_measure_model$locationStd)
  scale = list(mean = fit_measure_model$scale, sd = fit_measure_model$logdispStd)
  tail = list(mean = fit_measure_model$tailindex, sd = fit_measure_model$tailindexStd)
  print(paste('Measurements + Model estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  model_est <- qgp_distrib(p = 1- timestep_POT_measure / X, loc = loc, scale = scale, tail = tail,t=length(pot_measure$pot)/l0_measure)
  
  
  
  # Create a data frame for plotting
  df <- data.frame(
    return = X,
    original = measure_est$mean,
    lb_o = measure_est$lb,
    ub_o = measure_est$ub,
    sd_o = measure_est$sd,
    model_est = model_est$mean,
    lb = model_est$lb,
    ub = model_est$ub,
    sd = model_est$sd,
    proba = timestep_POT_measure/X
  )
  
  
  
  # Create a data frame for observed values
  
  points <- data.frame(speed = sort(pot_measure$pot))
  points$p <- seq(0,1-1/length(points$speed),length.out=length(points$speed))     # observed Values
  points$used <- 1-points$p <= l0_measure/length(pot_measure$pot)
  points$return <-  timestep_POT_measure / (1- points$p)
  
  
  
  
  # Create the ggplot
  custom_color <- "blue"
  
  gg <- ggplot(df, aes(x = return)) +
    geom_line(aes(y = model_est, color = 'shape estimated with model')) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval, 95%'), alpha = 0.2) +
    geom_line(aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
    geom_line(aes(y = lb_o), lty = 'dashed', color = "red") +
    geom_line(aes(y = ub_o), lty = 'dashed', color = "red") +
    geom_point(data = points, aes(y = speed, color=used,shape = 'Measurements observed')) +
    scale_x_log10(labels = scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\nPoT (GP), p=',peak_frac,', th=',th_measure,', on measurements')) +
    scale_color_manual(values = c('shape estimated with model' = custom_color, 'all parameters estimated\nfrom measurements' = 'red','TRUE'='black','FALSE'='lightgrey')) +
    scale_fill_manual(values = c('Confidence Interval, 95%' = custom_color)) +
    theme(legend.position = "right")
  
  last_values <- tail(df, 1)  # Get the last row of df_measure
  print(last_values)
  # Annotate the last values
  gg <- gg +
    geom_text(data = last_values, aes(label = sprintf("%.2f", original), 
                                      x = max(df$return), y = original ), 
              color = 'red', size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", model_est), 
                                      x = max(df$return), y = model_est), 
              color = custom_color, size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", lb), 
                                      x = max(df$return), y = lb), 
              color = custom_color, size = 3, hjust = 0) +
    geom_text(data = last_values, aes(label = sprintf("%.2f", ub), 
                                      x = max(df$return), y = ub), 
              color = custom_color, size = 3, hjust = 0)
  
  closest_to_50 <- df %>% 
    slice(which.min(abs(return - 50)))
  
  gg <- gg +  
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "red", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = model_est), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "red", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", model_est), 
                                        x = return, y = model_est + 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", lb), 
                                        x = return, y = lb - 1), 
              color = "blue", size = 3, hjust = 0) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", ub), 
                                        x = return, y = ub + 1), 
              color = "blue", size = 3, hjust = 0)
  # Return a list containing the ggplot object, data frame, and observed points
  parameter_distributions = list(
    model = list(tail = list(mean = fit_model$tailindex,sd= fit_model$tailindexStd),
                 scale = list(mean = fit_model$scale,sd= fit_model$logdispStd),
                 loc = list(mean = fit_model$location,sd= fit_model$locationStd)),
    measure = fit_measure,
    measure_model = list(tail=tail,scale=scale,loc=loc)
  )
  return(list(gg = gg, df = df,timestep = timestep_POT_measure,points=points, parameter_distributions = parameter_distributions))
}
