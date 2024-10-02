# wind speed estimation for Block maxima
# out should be the result of a bootstrap
WSEst_BM <- function(out, trueTail=NULL, bm=NULL,data=NULL,block_length=365.25, alpha = 0.05,plot_all_lines=F,col = 'F010',length.out=30){ # we consider we discard feb 29 from datasets so no leap year
  # block_length <- as.numeric(days(block_length),'days') 
  
  if (is.null(bm)){bm <- BM_select(data,block_length,height = col)}
  if (is.null(trueTail)){trueTail <- BM_fit(bm,block_length,plot=FALSE)}
  
  
  bm <- bm$max
  
  # print(bm$max)
  z <- qnorm(1 - alpha / 2)
  
  X <- 10^seq(log10(2),log10(10000),length.out = length.out)*365.25/block_length
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  return_values <- matrix(NA, nrow = nrow(out$df), ncol = length(X))
  # Compute the return values
  for (i in 1:nrow(out$df)){
    return_values[i,] <- qgev(p=exp(- 1 / X), loc = out$df$loc[i], scale = out$df$scale[i], tail = out$df$tail[i])
    if (any(return_values[i,]>3000)){
      return_values[i,]<-NA
    }
  }
  
  
  points <- data.frame(speed = sort(bm),return = (block_length/365.25)*1/(1-seq(1/length(bm),1-1/length(bm),length.out=length(bm))))            # observed Values
  
  
  
  
  
  if (plot_all_lines){
    
    data = data.frame(
      return = X*block_length/365.25,
      original = qgev(p=exp( - 1 / X), loc = trueTail['loc'], scale = trueTail['scale'], tail = trueTail['tail']), # fit for original dataset
      boot = t(return_values)
    )    
    
    df_plot = cbind(return = X*block_length/365.25,melt(data.frame(boot=t(return_values))))
    original <- data.frame(
      return =  X*block_length/365.25, 
      original = qgev(p=exp( - 1 / X), loc = trueTail['loc'], scale = trueTail['scale'], tail = trueTail['tail']) # fit for original dataset
    )
    print(df_plot)
    
    gg = ggplot(df_plot,aes(x=return)) +
      geom_line(aes(y=value,line=variable),col='gray') +
      geom_line(data = original,aes(x=return,y=original),col='red') +
      scale_x_log10()
    
  }else{
    # Create a data frame for plotting
    original = qgev(p=exp(- 1 / X), loc = trueTail['loc'], scale = trueTail['scale'], tail = trueTail['tail']) # fit for original dataset
    sd = apply(log(return_values),2,sd,na.rm=TRUE)

    data <- data.frame(
      return = X*block_length/365.25,
      original = original,
      lb = original*exp(-z*sd),
      ub = original*exp(z*sd),
      sd = sd
    )
    
    # Plot using ggplot2
    gg = ggplot(data, aes(x = return)) +
      geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval, 95% bootstrap'), alpha = 0.2) +
      # geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval,  95% log(bootstrap)'), alpha = 0.2) +
      geom_line(aes(x=return,y=original,linetype = 'Original estimation'),color='red') +
      geom_point(data=points,aes(x = return,y = speed)) +
      scale_x_log10(labels= scales::comma) +
      labs(x = 'Return Period (years)', y = 'Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years, ', block_length ,'d BM (GEV)')) +
      scale_linetype_manual(name="", values = c('Original estimation' = "solid")) +
      scale_fill_manual(name = "", values = c('Confidence Interval, 95% bootstrap' = 'blue','Confidence Interval,  95% log(bootstrap)'='red')) +
      theme_minimal() +
      theme(legend.position = "right")
    
    last_values <- tail(data, 1)  # Get the last row of df_measure
    # print(last_values)
    # Annotate the last values
    gg <- gg +
      geom_text(data = last_values, aes(label = sprintf("%.2f", original), 
                                        x = max(data$return), y = original ), 
                color = 'red', size = 3, hjust = 0) +
      geom_text(data = last_values, aes(label = sprintf("%.2f", lb), 
                                        x = max(data$return), y = lb), 
                color = 'red', size = 3, hjust = 0) +
      geom_text(data = last_values, aes(label = sprintf("%.2f", ub), 
                                        x = max(data$return), y = ub), 
                color = 'red', size = 3, hjust = 0)
    
    closest_to_50 <- data %>% 
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
    
  }
  
  return(list(gg=gg, pred=data))
}



WSEst_model_to_measure_BM <- function( data_model,data_measure, col_model='PF010', col_measure='F010', Nb_model=100, Nb_measure=100, tag=NULL,length.out=30, fixtailto0=F) {
  if (!fixtailto0){
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Bootstrap on model\n")
  flush.console()
  fit_model <- bootstrap(data_model, Nb = Nb_model, method = 'BM', column = col_model)
  }else{
    fit_model = list()
    fit_model$original['tail'] = 0
    fit_model$df$tail = rep(0,Nb_model)
    
  }
  # Get the block maxima from measurement dataset
  bm <- BM_select(data_measure, height = col_measure)
  
  # Function to constrain the tail index parameter in the estimation
  fpar <- function(p, xpar) {
    loc <- matrix(p[1], xpar[1], xpar[2])
    scale <- matrix(p[2], xpar[1], xpar[2])
    shape <- matrix(fit_model$original['tail'], xpar[1], xpar[2])
    list(loc = loc, scale = scale, shape = shape)
  }
  
  # Arguments for the BM_fit function
  kwargs <- list(start = c(mean(bm$max), sd(bm$max)), fpar = fpar, xpar = c(length(bm$max), 1))
  
  # Compute scale and loc for original and tail_index model estimated
  trueTail <- BM_fit(BM_select(data_measure, height = col_measure), plot = FALSE)
  modelTail <- BM_fit(BM_select(data_measure, height = col_measure), kwargs = kwargs, plot = FALSE)
  
  # Quantile for confidence interval (CI)
  z <- qnorm(1 - 0.05 / 2)
  # Return periods
  X <- 10^seq(log10(2), log10(10000), length.out = length.out)
  X <- sort(unique(c(c(50,100,1000,10000),X)))
  
  # Initialize matrix to store return values
  return_values <- matrix(NA, nrow = Nb_model * Nb_measure, ncol = length(X))
  
  for (i in 1:Nb_model) {
    cat('Bootstrap on Measurements:', i, '/', Nb_model, '\n')
    flush.console()
    
    # Update fpar function for each bootstrap iteration
    fpar <- function(p, xpar) {
      loc <- matrix(p[1], xpar[1], xpar[2])
      scale <- matrix(p[2], xpar[1], xpar[2])
      shape <- matrix(fit_model$df$tail[i], xpar[1], xpar[2])
      list(loc = loc, scale = scale, shape = shape)
    }
    
    # Arguments for the constrained bootstrap on measurements
    kwargs <- list(start = c(mean(bm$max), sd(bm$max)), fpar = fpar, xpar = c(length(bm$max), 1))
    
    # Perform bootstrap on measurement data
    fit_measure <- bootstrap(data_measure, Nb = Nb_measure,column = col_measure, method = 'BM', kwargs = kwargs)
    
    # Compute return values for each bootstrap sample
    for (j in 1:Nb_measure) {
      return_values[(i-1) * Nb_measure + j, ] <- qgev(p = exp( - 1 / X), loc = fit_measure$df$loc[j], scale = fit_measure$df$scale[j], tail = fit_measure$df$tail[j])
    }
  }
  
  # Compute standard deviation of the log of the return values
  sd <- apply(log(return_values), 2, sd,na.rm=TRUE)
  # Estimate return values using the model-estimated tail index
  model_est <- qgev(p = exp(- 1 / X), loc = modelTail['loc'], scale = modelTail['scale'], tail = modelTail['tail'])
  
  # Create a data frame for observed values
  points <- data.frame(speed = sort(bm$max), return = 1 / (1 - seq(1 / length(bm$max), 1 - 1 / length(bm$max), length.out = length(bm$max))))
  
  # Compute prediction intervals from the bootstrap measurements
  measure_est <- WSEst_BM(bootstrap(data_measure, Nb = Nb_measure, column = col_measure, method = 'BM'),col = col_measure, data = data_measure,length.out=length.out)$pred
  
  
  # Create a data frame for plotting
  df <- data.frame(
    return = X,
    original = measure_est$original,
    lb_o = measure_est$lb,
    ub_o = measure_est$ub,
    sd_o = measure_est$sd,
    model_est = model_est,
    lb = model_est * exp(-z * sd),
    ub = model_est * exp(z * sd),
    sd = sd,
    proba = 1/X
  )


  
  # Create the ggplot
  custom_color <- "blue"
  
  gg <- ggplot(df, aes(x = return)) +
    geom_line(aes(y = model_est, color = 'shape estimated with model')) +
    geom_line(aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
    geom_ribbon(aes(ymin = lb, ymax = ub, fill = 'Confidence Interval, 95%'), alpha = 0.2) +
    geom_line(aes(y = lb_o), lty = 'dashed', color = "red") +
    geom_line(aes(y = ub_o), lty = 'dashed', color = "red") +
    geom_point(data = points, aes(y = speed, shape = 'Measurements observed')) +
    scale_x_log10(labels = scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\n1y BM (GEV), on measurements')) +
    scale_color_manual(values = c('shape estimated with model' = custom_color, 'all parameters estimated\nfrom measurements' = 'red')) +
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
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = model_est), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "blue", size = 3, hjust = 0) +
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
  return(list(gg = gg, df = df))
}




# wind speed estimation for POT method with the first POT selection (deprecated, use second version)
# timestep in hours
# threshold : the percentage of data w want to use for the estimation
WSEst_GP <- function(data,timestep,threshold,col = "F010"){
  
  pot <- PoTselect(na.approx(data[[col]]),0.1,12/timestep)
  print(length(pot$pot))
  l0 <- round(threshold*length(pot$pot))-1
  pot_fit <- FitGP_MLE2(pot$pot,1, N= 0, r11= 1, fixedpar= NULL, l0= l0, sigma= Inf)
  
  # plot_extremes(pot,data)
  # plot(ecdf(pot$pot))
  # lines(seq(min(pot$pot),max(pot$pot),length.out=100),cgp((seq(min(pot$pot),max(pot$pot),length.out=100)-pot_fit$location)/pot_fit$scale,gamma = pot_fit$tailindex),col='red')
  
  
  timestep_selected_peaks = mean(diff(pot$ind)*timestep/(24*365.2425)) # mean of the timediff between each peak
  
  
  points <- data.frame(speed = sort(pot$pot))
  points$p <- seq(0,1-1/length(points$speed),length.out=length(points$speed))     # observed Values
  points$used <- 1-points$p <= l0/length(pot$pot)
  points$return <- timestep_selected_peaks / (1- points$p)

  
  # print(ggplot(points,aes(x=speed,y=p,color=used))+geom_line()+geom_point()+scale_color_manual(name="",values = c('TRUE'='black','FALSE'='lightgrey')))
  
  
  tail <- list(mean=pot_fit$tailindex,sd=pot_fit$tailindexStd)
  scale <- list(mean=pot_fit$scale,sd=pot_fit$logdispStd)
  loc <- list(mean=pot_fit$location,sd=pot_fit$locationStd)
  
  X <- 10^seq(log10(timestep_selected_peaks*length(pot$pot)/l0),log10(10000),length.out = 30)
  

  
  distrib = qgp_distrib(1-timestep_selected_peaks/X,tail = tail,scale = scale,loc = loc,t=length(pot$pot)/l0)
  
  
  df <- data.frame(
    return = X,
    original = distrib$mean,
    lb = distrib$lb,
    ub = distrib$ub
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





WSEst_model_to_measure_GP <- function(data_model, data_measure, col_model='PF010',col_measure='F010',timestep_model=1,timestep_measure=1, th_model=NULL, th_measure=NULL, tag=NULL) {
  
  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Estimation on model\n")
  flush.console()
  pot_model <- PoTselect(data_model[[col_model]],0.1,12/timestep_model)
  l0_model <- round(th_model*length(pot_model$pot))-1
  fit_model <- FitGP_MLE2(pot_model$pot, 1, N= nrow(data_model), r11= 1, l0= l0_model)
  
  
  
  
  # Function to constrain the tail index parameter in the estimation
  fixedpar = list(gamma0 = fit_model$tailindex,gamma0Std = fit_model$tailindexStd)
  
  # Compute scale and loc for original and tail_index model estimated
  
  
  pot_measure <- PoTselect(data_measure[[col_measure]],0.1,12/timestep_measure)
  print(sum(is.na(pot_measure$pot)))
  l0_measure <- round(th_measure*length(pot_measure$pot))-1
  fit_measure <- FitGP_MLE2(pot_measure$pot, 1, N= nrow(data_measure), r11= 1, l0= l0_measure)
  fit_measure_model <- FitGP_MLE2(pot_measure$pot, 1, N= nrow(data_measure),fixedpar = fixedpar, r11= 1, l0= l0_measure)
  
  # plot ECDFs
  xrange = c(min(c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model])),max(c(pot_measure$pot,pot_model$pot)))
  xspan = seq(xrange[1],xrange[2],length.out=100)
  
  plot(ecdf(-sort(-pot_model$pot)),pch=1,col='black',
       xlim = xrange,
       ylim = c(1-max(th_model,th_measure),1),
       xlab = 'wind speed (m/s)',
       ylab = 'probability', 
       main = 'CDF of POT on measurements and model')
  lines(ecdf(-sort(-pot_measure$pot)),cex = 0.1,col='red')
  lines(xspan,cgp(xspan, loc = fit_model$location,scale = fit_model$scale,tail = fit_model$tailindex, t = length(pot_model$pot)/l0_model),lty='dashed',col='green')
  lines(xspan,cgp(xspan,loc = fit_measure$location, scale = fit_measure$scale, tail = fit_measure$tailindex, t = length(pot_measure$pot)/l0_measure),lty='dashed',col='blue')
  lines(xspan,cgp(xspan, loc = fit_measure_model$location, scale = fit_measure_model$scale,tail = fit_measure_model$tailindex, t = length(pot_measure$pot)/l0_measure),lty=3,col='orange')
  abline(v=c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model]),col=c('red','black'))
  legend('bottomright',cex = 0.6,legend=c("ecdf PoT model","ecdf PoT measurements", "cgp fit model", "cgp fit measurements","cgp fit model\n+ measurements"),
         col = c('black','red','green','blue','orange'),
         pch = c(1,1,NA,NA,NA),
         lty = c(1,1,2,2,3))
  
  print("")

  

  # mean timestep for POT select on measurements
  timestep_POT_measure <- mean(diff(pot_measure$ind)*timestep_measure/(24*365.2425)) # mean of the timediff between each peak
  
  # Return periods
  X <- 10^seq(log10(timestep_POT_measure*length(pot_measure$pot)/l0_measure), log10(10000), length.out = 30)
  
  
  print(paste('Model estimation : tail =',fit_model$tailindex,'sd',fit_model$tailindexStd,'scale =',fit_model$scale,' loc =',fit_model$location))
  # measurement original estimation
  loc = list(mean = fit_measure$location, sd = fit_measure$locationStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$logdispStd)
  tail = list(mean = fit_measure$tailindex, sd = fit_measure$tailindexStd)
  print(paste('Original estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  measure_est = qgp_distrib(1-timestep_POT_measure/X,tail = tail,scale = scale,loc = loc,t=length(pot_measure$pot)/l0_measure)
  
  df_measure <- data.frame(
    return = X,
    original = measure_est$mean,
    lb = measure_est$lb,
    ub = measure_est$ub
  )
  
  
  # Estimate return values using the model-estimated tail index
  loc = list(mean = fit_measure_model$location, sd = fit_measure_model$locationStd)
  scale = list(mean = fit_measure_model$scale, sd = fit_measure_model$logdispStd)
  tail = list(mean = fit_measure_model$tailindex, sd = fit_measure_model$tailindexStd)
  print(paste('Measurements + Model estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  model_est <- qgp_distrib(p = 1 - timestep_POT_measure / X, loc = loc, scale = scale, tail = tail,t=length(pot_measure$pot)/l0_measure)
  
  
  
  # Create a data frame for plotting
  df <- data.frame(
    return = X,
    original = df_measure$original,
    model_est = model_est$mean,
    lb = model_est$lb,
    ub = model_est$ub
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
    geom_point(data = points, aes(y = speed, color=used,shape = 'Measurements observed')) +
    geom_line(data = df_measure, aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
    geom_line(data = df_measure, aes(y = lb), lty = 'dashed', color = "red") +
    geom_line(data = df_measure, aes(y = ub), lty = 'dashed', color = "red") +
    scale_x_log10(labels = scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\nPoT, GP, on measurements')) +
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
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = model_est), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "blue", size = 3, hjust = 0) +
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
  return(list(gg = gg, df = df))
}


# wind speed estimation for POT method
# timestep in hours
# threshold : the percentage of data w want to use for the estimation
WSEst_GW <- function(data,timestep,threshold,col="F010"){
  
  pot <- PoTselect(na.approx(data[[col]]),0.1,12/timestep)
  l0 <- round(threshold*length(pot$pot))-1
  pot_fit <- FitGW_iHill(pot$pot,p=1,l0 = l0)
  
  
  timestep_selected_peaks = mean(diff(pot$ind)*timestep/(24*365.2425)) # mean of the timediff between each peak
  
  
  points <- data.frame(speed = sort(pot$pot))
  points$p <- seq(0,1-1/length(points$speed),length.out=length(points$speed))     # observed Values
  points$used <- 1-points$p <= l0/length(pot$pot)
  points$return <- timestep_selected_peaks / (1- points$p)
  
  
  print(ggplot(points,aes(x=speed,y=p,color=used))+geom_line()+geom_point()+scale_color_manual(name="",values = c('TRUE'='black','FALSE'='lightgrey')))
  
  
  tail <- list(mean=pot_fit$tailindex,sd=pot_fit$tailindexStd)
  scale <- list(mean=pot_fit$scale,sd=pot_fit$logdispStd)
  loc <- list(mean=pot_fit$location,sd=pot_fit$locationStd)
  
  X <- 10^seq(log10(timestep_selected_peaks*length(pot$pot)/l0),log10(10000),length.out = 30)
  
  
  
  distrib = qgw_distrib(1-timestep_selected_peaks/X,tail = tail,scale = scale,loc = loc,y=pot_fit$y)
  
  
  df <- data.frame(
    return = X,
    original = distrib$mean,
    lb = distrib$lb,
    ub = distrib$ub
  )
  
  
  gg = ggplot(df, aes(x = return)) +
    geom_line(aes(x=return,y=original),col='red') +
    geom_ribbon(aes(ymin=lb,ymax=ub),col='red',alpha=0.2)+
    geom_point(data=points,aes(x = return,y = speed,color=used)) +
    geom_vline(xintercept = timestep_selected_peaks*length(pot$pot)/l0,lty='dashed',alpha=0.3)+
    scale_x_log10(labels= scales::comma) +
    labs(x = 'Return Period (years)',
         y = 'Wind Speed (m/s)',
         title = paste0('Estimated wind speed for return period until 10000 years,\nPoT: ',threshold*100,'% of the peaks, GW')) +
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





WSEst_model_to_measure_GW <- function(data_model, data_measure, col_model='PF010',col_measure='F010',timestep_model=1,timestep_measure=1, th_model=NULL, th_measure=NULL, tag=NULL) {

  # Bootstrap on model to get uncertainties of tail index estimation
  cat("Estimation on model\n")
  flush.console()
  pot_model <- PoTselect(data_model[[col_model]],0.1,12/timestep_model)
  l0_model <- round(th_model*length(pot_model$pot))-1
  fit_model <- FitGW_iHilli(pot_model$pot,p= 1, l0= l0_model)
  
  
  
  
  # Function to constrain the tail index parameter in the estimation
  fixedpar = list(theta0 = fit_model$tailindex,theta0Std = fit_model$tailindexStd)
  
  # Compute scale and loc for original and tail_index model estimated
  
  
  pot_measure <- PoTselect(na.approx(data_measure[[col_measure]]),0.1,12/timestep_measure)
  l0_measure <- round(th_measure*length(pot_measure$pot))-1
  fit_measure <- FitGW_iHilli(pot_measure$pot,p= 1, l0= l0_measure)
  fit_measure_model <- FitGW_iHilli(pot_measure$pot,p= 1,fixedpar = fixedpar, l0= l0_measure)
  
  # plot ECDFs
  xrange = c(min(c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model])),max(c(pot_measure$pot,pot_model$pot)))
  xspan = seq(xrange[1],xrange[2],length.out=100)
  
  if (!is.null(tag)){
  plot(ecdf(-sort(-pot_model$pot)),pch=1,col='black',
       xlim = xrange,
       ylim = c(1-max(th_model,th_measure),1),
       xlab = 'wind speed (m/s)',
       ylab = 'probability', 
       main = 'CDF of POT on measurements and model')
  lines(ecdf(-sort(-pot_measure$pot)),cex = 0.1,col='red')
  lines(xspan,cgw(xspan, loc = fit_model$location,scale = fit_model$scale,tail = fit_model$tailindex, y= fit_model$y),lty='dashed',col='green')
  lines(xspan,cgw(xspan,loc = fit_measure$location, scale = fit_measure$scale, tail = fit_measure$tailindex,  y= fit_measure$y),lty='dashed',col='blue')
  lines(xspan,cgw(xspan, loc = fit_measure_model$location, scale = fit_measure_model$scale,tail = fit_measure_model$tailindex,  y= fit_measure_model$y),lty=3,col='orange')
  abline(v=c(-sort(-pot_measure$pot)[l0_measure],-sort(-pot_model$pot)[l0_model]),col=c('red','black'))
  legend('bottomright',cex = 0.6,legend=c("ecdf PoT model","ecdf PoT measurements", "cgw fit model", "cgw fit measurements","cgw fit model\n+ measurements"),
         col = c('black','red','green','blue','orange'),
         pch = c(1,1,NA,NA,NA),
         lty = c(1,1,2,2,3))
  }
  print("")
  
  
  
  
  # mean timestep for POT select on measurements
  timestep_POT_measure <- mean(diff(pot_measure$ind)*timestep_measure/(24*365.2425)) # mean of the timediff between each peak
  
  # Return periods
  X <- 10^seq(log10(timestep_POT_measure*length(pot_measure$pot)/l0_measure), log10(10000), length.out = 30)
  
  
  print(paste('Model estimation : tail =',fit_model$tailindex,'sd',fit_model$tailindexStd,'scale =',fit_model$scale,' loc =',fit_model$location))
  # measurement original estimation
  loc = list(mean = fit_measure$location, sd = fit_measure$locationStd)
  scale = list(mean = fit_measure$scale, sd = fit_measure$logdispStd)
  tail = list(mean = fit_measure$tailindex, sd = fit_measure$tailindexStd)
  print(paste('Original estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  measure_est = qgw_distrib(1-timestep_POT_measure/X,tail = tail,scale = scale,loc = loc, y= fit_measure$y)
  
  
  
  # Estimate return values using the model-estimated tail index
  loc = list(mean = fit_measure_model$location, sd = fit_measure_model$locationStd)
  scale = list(mean = fit_measure_model$scale, sd = fit_measure_model$logdispStd)
  tail = list(mean = fit_measure_model$tailindex, sd = fit_measure_model$tailindexStd)
  print(paste('Measurements + Model estimation : tail =',tail$mean,'sd',tail$sd,'scale =',scale$mean,' loc =',loc$mean))
  
  model_est <- qgw_distrib(p = 1 - timestep_POT_measure / X, loc = loc, scale = scale, tail = tail, y= fit_measure_model$y)
  
  
  
  # Create a data frame for plotting
  df <- data.frame(
    return = X,
    original = measure_est$mean,
    lb_o = measure_est$lb,
    ub_o = measure_est$ub,
    model_est = model_est$mean,
    lb = model_est$lb,
    ub = model_est$ub
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
    geom_point(data = points, aes(y = speed, color=used,shape = 'Measurements observed')) +
    geom_line(aes(y = original, color = 'all parameters estimated\nfrom measurements')) +
    geom_line(aes(y = lb_o), lty = 'dashed', color = "red") +
    geom_line(aes(y = ub_o), lty = 'dashed', color = "red") +
    scale_x_log10(labels = scales::comma) +
    ylim(range(c(df$original, df$model_est, df$lb, df$ub, points$speed))) +
    labs(x = 'Return Period (years)', y = 'Return Wind Speed (m/s)', title = paste0('Estimated wind speed for return period until 10000 years,\nPoT + GW, on measurements')) +
    scale_color_manual(values = c('shape estimated with model' = custom_color, 'all parameters estimated\nfrom measurements' = 'red','TRUE'='black','FALSE'='lightgrey')) +
    scale_fill_manual(values = c('Confidence Interval, 95%' = custom_color)) +
    theme(legend.position = "right")
  
  last_values <- tail(df, 1)  # Get the last row of df_measure
  # print(last_values)
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
    geom_point(data = closest_to_50, aes(x = return, y = original), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = model_est), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = lb), color = "blue", size = 1) +
    geom_point(data = closest_to_50, aes(x = return, y = ub), color = "blue", size = 1) +
    geom_text(data = closest_to_50, aes(label = sprintf("%.2f", original), 
                                        x = return, y = original + 1), 
              color = "blue", size = 3, hjust = 0) +
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
  return(list(gg = gg, df = df,points = points))
}
