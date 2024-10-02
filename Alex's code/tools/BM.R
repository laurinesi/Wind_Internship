#With this function, we can fit the gev distribution and plot it over the ECDF and hist of the block maxima
BM_fit <- function(bm, block_length=365.25,plot=T, kwargs=NULL) {
  temp <- matrix(bm$max)
  
  if (is.null(kwargs)){
    start <- c(mean(temp),sd(temp),0) 
    xpar <- c(length(temp),1)
    fpar <- function(p, xpar) {
      loc <- matrix(p[1], xpar[1], xpar[2])
      scale <- matrix(p[2], xpar[1], xpar[2])
      shape <- matrix(p[3], xpar[1], xpar[2])
      list(loc = loc, scale = scale, shape = shape)
    }
  }else{
    start <- kwargs$start
    xpar <- kwargs$xpar
    fpar <- kwargs$fpar
  }
  # fit the GEV distribtuion and estimate tail index
  out <- fgev.flex(temp,start,fpar,xpar)
  if (length(out$estimate)!=3){
    out$estimate <- unlist(fpar(out$estimate,c(1,1)))
  }
  trueTail = c(out$estimate[3],out$estimate[1],out$estimate[2])
  names(trueTail) = c('tail','loc','scale')
  # print(paste('Tail index =',trueTail[1]))
  
  # plot the result
  xspan <- seq(min(temp), max(temp), length.out = 100)
  
  if(plot){
    ## density distribution
    hist(temp,breaks = 10,freq = F, main=paste0('extreme wind speed distribution, 10m,\n',block_length,'d BM'), xlab='wind speed (m/s)')
    lines(xspan,gev(xspan,loc = trueTail['loc'], scale = trueTail['scale'],tail = trueTail['tail']),col='red')
    
    ## cumulative distribution
    plot(ecdf(temp), main=paste0('extreme wind speed cumulative distribution, 10m,\n',block_length,'d BM'), xlab='wind speed (m/s)')
    lines(xspan,cgev(xspan,loc = trueTail['loc'], scale = trueTail['scale'],tail = trueTail['tail']),col='red')
  }
  return(trueTail)
}

BM_fit_cov <- function(bm, block_length=365.25,plot=T, kwargs=NULL,cov_type='num') {
  temp <- matrix(bm$max)

  if (is.null(kwargs)){
    start <- c(mean(temp),sd(temp),0) 
    xpar <- c(length(temp),1)
    fpar <- function(p, xpar) {
      loc <- matrix(p[1], xpar[1], xpar[2])
      scale <- matrix(p[2], xpar[1], xpar[2])
      shape <- matrix(p[3], xpar[1], xpar[2])
      list(loc = loc, scale = scale, shape = shape)
    }
  }else{
    start <- kwargs$start
    xpar <- kwargs$xpar
    fpar <- kwargs$fpar
  }
  # fit the GEV distribtuion and estimate tail index
  out <- fgev.flex(temp,start,fpar,xpar)
  
  if (length(out$estimate)!=3){
    out$estimate <- unlist(fpar(out$estimate,c(1,1)))
  }
  trueTail = c(out$estimate[1],out$estimate[2],out$estimate[3])
  names(trueTail) = c('loc','scale','tail')
  
  
  if (cov_type == 'analytic'){
      out$cov = gevinfom(trueTail['tail'],trueTail['scale'],length(temp))$cov
    }
  
  if (!is.null(kwargs)){
    cov = matrix(0,3,3)
    colnames(cov) <- c('loc','scale','tail')
    rownames(cov) <- c('loc','scale','tail')
    cov[c('loc','scale'),c('loc','scale')] = out$cov[c('loc','scale'),c('loc','scale')]
    cov['tail','tail'] = kwargs$tailStd
  }else{
    cov = out$cov
    colnames(cov) <- c('loc','scale','tail')
    rownames(cov) <- c('loc','scale','tail')
  }
  
  
  print(cov)
  
  if(plot){
    # plot the result
    xspan <- seq(min(temp), max(temp), length.out = 100)
    ## density distribution
    hist(temp,breaks = 10,freq = F, main=paste0('extreme wind speed distribution, 10m,\n',block_length,'d BM'), xlab='wind speed (m/s)')
    lines(xspan,gev(xspan,loc = trueTail['loc'], scale = trueTail['scale'],tail = trueTail['tail']),col='red')
    
    ## cumulative distribution
    plot(ecdf(temp), main=paste0('extreme wind speed cumulative distribution, 10m,\n',block_length,'d BM'), xlab='wind speed (m/s)')
    lines(xspan,cgev(xspan,loc = trueTail['loc'], scale = trueTail['scale'],tail = trueTail['tail']),col='red')
  }
  
  return(list(trueTail=trueTail,cov=cov))
}

#This function is to selct the block maxima
BM_select <- function(data, block_length=365.25, height='F010') {
  if (is.Date(data$Year)){data$Year <- year(data$Year)}
  if (block_length==365.25){
    bm <- data %>%
      mutate(row_index = row_number()) %>%
      group_by(year = Year) %>%
      summarise(
        max = max(!!sym(height), na.rm = TRUE),
        ind = row_index[which.max(!!sym(height))]
      ) %>%
      ungroup()
  }
  else{
    bm <- data %>%
      # Remove February 29th entries
      filter(!(month(DateTime) == 2 & day(DateTime) == 29)) %>%
      # Ensure DateTime is in Date format, not DateTime
      # Calculate the block number
      mutate(
        block = (as.numeric(as.Date(DateTime, format = "%Y-%m-%d") - ymd(Year)) %/% block_length) + 1
      ) %>%
      # Group by year and block number
      group_by(year = year(Year), block) %>%
      # Summarize to get the maximum value of F010
      summarise(max = max(!!sym(height), na.rm = TRUE), .groups = 'drop') %>%
      ungroup()
  }
  return(bm)
}



# BM estimation plot different heights with CI


BM_TIplot <- function(data,Nbs=c(500),heights=NULL, title=NULL,parameter='tail'){
  if (is.null(heights)) {heights = colnames(data)[grepl("^[Ffw]", names(data))]}
  
  datasets <- list()
  data_desc <- paste(deparse(substitute(data)),'|',difftime(data$DateTime[2],data$DateTime[1],units = 'min'),'min')
  
  
  for (Nb in Nbs) {
    outs <- list()
    for (height in heights) {
      print(height)
      # Bootstrap
      tmp <- bootstrap(data, Nb=Nb, column=height, method='BM',data_desc = data_desc)
      outs[[as.character(height)]] <- tmp$df[[parameter]]
      if (is.null(title)){
        title <- tmp$desc$fulltext
      }
    }
    datasets[[as.character(Nb)]] <- data.frame(outs)
  }
  
  
  # Colors for each Nb
  if (length(Nbs)>2){
    colors <- brewer.pal(n = length(Nbs), name = "Set1")
  }else{colors = c('red','blue')}
  
  offset <- 0.15
  alpha <- 0.05
  z <- qnorm(1 - alpha / 2)
  
  
  tailindex <- c()
  
  for (height in heights){
    tailindex <- c(tailindex,BM_fit(BM_select(data,height = height),plot=F)[parameter])
  }
  
  # Set up plot area with extra space on the right for the legend
  par(mar = c(5, 4, 4, 15) + 0.1)
  
  # Initialize plot
  plot(x = 1:length(heights), y = tailindex, col = "blue", pch = 19, cex = 1,
       xlab = "heights", ylim = c(min(unlist(lapply(datasets, function(x) tailindex - z * apply(x, 2, sd)))), 
                                  max(unlist(lapply(datasets, function(x) tailindex + z * apply(x, 2, sd))))),
       xlim = c(1-(length(Nbs)-1)*offset,length(heights)),
       xaxt = "n",ylab=parameter)
  axis(1, at = 1:length(heights), labels = heights, tick = TRUE)
  
  # Loop over datasets
  for (i in seq_along(Nbs)) {
    Nb <- Nbs[i]
    dataset <- datasets[[as.character(Nb)]]
    color <- colors[i]
    
    # Calculate the 95% confidence interval
    lower_bounds <- tailindex - z * apply(dataset, 2, sd)
    upper_bounds <- tailindex + z * apply(dataset, 2, sd)
    
    # Add points and confidence intervals
    arrows(x0 = 1:ncol(dataset) - offset * (i-1), y0 = lower_bounds, x1 = 1:ncol(dataset) - offset * (i-1), 
           y1 = upper_bounds, code = 3, angle = 90, length = 0.1, col = color)
  }
  
  # Add legend outside the plot
  legend("topright", inset = c(-0.8, 0), legend = c("Value without bootstrap", paste("Nb =", Nbs)), 
         col = c("blue", colors), pch = c(19, rep(1, length(Nbs))), 
         bty = 'n', xpd = TRUE)
  title(title)
  
  return(datasets)
}

