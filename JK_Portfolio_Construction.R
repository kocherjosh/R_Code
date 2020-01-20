# Local_Returns <- returns[,Symbols]
# Freq = "months"
Equal_Weight <- function(Local_Returns,Freq = "months") {
  
  
  EW_Wts <- function(r) {
    x <- ifelse(is.na(r),0,1)
    return(x/sum(x, na.rm = T))
  }
  
  EW_Weights_Xts <- xts(x = t(apply(Local_Returns,1, FUN = EW_Wts)), order.by = index(Local_Returns))
  EW_Weights_Xts[is.na(EW_Weights_Xts)] <- 0
  
  print(paste0("Weights at each time point are correct: ",sum(rowSums(EW_Weights_Xts)) == nrow(EW_Weights_Xts)))
  
  Local_Returns[is.na(Local_Returns)] <- 0
  EW_Weights_Xts <- EW_Weights_Xts[endpoints(EW_Weights_Xts, on = Freq)]
  Port_EW <- Return.portfolio(R = Local_Returns, weights = EW_Weights_Xts, 
                              verbose = TRUE, rebalance_on = "months")
  
  charts.PerformanceSummary(Port_EW$returns, main = "Equal Weight")
  
  print(SummaryStats(Port_EW$returns))
  apply(Port_EW$contribution, 2, sum)
  
  return(Port_EW)
}

# end section 4.1

# 4.2 Momentum
# momentum weighted by cross sectional momentum itself
# Local_Returns <- returns[,Symbols]
# Freq <- c("months")
# Lookback <- 6

Momentum <- function(Local_Returns, Lookback = 6, Freq = "months") {
  RefDates <- index(Local_Returns)
  RefDates2 <- RefDates[endpoints(RefDates, on = Freq)]
  i <- 7
  tank <- list()
  for (i in (Lookback+1):length(RefDates2)) {
    Train <- Local_Returns[paste0(RefDates2[i-Lookback],"/",RefDates2[i])]
    tank[[i-Lookback]] <- data.frame(Date = RefDates2[i], t(apply(Train, 2, function(c) {
      return((prod(1+c, na.rm = T)-1))})))
  }
  Rolling_Returns_Xts <- do.call(rbind, tank)
  Rolling_Returns_Xts <- xts(x = Rolling_Returns_Xts[,-1], order.by = Rolling_Returns_Xts$Date)
  # Rolling_Returns_Xts[Rolling_Returns_Xts == 0] <- NA
  
  Mom_Wts <- function(r) {
    x <- (r - min(r, na.rm = T))
    wts <- x/sum(x, na.rm = T)
    return(wts)
  }
  
  Momentum_Weights <- xts(order.by = index(Rolling_Returns_Xts),
                          x = t(apply(Rolling_Returns_Xts, 1, function (r){
                            Mom_Wts(r)
                          })))
  
  print(paste0("Weights at each time point are correct: ",sum(rowSums(Momentum_Weights)) == nrow(Momentum_Weights)))
  start(Local_Returns)
  start(Momentum_Weights)
  Local_Returns[is.na(Local_Returns)] <- 0
  Momentum_Weights[is.na(Momentum_Weights)] <- 0
  
  Port_Momentum <- Return.portfolio(R = Local_Returns[paste0(start(Momentum_Weights),"/")],
                                    weights = Momentum_Weights, verbose = T, rebalance_on = "months")
  
  charts.PerformanceSummary(Port_Momentum$returns, main = "Momentum")
  print(SummaryStats(Port_Momentum$returns))
  apply(Port_Momentum$contribution, 2, sum)
  return(Port_Momentum)
}

# end section 4.2

# 4.3 Volatility inverse weight
# Local_Returns <- returns[,Symbols]
# Freq <- c("months")
# Lookback <- 12

Volatility <- function(Local_Returns, Lookback = 6, Freq = "months") {
  RefDates <- index(Local_Returns)
  RefDates2 <- RefDates[endpoints(RefDates, on = Freq)]
  i <- 7
  tank <- list()
  for (i in (Lookback+1):length(RefDates2)) {
    Train <- Local_Returns[paste0(RefDates2[i-Lookback],"/",RefDates2[i])]
    tank[[i-Lookback]] <- data.frame(Date = RefDates2[i], t(apply(Train, 2, function(c) {
      return((sd(c, na.rm = T)))})))
  }
  Rolling_Volatility_Xts <- do.call(rbind, tank)
  Rolling_Volatility_Xts <- xts(x = Rolling_Volatility_Xts[,-1], order.by = Rolling_Volatility_Xts$Date)
  Rolling_Volatility_Xts[Rolling_Volatility_Xts == 0] <- NA
  
  Vol_Wts <- function(r) {
    x <- 1/r
    wts <- x / sum(x, na.rm = T)
    return(wts)
  }
  
  InverseVol_Weights <- xts(order.by =  index(Rolling_Volatility_Xts),
                            x = t(apply(Rolling_Volatility_Xts,1, function(r) {
                              Vol_Wts(r)
                            })))
  InverseVol_Weights[is.na(InverseVol_Weights)] <- 0
  rowSums(InverseVol_Weights)
  
  print(paste0("Weights at each time point are correct: ",sum(rowSums(InverseVol_Weights)) == nrow(InverseVol_Weights)))
  
  Local_Returns[is.na(Local_Returns)] <- 0
  
  Port_InverseVol <- Return.portfolio(R = Local_Returns[paste0(start(InverseVol_Weights),"/")],
                                      weights = InverseVol_Weights,
                                      verbose = T,
                                      rebalance_on = "months")
  
  charts.PerformanceSummary(Port_InverseVol$returns, main = "Inverse Vol")
  
  print(SummaryStats(Port_InverseVol$returns))
  apply(Port_InverseVol$contribution,2,sum)
  return(Port_InverseVol) 
}



# end section 4.3

# Section 4.4 Risk Parity
# Local_Returns <- returns[,Symbols]
# Freq <- c("months")

RiskParity <- function(Local_Returns, Lookback = 6) {
  # lookback refers to monthly frequency  
  eval_f <- function(w,cov.mat,vol.target) {
    vol <- sqrt(as.numeric(t(w) %*% cov.mat %*% w))
    marginal.contribution <- cov.mat %*% w / vol
    return( sum((vol/length(w) - w * marginal.contribution)^2) )
  }
  
  # numerical gradient approximation for solver
  eval_grad_f <- function(w,cov.mat,vol.target) {
    out <- w
    for (i in 0:length(w)) {
      up <- dn <- w
      up[i] <- up[i]+.0001
      dn[i] <- dn[i]-.0001
      out[i] = (eval_f(up,cov.mat=cov.mat,vol.target=vol.target) - eval_f(dn,cov.mat=cov.mat,vol.target=vol.target))/.0002
    }
    return(out)
  }
  
  RiskParityCalc <- function(r) {
    std <- apply(r,2,sd)
    cov.mat <- cov(r)
    x0 <- 1/std/sum(1/std)
    Cnames <- names(std)
    res <- nloptr( x0=x0,
                   eval_f=eval_f,
                   eval_grad_f=eval_grad_f,
                   eval_g_eq=function(w,cov.mat,vol.target) { sum(w) - 1 },
                   eval_jac_g_eq=function(w,cov.mat,vol.target) { rep(1,length(std)) },
                   lb=rep(0,length(std)),ub=rep(1,length(std)),
                   opts = list("algorithm"="NLOPT_LD_SLSQP","print_level" = 3,"xtol_rel"=1.0e-8,"maxeval" = 1000),
                   cov.mat = cov.mat,vol.target=.10 )
    names(res$solution) <- Cnames
    # print(res$solution)
    dim(res)
    return(res$solution)
  }
  
  tank <- list()
  i <- 206
  Local_Returns[is.na(Local_Returns)] <- 0
  RefDates <- index(Local_Returns)
  RefDates2 <- RefDates[endpoints(RefDates, on = Freq)]
  
  for (i in ((Lookback+1):length(RefDates2))) {
    StartDate <- index(Local_Returns)[i-Lookback]
    EndDate <- index(Local_Returns)[i]
    # Return_subset <- Local_Returns[paste0(StartDate,"/",EndDate)]
    
    Return_subset <- Local_Returns[paste0(RefDates2[i-Lookback],"/",RefDates2[i])]
    
    # capture the original column order
    Cnames <- names(Return_subset)
    # find columns with no returns
    Zeros <- data.frame(Values = which(colSums(Return_subset) == 0))
    if (nrow(Zeros) > 0 ) {
      Zeros$Ticker <- row.names(Zeros)
      Zeros$Values <- 0
      Zeros <- Zeros[,c("Ticker","Values")]
      Return_subset2 <- Return_subset[,-which(colSums(Return_subset) == 0)]
      dummy <- data.frame(t(RiskParityCalc(Return_subset2)))
      dummy <- gather(dummy,key = "Ticker",value = "Values")
      dummy2 <- rbind(Zeros,dummy)
      
    } else {
      Return_subset2 <- Return_subset
      dummy <- data.frame(t(RiskParityCalc(Return_subset2)))
      dummy <- gather(dummy,key = "Ticker",value = "Values")
      dummy2 <- dummy
    }
    
    # run risk parity script for columns with returns
    
    # combine with zero return columns
    dummy2_wide <- spread(dummy2,Ticker, Values)
    dummy2_wide <- dummy2_wide[,Cnames]
    
    tank[[i-Lookback]] <- data.frame(Date = EndDate, dummy2_wide)
  }
  
  RiskParity_Weights <- do.call(rbind, tank)
  RiskParity_Weights <- xts(x = RiskParity_Weights[,-1], order.by = RiskParity_Weights$Date)
  apply(RiskParity_Weights, 2, function(c) mean(c, na.rm = T))
  RiskParity_Weights[is.na(RiskParity_Weights)] <- 0
  Port_RiskParity <- Return.portfolio(R = Local_Returns[paste0(start(RiskParity_Weights),"/")],
                                      weights = RiskParity_Weights,
                                      verbose = T,
                                      rebalance_on = "months")
  
  # charts.PerformanceSummary(Port_RiskParity$returns, main = "Risk Parity")
  # print(SummaryStats(Port_RiskParity$returns))
  
  return(Port_RiskParity)
}


# end of section 4.4
# 4.5 Momentum TopN
# momentum, take top n, weight equally
# Local_Returns <- returns[,Symbols]
# Freq <- c("months")
# Lookback <- 6
# N <- 3

Momentum_TopN <- function(Local_Returns, Lookback = 6, Freq = "months",N = 3) {
  
  RefDates <- index(Local_Returns)
  RefDates2 <- RefDates[endpoints(RefDates, on = Freq)]
  i <- 7
  tank <- list()
  for (i in (Lookback+1):length(RefDates2)) {
    Train <- Local_Returns[paste0(RefDates2[i-Lookback],"/",RefDates2[i])]
    tank[[i-Lookback]] <- data.frame(Date = RefDates2[i], t(apply(Train, 2, function(c) {
      return((prod(1+c, na.rm = T)-1))})))
  }
  Rolling_Returns_Xts <- do.call(rbind, tank)
  Rolling_Returns_Xts <- xts(x = Rolling_Returns_Xts[,-1], order.by = Rolling_Returns_Xts$Date)
  
  Mom_Wts <- function(r) {
    # x <- (r - min(r, na.rm = T))
    # wts <- x/sum(x, na.rm = T)
    x <- rank(x= -r,ties.method = "random")
    x <- ifelse(x <= N, 1/N,0)
    return(x)
  }
  
  Momentum_Weights <- xts(order.by = index(Rolling_Returns_Xts),
                          x = t(apply(Rolling_Returns_Xts, 1, function (r){
                            Mom_Wts(r)
                          })))
  
  print(paste0("Weights at each time point are correct: ",sum(rowSums(Momentum_Weights)) == nrow(Momentum_Weights)))
  start(Local_Returns)
  start(Momentum_Weights)
  Local_Returns[is.na(Local_Returns)] <- 0
  Momentum_Weights[is.na(Momentum_Weights)] <- 0
  
  Port_Momentum_TopN <- Return.portfolio(R = Local_Returns[paste0(start(Momentum_Weights),"/")],
                                         weights = Momentum_Weights, verbose = T, rebalance_on = "months")
  
  # charts.PerformanceSummary(Port_Momentum_TopN$returns, main = "Momentum")
  # print(SummaryStats(Port_Momentum_TopN$returns))
  
  return(Port_Momentum_TopN)
}

# end section 4.5

# 4.6 Sharpe Weighted
# risk adjusted  weighted by cross sectional risk adjusted ratio
# Local_Returns <- returns[,Symbols]
# Freq <- c("months")
# Lookback <- 6

Sharpe_TopN <- function(Local_Returns, Lookback = 6, Freq = "months",N = 3) {
  
  RefDates <- index(Local_Returns)
  RefDates2 <- RefDates[endpoints(RefDates, on = Freq)]
  i <- 7
  tank <- list()
  for (i in (Lookback+1):length(RefDates2)) {
    Train <- Local_Returns[paste0(RefDates2[i-Lookback],"/",RefDates2[i])]
    tank[[i-Lookback]] <- data.frame(Date = RefDates2[i], t(apply(Train, 2, function(c) {
      return((prod(1+c, na.rm = T)-1)/(sd(c, na.rm = T)*sqrt(252)))})))
  }
  Rolling_Sharpe_Xts <- do.call(rbind, tank)
  Rolling_Sharpe_Xts <- xts(x = Rolling_Sharpe_Xts[,-1], order.by = Rolling_Sharpe_Xts$Date)
  
  Sharpe_Wts <- function(r) {
    # x <- (r - min(r, na.rm = T))
    # wts <- x/sum(x, na.rm = T)
    x <- rank(x= -r,ties.method = "random")
    x <- ifelse(x <= N, 1/N,0)
    return(x)
  }
  Sharpe_Weights <- xts(order.by = index(Rolling_Sharpe_Xts),
                        x = t(apply(Rolling_Sharpe_Xts, 1, function (r){
                          Sharpe_Wts(r)
                        })))
  rowSums(Sharpe_Weights)
  
  print(paste0("Weights at each time point are correct: ",sum(rowSums(Sharpe_Weights)) == nrow(Sharpe_Weights)))
  start(Local_Returns)
  start(Sharpe_Weights)
  Local_Returns[is.na(Local_Returns)] <- 0
  Sharpe_Weights[is.na(Sharpe_Weights)] <- 0
  
  Port_Sharpe <- Return.portfolio(R = Local_Returns[paste0(start(Sharpe_Weights),"/")],
                                  weights = Sharpe_Weights, verbose = T, rebalance_on = "months")
  
  # charts.PerformanceSummary(Port_Sharpe$returns, main = "Sharpe")
  # print(SummaryStats(Port_Sharpe$returns))
  # print(paste0("Contribution: ",apply(Port_Sharpe$contribution, 2, sum)))
  return(Port_Sharpe)
}

# end section 4.6

# 4.7 Downside volatility weighted
# Local_Returns <- returns[,Symbols]
# Freq <- c("months")
# Lookback <- 6

Semi_Variance_Weighted <- function(Local_Returns, Lookback = 6, Freq = "months") {
  # this roll apply is not calculating correctly
  RefDates <- index(Local_Returns)
  RefDates2 <- RefDates[endpoints(RefDates, on = Freq)]
  i <- 7
  tank <- list()
  for (i in (Lookback+1):length(RefDates2)) {
    Train <- Local_Returns[paste0(RefDates2[i-Lookback],"/",RefDates2[i])]
    tank[[i-Lookback]] <- data.frame(Date = RefDates2[i], SemiVariance(Train))
    
  }
  
  Rolling_Semi_Xts <- do.call(rbind, tank)
  Rolling_Semi_Xts <- xts(x = Rolling_Semi_Xts[,-1], order.by = Rolling_Semi_Xts$Date)
  Rolling_Semi_Xts[Rolling_Semi_Xts == 0] <- NA
  
  Vol_Wts <- function(r) {
    x <- 1/r
    wts <- x / sum(x, na.rm = T)
    return(wts)
  }
  
  InverseVol_Weights <- xts(order.by =  index(Rolling_Semi_Xts),
                            x = t(apply(Rolling_Semi_Xts,1, function(r) {
                              Vol_Wts(r)
                            })))
  InverseVol_Weights[is.na(InverseVol_Weights)] <- 0
  rowSums(InverseVol_Weights)
  apply(InverseVol_Weights,2,mean)
  write.csv(x = data.frame(Date = index(InverseVol_Weights),coredata(InverseVol_Weights), stringsAsFactors = F),
            file = "Josh_Kocher_Weight_Signal_File.csv", row.names = F)
  
  print(paste0("Weights at each time point are correct: ",sum(rowSums(InverseVol_Weights)) == nrow(InverseVol_Weights)))
  
  Local_Returns[is.na(Local_Returns)] <- 0
  
  Port_Semi <- Return.portfolio(R = Local_Returns[paste0(start(InverseVol_Weights),"/")],
                                weights = InverseVol_Weights,
                                verbose = T,
                                rebalance_on = "months")
  
  charts.PerformanceSummary(Port_Semi$returns, main = "Inverse Downside Variance")
  
  print(SummaryStats(Port_Semi$returns))
  apply(Port_Semi$contribution,2,sum)
  
  Port <- Port_Semi
  # print("Summary Stats:")
  # print(SummaryStats(Port_Semi$returns))
  # print("Contribution:")
  # print(apply(Port_Semi$contribution, 2, function(s) sum(s, na.rm = T)))
  # print("Average Weight:")
  # print(apply(Port_Semi$BOP.Weight, 2, function(s) mean(s, na.rm = T)))
  # print("SD Weight:")
  # print(apply(Port_Semi$BOP.Weight, 2, function(s) sd(s, na.rm = T)))
  # print("Risk adjusted contribution:")
  # print(apply(Port_Semi$contribution, 2, function(c) {
  #   mean(c, na.rm = T) / sd(c, na.rm = T)
  # }))
  # 
  beginWeights <- Port$BOP.Weight[endpoints(Port$BOP.Weight, on = "months"),]
  endWeights <- Port$EOP.Weight[endpoints(Port$EOP.Weight, on = "months"),]

  # txns <- beginWeights - lag(endWeights)
  # monthlyTO <- xts(rowSums(abs(txns[,1:9])), order.by=index(txns))
  # plot(monthlyTO)
  
  
  
  return(Port_Semi)
}

# end section 4.6

# Local_Returns <- returns
# Freq <- c("months")
# Lookback <- 6

Sharpe_Weighted <- function(Local_Returns, Lookback = 6, Freq = "months") {
  RefDates <- index(Local_Returns)
  RefDates2 <- RefDates[endpoints(RefDates, on = Freq)]
  i <- 7
  tank <- list()
  for (i in (Lookback+1):length(RefDates2)) {
    Train <- Local_Returns[paste0(RefDates2[i-Lookback],"/",RefDates2[i])]
    tank[[i-Lookback]] <- data.frame(Date = RefDates2[i], t(apply(Train, 2, function(c) {
      return( (prod(1+c, na.rm = T)-1) / (sd(c, na.rm = T)*sqrt(252)) )})))
  }
  
  Rolling_Sharpe_Xts <- do.call(rbind, tank)
  Rolling_Sharpe_Xts <- xts(x = Rolling_Sharpe_Xts[,-1], order.by = Rolling_Sharpe_Xts$Date)
  
  Shift1 <- xts(order.by = index(Rolling_Sharpe_Xts),
                x = t(apply(Rolling_Sharpe_Xts, 1, function (r){
                  return(r-min(r,na.rm = T))
                })))
  Sharpe_Weights2 <- xts(order.by = index(Shift1),
                         x = t(apply(Shift1, 1, function (r){
                           return(r / sum(r, na.rm = T))
                         })))
  # Sharpe_Weights2 <- xts(order.by = index(Shift1),
  #                        x = t(apply(Shift1, 1, function (r){
  #                          return(exp(r) / sum(exp(r), na.rm = T))
  #                        })))
  
  # Sharpe_Wts2 <- function(r) {
  #   x <- (r - min(r, na.rm = T))
  #   wts <- x/sum(x, na.rm = T)
  #   # x <- rank(x= -r,ties.method = "random")
  #   # x <- ifelse(x <= N, 1/N,0)
  #   return(x)
  # }
  # Sharpe_Weights2 <- xts(order.by = index(Rolling_Sharpe_Xts),
  #                       x = t(apply(Rolling_Sharpe_Xts, 1, function (r){
  #                         Sharpe_Wts2(r)
  #                       })))
  rowSums(Sharpe_Weights2)
  
  print(paste0("Weights at each time point are correct: ",sum(rowSums(Sharpe_Weights2)) == nrow(Sharpe_Weights2)))
  start(Local_Returns)
  start(Sharpe_Weights2)
  Local_Returns[is.na(Local_Returns)] <- 0
  Sharpe_Weights2[is.na(Sharpe_Weights2)] <- 0
  
  Port_Sharpe <- Return.portfolio(R = Local_Returns[paste0(start(Sharpe_Weights2),"/")],
                                  weights = Sharpe_Weights2, verbose = T, rebalance_on = "months")
  
  charts.PerformanceSummary(Port_Sharpe$returns, main = "Sharpe")
  print(SummaryStats(Port_Sharpe$returns))
  print("Contribution:")
  print(apply(Port_Sharpe$contribution, 2, sum))
  print("Average Weight:")
  print(apply(Port_Sharpe$BOP.Weight, 2, mean))
  return(Port_Sharpe) 
}

# 4.3 Volatility inverse weight with crash

# Local_Returns <- returns
# Freq <- c("months")
# Lookback <- 6

Volatility_Crash <- function(Local_Returns, Lookback = 6, Freq = "months") {
  
  RefDates <- index(Local_Returns)
  RefDates2 <- RefDates[endpoints(RefDates, on = Freq)]
  i <- 7
  VolTank <- list()
  RetTank <- list()
  
  for (i in (Lookback+1):length(RefDates2)) {
    Train <- Local_Returns[paste0(RefDates2[i-Lookback],"/",RefDates2[i])]
    VolTank[[i-Lookback]] <- data.frame(Date = RefDates2[i], t(apply(Train, 2, function(c) {
      return((sd(c, na.rm = T)))})))
    RetTank[[i-Lookback]] <- data.frame(Date = RefDates2[i], t(apply(Train, 2, function(c) {
      return((prod(1+c, na.rm = T)-1))})))
  }
  
  Rolling_Volatility <- do.call(rbind, VolTank)
  Rolling_Volatility_Xts <- xts(x = Rolling_Volatility[,-1], order.by = Rolling_Volatility$Date)
  Rolling_Volatility_Xts[Rolling_Volatility_Xts == 0] <- NA
  
  Rolling_Returns <- do.call(rbind, RetTank)
  # Rolling_Returns_Xts <- xts(x = Rolling_Returns_Xts[,-1], order.by = Rolling_Returns_Xts$Date)
  
  
  Rolling_Volatility_Long <- gather(Rolling_Volatility, key = "Asset",value = "Volatility",-Date)
  
  Rolling_Returns_Long <- gather(Rolling_Returns, key = "Asset",value = "Return",-Date)
  
  
  Vol_Wts <- function(r) {
    x <- 1/r
    wts <- x / sum(x, na.rm = T)
    return(wts)
  }
  
  InverseVol_Weights <- xts(order.by =  index(Rolling_Volatility_Xts),
                            x = t(apply(Rolling_Volatility_Xts,1, function(r) {
                              Vol_Wts(r)
                            })))
  InverseVol_Weights[is.na(InverseVol_Weights)] <- 0
  rowSums(InverseVol_Weights)
  
  InverseVol_Weights_Long <- data.frame(Date = index(InverseVol_Weights),coredata(InverseVol_Weights))
  InverseVol_Weights_Long <- gather(InverseVol_Weights_Long, key = "Asset",value = "Weight",-Date)
  InverseVol_Weights_Long %>% group_by(Date) %>% summarize(Total = sum(Weight, na.rm = T))
  
  All_Data <- merge(InverseVol_Weights_Long,Rolling_Returns_Long, by = c("Date","Asset"))
  All_Data$Adjusted_Weight <- ifelse(All_Data$Return <0 , 0 , 1)
  All_Data$Adjusted_Weight <- All_Data$Adjusted_Weight * All_Data$Weight
  
  Risky_Weights <- All_Data %>% group_by(Date) %>% summarize(Total = sum(Adjusted_Weight, na.rm = T))
  Risky_Weights$Cash <- 1-Risky_Weights$Total
  
  Cash <- xts(x = Risky_Weights$Cash, order.by = Risky_Weights$Date)
  
  InverseVol_Weights2 <- spread(All_Data[,c("Date","Asset","Adjusted_Weight")],Asset, Adjusted_Weight)
  InverseVol_Weights_Xts <- xts(x = InverseVol_Weights2[,-1], order.by = InverseVol_Weights2$Date)
  
  All_Weights_Xts <- merge.xts(InverseVol_Weights_Xts, Cash)
  
  
  
  print(paste0("Weights at each time point are correct: ",sum(rowSums(All_Weights_Xts)) == nrow(All_Weights_Xts)))
  
  Local_Returns[is.na(Local_Returns)] <- 0
  Local_Returns$Cash <- rep(0,length(index(Local_Returns)))
  
  # Local_Returns2 <- merge.xts(Local_Returns, Cash_Returns)
  
  All_Weights_Xts <- All_Weights_Xts[,names(Local_Returns)]
  names(Local_Returns) == names(All_Weights_Xts)
  start(Local_Returns)
  start(All_Weights_Xts)
  
  Port_InverseVol_Crash <- Return.portfolio(R = Local_Returns[paste0(start(All_Weights_Xts),"/")],
                                            weights = All_Weights_Xts,
                                            verbose = T,
                                            rebalance_on = "months")
  
  charts.PerformanceSummary(Port_InverseVol_Crash$returns, main = "Inverse Vol Crash")
  
  # Port <- Port_InverseVol_Crash
  # print("Summary Stats:")
  # print(SummaryStats(Port$returns))
  # print("Contribution:")
  # print(apply(Port$contribution, 2, sum))
  # print("Average Weight:")
  # print(apply(Port$BOP.Weight, 2, mean))
  # print("SD Weight:")
  # print(apply(Port$BOP.Weight, 2, sd))
  # print("Risk adjusted contribution:")
  # print(apply(Port$contribution, 2, function(c) {
  #   mean(c, na.rm = T) / sd(c, na.rm = T)
  # }))
  # # write.xlsx(x = data.frame(Date = index(Port_Inverse_Vol_Crash$BOP.Weight), coredata(Port_Inverse_Vol_Crash$BOP.Weight)),
  # #            file = "Inverse_Vol_Crash.xlsx",sheetName = "BOP_Weight")
  # # write.xlsx(x = data.frame(Date = index(Port_Inverse_Vol_Crash$EOP.Weight), coredata(Port_Inverse_Vol_Crash$EOP.Weight)),
  # #            file = "Inverse_Vol_Crash.xlsx",sheetName = "EOP_Weight",append = T)
  # # 
  # Weight_long <- gather(data.frame(Date = index(Port$BOP.Weight), coredata(Port$BOP.Weight)),key = "Ticker",value = "Weight",-Date)
  # ggplot(data = Weight_long, aes(x=Date,y = Weight, fill = Ticker)) + geom_area()
  # 
  # 
  # print(tail(Rolling_Returns,1))
  
  return(Port_InverseVol_Crash) 
}



Local_Returns <- returns
Freq <- c("days")
Lookback <- 126


Momentum_Crash <- function(Local_Returns, Lookback = 6, Freq = "months") {
  RefDates <- index(Local_Returns)
  RefDates2 <- RefDates[endpoints(RefDates, on = Freq)]
  i <- 7
  tank <- list()
  for (i in (Lookback+1):length(RefDates2)) {
    Train <- Local_Returns[paste0(RefDates2[i-Lookback],"/",RefDates2[i])]
    tank[[i-Lookback]] <- data.frame(Date = RefDates2[i], t(apply(Train, 2, function(c) {
      return((prod(1+c, na.rm = T)-1))})))
  }
  Rolling_Returns_Xts <- do.call(rbind, tank)
  Rolling_Returns_Xts <- xts(x = Rolling_Returns_Xts[,-1], order.by = Rolling_Returns_Xts$Date)
  # set any security with negaqtive return to zero so that the weighting function will set the weight
  # to zero as well
  
  Rolling_Returns_Xts[Rolling_Returns_Xts < 0] <- 0
  

  Mom_Wts <- function(r) {
    x <- (r - min(r, na.rm = T))
    wts <- x/sum(x, na.rm = T)
    return(wts)
  }
  
  Momentum_Weights <- xts(order.by = index(Rolling_Returns_Xts),
                          x = t(apply(Rolling_Returns_Xts, 1, function (r){
                            Mom_Wts(r)
                          })))
  
  
  Cash_Xts <- xts(x = sum(1 - rowSums(coredata(Momentum_Weights)), na.rm = T), order.by = index(Momentum_Weights))
  
  
  
  
  print(paste0("Weights at each time point are correct: ",sum(rowSums(Momentum_Weights)) == nrow(Momentum_Weights)))
  start(Local_Returns)
  start(Momentum_Weights)
  Local_Returns[is.na(Local_Returns)] <- 0
  Momentum_Weights[is.na(Momentum_Weights)] <- 0
  
  Port_Momentum_Crash <- Return.portfolio(R = Local_Returns[paste0(start(Momentum_Weights),"/")],
                                    weights = Momentum_Weights, verbose = T, rebalance_on = "months")
  
  charts.PerformanceSummary(Port_Momentum_Crash$returns, main = "Momentum")
  
  
  Port <- Port_Momentum_Crash
  print("Summary Stats:")
  print(SummaryStats(Port$returns))
  print("Contribution:")
  print(apply(Port$contribution, 2, sum))
  print("Average Weight:")
  print(apply(Port$BOP.Weight, 2, mean))
  print("SD Weight:")
  print(apply(Port$BOP.Weight, 2, sd))
  print("Risk adjusted contribution:")
  print(apply(Port$contribution, 2, function(c) {
    mean(c, na.rm = T) / sd(c, na.rm = T)
  }))
  # write.xlsx(x = data.frame(Date = index(Port_Inverse_Vol_Crash$BOP.Weight), coredata(Port_Inverse_Vol_Crash$BOP.Weight)),
  #            file = "Inverse_Vol_Crash.xlsx",sheetName = "BOP_Weight")
  # write.xlsx(x = data.frame(Date = index(Port_Inverse_Vol_Crash$EOP.Weight), coredata(Port_Inverse_Vol_Crash$EOP.Weight)),
  #            file = "Inverse_Vol_Crash.xlsx",sheetName = "EOP_Weight",append = T)
  # 
  Weight_long <- gather(data.frame(Date = index(Port$BOP.Weight), coredata(Port$BOP.Weight)),key = "Ticker",value = "Weight",-Date)
  ggplot(data = Weight_long, aes(x=Date,y = Weight, fill = Ticker)) + geom_area()
  
  
  return(Port_Momentum_Crash)
  
}
