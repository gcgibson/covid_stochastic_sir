###simulate data


library(rjags)
library(R2jags)
model.loc = ("rw_intercept.txt")
jagsscript =  cat("
                  model {
                  S[1] <- pop - infected
                  I[1]  <- infected
                  gamma ~ dgamma(1,14)
                  beta[1] ~ dbeta(1000*.2,1000 *(1-.2))
                  
                  for(i in 2:N) {
                  beta[i] ~ dnorm(beta[i-1],10000)
                  Sk1[i] <- -beta[i]/pop*S[i-1]*I[i-1]
                  Ik1[i] <- beta[i]/pop * S[i-1] * I[i-1]- gamma * I[i-1]
                  
                  Sk2[i] <- -beta[i]/pop*(S[i-1] +  .5*Sk1[i])*(I[i-1] + .5*Ik1[i])
                  Ik2[i] <- beta[i]/pop * (S[i-1] +  .5*Sk1[i]) * (I[i-1] + .5*Ik1[i]) - gamma * (I[i-1] + .5*Ik1[i])
                  
                  
                  
                  Sk3[i] <-  -beta[i]/pop*(S[i-1] + .5*Sk2[i])*(I[i-1] + .5*Ik2[i])
                  Ik3[i] <- beta[i]/pop*(S[i-1] +.5*Sk2[i])*(I[i-1] + .5*Ik2[i]) - gamma*(I[i-1] + .5*Ik2[i])
                  
                  Sk4[i] <- -beta[i]/pop*(S[i-1] + Sk3[i])*(I[i-1] + Ik3[i])
                  Ik4[i] <- beta[i]/pop*(S[i-1] + Sk3[i])*(I[i-1] + Ik3[i])-gamma*(I[i-1] + Ik3[i])
                  
                  mean_S[i] <- S[i-1] +  .16666* (Sk1[i] + 2 * Sk2[i] + 2*Sk3[i] +Sk4[i] )
                  mean_I[i] <- I[i-1] +  .166666* (Ik1[i] + 2 * Ik2[i] + 2*Ik3[i] +Ik4[i] )
                  
                  S[i]  <- mean_S[i]
                  I[i]  <- mean_I[i]
                  }
                  
                  
                  for(i in 1:N){
                  Y[i] ~dnorm(.8*I[i],100000)
                  }
                  
                  
                  }
                  ", 
                  file = model.loc)


# define location


## ILI Dataarima forecasts
library(cdcForecastUtils)
library(forecast)
state_flu_data <- download_and_preprocess_state_flu_data()
full_sub_file <- data.frame()
unique_locations <-unique(state_flu_data$region)
unique_locations <- unique_locations[unique_locations!="Florida"]
unique_locations <- unique_locations[unique_locations!="Commonwealth of the Northern Mariana Islands"]
unique_locations_sample <- sample(unique_locations,10)
for (location in unique_locations){
  tryCatch({
  current_season_sate_flu <- state_flu_data[state_flu_data$region %in% c(unique_locations_sample) & state_flu_data$season == "2019/2020",]
  ggplot(current_season_sate_flu,aes(x=season_week,y=total_patients,col="total patients")) +geom_line() +geom_line(aes(x=season_week,y=ilitotal*5,col='ilitotal*5'))+ facet_wrap(~region,scales='free')

  model_fit_numerator <- auto.arima(current_season_sate_flu$ilitotal,lambda = "auto")
  model_fit_denominator <- auto.arima(current_season_sate_flu$total_patients)
  plot(forecast(model_fit_denominator,h=6))
  ##### COVID CASE RELATED DATA
  covid_ts <- read.csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
  covid_ts_us <- covid_ts[covid_ts$Province.State == location,]
  covid_ts_us <- (colSums(covid_ts_us[,5:ncol(covid_ts_us)]))
  uspop <- 6.9*10^6
  covid_ts_us <- tail(covid_ts_us,8)
  
  jags.data = list(Y = c(covid_ts_us,rep(NA,7*27)),N=length(covid_ts_us)+7*27,infected=covid_ts_us[1],pop=uspop)
  jags.params = c("I","beta","gamma")
  mod_rw_intercept  = jags(jags.data, parameters.to.save = jags.params, 
                           model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 1, 
                           n.iter = 10000, DIC = TRUE)
  
  N <- ncol+7*27
  current_time <- length(covid_ts_us)
  library(ggplot2)
  #df_plot <- data.frame(y=c(t(mod_rw_intercept$BUGSoutput$sims.matrix[1:1000,paste0("I[",2:N,"]")])),x=rep(2:N,1000))
  #ggplot(df_plot,aes(x=x,y=y,group=1))+ geom_point(alpha=.1) + geom_line(data=data.frame(x=1:(length(covid_ts_us)),y=covid_ts_us),aes(x=x,y=y,col='truth'))
  
  ## set up submission 
  nsim <- 1000
  cdc_report_ew <- get_current_date_from_flu_data(state_flu_data)
  
  trajectory_matirx <- matrix(NA,nrow=nsim,ncol=27)
  
  for (h in 1:27){
    sequence_by_7 <- seq(1,28*7,by=7)
    # 1 week ahead
    k_wk_ahead_df <- df_plot[df_plot$x %in% seq(current_time+sequence_by_7[h],current_time+sequence_by_7[h+1]-1),]
    k_wk_ahead_sim_num <- rep(NA,1000)
    k_wk_ahead_sim_denom <- rep(NA,1000)
    
    for (i in 1:1000) {
      try(k_wk_ahead_sim_num[i] <- sarimaTD:::simulate.Arima(model_fit_numerator, nsim = h)[h], silent = TRUE)
      try(k_wk_ahead_sim_denom[i] <- sarimaTD:::simulate.Arima(model_fit_denominator, nsim = h)[h], silent = TRUE)
      
    }
    if (h <= 6){
      trajectory_matirx[,h] <- 100*(k_wk_ahead_sim_num+k_wk_ahead_sim_denom/uspop*k_wk_ahead_df$y[1:nsim])/k_wk_ahead_sim_denom
    } else{
      trajectory_matirx[,h] <- 100*(k_wk_ahead_sim_num+k_wk_ahead_sim_denom/uspop*k_wk_ahead_df$y[1:nsim])/k_wk_ahead_sim_denom
      
    }
    }
  
  library(cdcForecastUtils)
  traj_w_observed_data <- cbind(matrix(rep(1,nsim*1),nrow=nsim,ncol=1),trajectory_matirx)
  traj_w_observed_data[traj_w_observed_data>100] <- 100
  traj_w_observed_data[traj_w_observed_data <0] <-0
  traj <- trajectories_to_binned_distributions(traj_w_observed_data,targets=c("wk ahead"),
                                               season_start_ew = season_start_ew,
                                               season_end_ew =season_end_ew,
                                               cdc_report_ew = "2020-ew10",
                                               h_max=6,
                                               bins=   c(seq(0, 25, by = .1), 100))
  traj_list_inp <- list()
  traj_list_inp[[1]] <- traj_w_observed_data
  traj_list <- data.frame(location=location,trajectories=I(traj_list_inp))
  traj$location <- location
  
  full_sub_file <- rbind(full_sub_file,traj)
  })
  #get_viz_from_submission_df(traj)
  
  #plot_trajectories_and_intervals(flu_data = state_flu_data,
                      #            submission = traj,
                     #             trajectories_by_location = traj_list,
                       #           season_start_ew = "2020-ew10",
                       #           season_end_ew = "2020-ew35",
                       #           target_variable = "unweighted_ili",
                        #          cdc_report_ew = cdc_report_ew)
}
