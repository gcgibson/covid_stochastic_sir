library(R2jags)
model.loc = ("rw_intercept.txt")
jagsscript =  cat("
        model {

                  
                  gamma_tot ~ dgamma(1,7)
                  beta_tot ~ dbeta(1000*.2,1000 *(1-.2))
                  for (j in 1:3){
                    gamma[j] ~ dbeta(10000*gamma_tot,10000*(1-gamma_tot))
                    beta[j] ~ dbeta(10000*beta_tot,10000*(1-beta_tot))
                    

                    S[1,j] <- pop[j] - round(infected[j])
                    I[1,j]  <- infected[j]
                    
                    
                    for(i in 2:N) {
                      num_contact[i,j] ~ dbinom(min(I[i-1,j]/(max(S[i-1,j],1)),.9999), max(S[i-1,j],1))
                      num_transmit[i,j] ~ dbinom(beta[j],num_contact[i,j])
                      num_recovered[i,j] ~ dbinom(gamma[j],I[i-1,j])
                      I[i,j] <- I[i-1,j] + num_transmit[i,j] - num_recovered[i,j]
                      S[i,j] <- S[i-1,j] - num_transmit[i,j]
  
                    }

                    I_obs[1,j] <- I[1,j]* num_presented_to_network[1,j]/pop[j]
                    Y[1,j] ~dnorm(I_obs[1,j],100000)

                    for(i in 2:N){
                      num_presented_to_network[i,j] ~ dnorm(num_presented_to_network[i-1,j],1000000000)
                      I_obs[i,j] <-  I[i,j]*num_presented_to_network[i,j]/pop[j]
                      Y[i,j] ~dnorm(I_obs[i,j],100000)
                    }
                    
                  }
      }
                  ", 
                  file = model.loc)


# define location


## ILI Dataarima forecasts
library(cdcForecastUtils)
library(forecast)
library(stringr)
census_2010 <- read.csv("/Users/gcgibson/census_2010.csv")
location <- c("New York","Massachusetts","Connecticut")
    
  ##### COVID CASE RELATED DATA
    covid_ts <- read.csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
    covid_ts_us <- covid_ts[covid_ts$Province.State %in% location,]
    covid_ts_us <- unname(covid_ts_us)
    covid_ts_us <- covid_ts_us[,(ncol(covid_ts_us)-9):ncol(covid_ts_us)]
    rownames(covid_ts_us) <- NULL
    covid_ts_us <- matrix(unlist(c(covid_ts_us)),nrow=3)
    uspop <-census_2010[census_2010$state %in%str_replace(location," ","_"),]$year_2010
    h <- 3*7
    covid_ts_us <- cbind(covid_ts_us,matrix(rep(NA,h),nrow=3))

    jags.data = list(Y = t(covid_ts_us),N=ncol(covid_ts_us),infected=covid_ts_us[,1],pop=uspop)
    jags.params = c("I","beta","gamma","I_obs")
    mod_rw_intercept  = jags(jags.data, parameters.to.save = jags.params, 
                             model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 1, 
                             n.iter = 10000, DIC = TRUE)
    N <- ncol(covid_ts_us)
    library(ggplot2)
    library(jagstools)
    I_res <- mod_rw_intercept$BUGSoutput$sims.matrix[1:1000,1:N]
    
    df_plot <- data.frame(y=c(t(I_res)),x=rep(1:N,1000))
    ggplot(df_plot,aes(x=x,y=y,group=1))+ geom_point(alpha=.1) + geom_line(data=data.frame(x=1:(ncol(covid_ts_us)),y=unname(covid_ts_us[1,])),aes(x=x,y=y,col='truth')) 
    
    I_res_2 <- mod_rw_intercept$BUGSoutput$sims.matrix[1:1000,(N+1):(2*N)]
    df_plot <- data.frame(y=c(t(I_res_2)),x=rep(1:N,1000))
    ggplot(df_plot,aes(x=x,y=y,group=1))+ geom_point(alpha=.1) + geom_line(data=data.frame(x=1:(ncol(covid_ts_us)),y=unname(covid_ts_us[2,])),aes(x=x,y=y,col='truth')) 
    
    I_res_3 <- mod_rw_intercept$BUGSoutput$sims.matrix[1:1000,(2*N+1):(3*N)]
    df_plot <- data.frame(y=c(t(I_res_3)),x=rep(1:N,1000))
    ggplot(df_plot,aes(x=x,y=y,group=1))+ geom_point(alpha=.1) + geom_line(data=data.frame(x=1:(ncol(covid_ts_us)),y=unname(covid_ts_us[3,])),aes(x=x,y=y,col='truth')) 
    
  #### ILI data 
  current_season_sate_flu <- state_flu_data[state_flu_data$region %in% location & state_flu_data$season == "2019/2020" & state_flu_data$season_week <= 32,]$ilitotal
  current_season_sate_flu_matrix <- matrix(current_season_sate_flu,nrow=length(location),byrow=F)
  
  current_season_sate_flu_patients <- state_flu_data[state_flu_data$region %in% location & state_flu_data$season == "2019/2020" & state_flu_data$season_week <= 32,]$total_patients
  current_season_sate_flu_matrix_patients <- matrix(current_season_sate_flu_patients,nrow=length(location),byrow=F)
  ggplot( state_flu_data[state_flu_data$region %in% location & state_flu_data$season == "2019/2020" & state_flu_data$season_week <= 32,],aes(x=season_week,y=ilitotal)) + geom_line() +facet_wrap(~region)
  h <- 3*7
  current_season_sate_flu_matrix <- cbind(current_season_sate_flu_matrix,matrix(rep(NA,h),nrow=3))
  current_season_sate_flu_matrix_patients <- cbind(current_season_sate_flu_matrix_patients,matrix(rep(NA,h),nrow=3))
  
  jags.data = list(Y = t(current_season_sate_flu_matrix),N=ncol(current_season_sate_flu_matrix),infected=current_season_sate_flu_matrix[,1],pop=round(uspop),num_presented_to_network=t(current_season_sate_flu_matrix_patients))
  jags.params = c("I","beta","gamma","I_obs")
  mod_rw_intercept  = jags(jags.data, parameters.to.save = jags.params, 
                           model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 1, 
                           n.iter = 10000, DIC = TRUE)
  N <- ncol(current_season_sate_flu_matrix)
  I_res <- mod_rw_intercept$BUGSoutput$sims.matrix[1:1000,118:(118+N-1)]
  
  df_plot <- data.frame(y=c(t(I_res)),x=rep(1:N,1000))
  ggplot(df_plot,aes(x=x,y=y,group=1))+ geom_point(alpha=.1) + geom_line(data=data.frame(x=1:(ncol(current_season_sate_flu_matrix)),y=unname(current_season_sate_flu_matrix[1,])),aes(x=x,y=y,col='truth')) + ylim(0,5000)
  
  
  I_res_2 <- mod_rw_intercept$BUGSoutput$sims.matrix[1:1000,(118+N+1):(118+2*N )]
  df_plot <- data.frame(y=c(t(I_res_2)),x=rep(1:N,1000))
  ggplot(df_plot,aes(x=x,y=y,group=1))+ geom_point(alpha=.1) + geom_line(data=data.frame(x=1:(ncol(current_season_sate_flu_matrix)),y=unname(current_season_sate_flu_matrix[2,])),aes(x=x,y=y,col='truth')) 
     
  