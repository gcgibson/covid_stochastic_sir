

library(rjags)
model <- "
model {
    
    # overall  random walk 
    X_season[1] ~ dnorm(0,1)
    for (i in 2:33){
      X_season[i] ~ dnorm(X_season[i-1],10000)
    }

    # region specific random walk 
    X_region[1] ~ dnorm(0,1)
    for (i in 2:54){
          X_region[i] ~ dnorm(X_region[i-1],10000)
    }
  
    

    I[1] ~ dbeta(100*.02,100*(1-.98))
    S[1] ~ dbeta(100*.98,100*(.02))
    for(i in 2:N) {
        Ik1[i] <- beta * S[i-1] * I[i-1]- gamma * I[i-1]
        Sk1[i] <- -beta*S[i-1]*I[i-1]

        Ik2[i] <- beta * (S[i-1] +  .5*Sk1[i]) * (I[i-1] + .5*Ik1[i]) - gamma * (I[i-1] + .5*Ik1[i])
        Sk2[i] <- -beta*(S[i-1] +  .5*Sk1[i])*(I[i-1] + .5*Ik1[i])

        Ik3[i] <- beta*(S[i-1] +.5*Sk2[i])*(I[i-1] + .5*Ik2[i]) - gamma*(I[i-1] + .5*Ik2[i])
        Sk3[i] <-  -beta*(S[i-1] + .5*Sk2[i])*(I[i-1] + .5*Ik2[i])

        Ik4[i] <- beta*(S[i-1] + Sk3[i])*(I[i-1] + Ik3[i])-gamma*(I[i-1] + Ik3[i])
        Sk4[i] <- -beta*(S[i-1] + Sk3[i])*(I[i-1] + Ik3[i])


        mean_I[i] <- I[i-1] +  .166666* (Ik1[i] + 2 * Ik2[i] + 2*Ik3[i] +Ik4[i] )
        mean_S[i] <- S[i-1] +  .16666* (Sk1[i] + 2 * Sk2[i] + 2*Sk3[i] +Sk4[i] )
        S[i]  <- mean_S[i]
        I[i]  <- mean_I[i]
    }


    I_corona[1] ~ dbeta(100*.02,100*(1-.98))
    S_corona[1] ~ dbeta(100*.98,100*(.02))
    for(i in 2:N) {
      Ik1_corona[i] <- beta_corona * S_corona[i-1] * I_corona[i-1]- gamma_corona * I_corona[i-1]
      Sk1_corona[i] <- -beta_corona*S_corona[i-1]*I_corona[i-1]
      
      Ik2_corona[i] <- beta_corona * (S_corona[i-1] +  .5*Sk1_corona[i]) * (I_corona[i-1] + .5*Ik1_corona[i]) - gamma_corona * (I_corona[i-1] + .5*Ik1_corona[i])
      Sk2_corona[i] <- -beta_corona*(S_corona[i-1] +  .5*Sk1_corona[i])*(I_corona[i-1] + .5*Ik1_corona[i])
      
      Ik3_corona[i] <- beta_corona*(S_corona[i-1] +.5*Sk2_corona[i])*(I_corona[i-1] + .5*Ik2_corona[i]) - gamma_corona*(I_corona[i-1] + .5*Ik2_corona[i])
      Sk3_corona[i] <-  -beta_corona*(S_corona[i-1] + .5*Sk2_corona[i])*(I_corona[i-1] + .5*Ik2_corona[i])
      
      Ik4_corona[i] <- beta_corona*(S_corona[i-1] + Sk3_corona[i])*(I_corona[i-1] + Ik3_corona[i])-gamma_corona*(I_corona[i-1] + Ik3_corona[i])
      Sk4_corona[i] <- -beta_corona*(S_corona[i-1] + Sk3_corona[i])*(I_corona[i-1] + Ik3_corona[i])
      
      
      mean_I_corona[i] <- I_corona[i-1] +  .166666* (Ik1_corona[i] + 2 * Ik2_corona[i] + 2*Ik3_corona[i] +Ik4_corona[i] )
      mean_S_corona[i] <- S_corona[i-1] +  .16666* (Sk1_corona[i] + 2 * Sk2_corona[i] + 2*Sk3_corona[i] +Sk4_corona[i] )
      S_corona[i]  <- mean_S_corona[i]
      I_corona[i]  <- mean_I_corona[i]
    }


    
    X_interaction[1] ~ dnorm(0,1)

    for (i in 2:N){
        X_interaction[i] ~ dnorm(X_interaction[i-1],10000000)
        I_dis[i] <- ifelse(d_idx[i] == 1,I[s_idx[i]],I_corona[s_idx[i]])
        y_mean[i] <- X_interaction[i] + X_season[s_idx[i]] +I_dis[i] + X_region[r_idx[i]]
        y[i] ~ dnorm(y_mean[i],10000000)
    }


}
"

library(cdcForecastUtils)
#epi parameters
flu_beta <- 2
flu_gamma <-  1.4




# get unique regions
unique_regions <- unique(state_data$region)
unique_regions <- unique_regions[unique_regions != "Florida"]


#state_data <- download_and_preprocess_state_flu_data(latest_year = 2020)

data_frame_for_fit <- data.frame()

for (location_itr in 1:length(unique_regions)){
    location <- unique_regions[location_itr]
    state_data_ny <- state_data[state_data$region == location ,]
    state_data_ny <- state_data_ny[state_data_ny$season != "2019/2020",]
    state_data_ny <- state_data_ny[state_data_ny$week %in% c(40:52,1:20),]
    state_data_ny <- state_data_ny[ state_data_ny$season >= "2012/2013",]
    
    current_season <- state_data[state_data$region == location & state_data$year == "2020" & state_data$week >= 10,]$unweighted_ili
    s_idx <- c(rep(1:33,length.out=length(state_data_ny$unweighted_ili)),1:6)
    y <- c(state_data_ny$unweighted_ili,current_season,rep(NA,4))
    d_idx <- c(rep(1,length(state_data_ny$unweighted_ili)),rep(2,6))
    data_frame_for_fit <- rbind(data_frame_for_fit,data.frame(s_idx=s_idx,y=y,r_idx=location_itr,d_idx=d_idx))
}


dat <- list(beta=flu_beta,gamma=flu_gamma,N=nrow(data_frame_for_fit),beta_corona=.2,gamma_corona=1/14,
            y=data_frame_for_fit$y,s_idx=data_frame_for_fit$s_idx,r_idx=data_frame_for_fit$r_idx,d_idx=data_frame_for_fit$d_idx)
            

## add to dataframe
jgs <- jags.model(file = textConnection(model), data = dat, n.adapt = 1000)
update(jgs, 1000)
out_jags <- jags.samples(jgs, c('y','y_mean'), 3000, 3)
data_frame_for_fit$y_mean <- out_jags$y_mean



