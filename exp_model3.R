sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 1.4247, gamma = 0.14286)
## Time frame
times      <- seq(0, 70, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
head(out, 10)

plot(cumsum(diff(out$I)))


library(rjags)
model <- "
model {
    
    # overall specific random walk 
    X[1] ~ dnorm(0,1)
    for (i in 2:33){
      X[i] ~ dnorm(X[i-1],1000)
    }
  
    #gamma ~ dgamma(1,14)
    #beta ~ dbeta(1000*.2,1000 *(1-.2))
    

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


    
    Z[1] ~ dnorm(0,1)

    for (i in 2:N){
        Z[i] ~ dnorm(Z[i-1],1000)
        I_dis[i] <- ifelse(disease[i] == 1,I[s_idx[i]],I_corona[s_idx[i]])
        y_mean[i] <- Z[i] + X[s_idx[i]] +I_dis[i]
        y[i] ~ dnorm(y_mean[i],1000)
    }


}
"

library(cdcForecastUtils)

#state_data <- download_and_preprocess_state_flu_data(latest_year = 2020)
location <- "New York"
flu_beta <- 2
flu_gamma <-  1.4
state_data_ny <- state_data[state_data$region == location ,]
state_data_ny <- state_data_ny[state_data_ny$season != "2019/2020",]
state_data_ny <- state_data_ny[state_data_ny$week %in% c(40:52,1:20),]
state_data_ny <- state_data_ny[ state_data_ny$season >= "2012/2013",]

ggplot(state_data_ny,aes(x=rep(1:33,length.out=length(unweighted_ili)),y=unweighted_ili,col=season)) + geom_line() + facet_grid(~season)
  
current_season <- state_data[state_data$region == location & state_data$year == "2020" & state_data$week >= 10,]$unweighted_ili-3


x <- c(rep(1:33,length.out=length(state_data_ny$unweighted_ili)),1:6)
y <- c(state_data_ny$unweighted_ili,current_season,rep(NA,4))
dat <- list(beta=flu_beta,gamma=flu_gamma,N=length(x),beta_corona=.2,gamma_corona=1/14,
            y=y,s_idx=x,disease=c(rep(1,length(state_data_ny$unweighted_ili)),rep(2,length(current_season)+4)))
jgs <- jags.model(file = textConnection(model), data = dat, n.adapt = 1000)
update(jgs, 1000)
out_jags <- jags.samples(jgs, c('y','y_mean'), 3000, 3)

adjustment <- c(rep(0,nrow(state_data_ny)),rep(2,length(current_season)+4))
plot(rowMeans(out_jags$y_mean) +adjustment ,type='l')
lines(y+adjustment,col='red')
