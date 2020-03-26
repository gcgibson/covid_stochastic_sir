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


library(rjags)
model <- "model {
  X[1] ~ dnorm(3.8,100)
  for (i in 2:N){
    X[i] ~ dnorm(X[i-1],1000)
  }
  for (i in 1:N) {
    S[i] <- 1-1/(1+exp(-.2*i+4))
    #mu0[i] = ifelse(S[i] < .2/3,  mu1*exp(-i*beta1),  mu2*log(i*beta2))
    #mu0[i] = ifelse(i < cp,  mu1*exp(-i*beta1),   mu2*exp(-i*beta2))
  
    mu0[i] =  mu1*exp(-i*beta1)
  
    x[i] ~ dnorm(mu0[i], 1000)
    ili[i] ~ dnorm(X[i] +mu0[i],100000)
  }  
mu1 ~ dnorm(0, 1)
mu2 ~ dnorm(0, 1)
beta1 ~ dnorm(0, 1)
beta2 ~ dnorm(1, 1)
cp ~ dunif(0,N)
}"

library(cdcForecastUtils)

state_data <- download_and_preprocess_state_flu_data(latest_year = 2020)
location <- "New York"
state_data_ny <- state_data[state_data$year == "2020" & state_data$season_week >=10 & state_data$region == location,]$unweighted_ili 

dat <- list('x'=c(covid_ts_us/10000,rep(NA,4)),'ili'=c(tail(state_data_ny,10),rep(NA,4)), 'N'=length(covid_ts_us)+4)
jgs <- jags.model(file = textConnection(model), data = dat, n.adapt = 1000)
update(jgs, 1000)
out_jags <- jags.samples(jgs, c('ili','x','mu1','mu2','mu0','S'), 3000, 3)

col_nums <- nrow(out_jags$ili)
library(ggplot2)
out_jags_df <- data.frame(x=rep(1:col_nums,1000),y=c(t(matrix(out_jags$ili,nrow=col_nums,byrow = T))))
ggplot(out_jags_df,aes(x=x,y=y)) + geom_line(alpha=.1) + geom_point(data=data.frame(x=1:10,y=tail(state_data_ny,10)),aes(x=x,y=y,col="observed"))  +
  geom_line(data=data.frame(x=1:col_nums,y=rowMeans(out_jags$ili)),aes(x=x,y=y,col='posterior_mean')) + theme_bw()


