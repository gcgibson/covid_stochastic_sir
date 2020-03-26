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
    for (i in 1:N) {
        S[i] <- 1-1/(1+exp(-.2*i+4))
        mu0[i] = ifelse(S[i] < .2/3,  mu1*exp(-i*beta1),  mu2*log(i*beta2))

        x[i] ~ dnorm(mu0[i], 1000)
        
    }  
    mu1 ~ dnorm(0, 1)
    mu2 ~ dnorm(0, 1)
    beta1 ~ dnorm(0, 1)
    beta2 ~ dgamma(1, 1)
    cp ~ dunif(0, N)   
}"

covid_ts <- read.csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
covid_ts_us <- covid_ts[covid_ts$Province.State == "New York",]
covid_ts_us <- (colSums(covid_ts_us[,5:ncol(covid_ts_us)]))
covid_ts_us <- tail(covid_ts_us,10)


dat <- list('x'=c(covid_ts_us/10000,rep(NA,30)), 'N'=length(covid_ts_us)+30)
jgs <- jags.model(file = textConnection(model), data = dat, n.adapt = 1000)
update(jgs, 1000)
out_jags <- jags.samples(jgs, c('cp','mu1','mu2','mu0','S'), 3000, 3)
mean(out_jags$cp)

library(ggplot2)
out_jags_df <- data.frame(x=rep(1:40,1000),y=c(t(matrix(out_jags$mu0,nrow=40,byrow = T))))
ggplot(out_jags_df,aes(x=x,y=y)) + geom_line(alpha=.1) + geom_line(data=data.frame(x=1:10,y=covid_ts_us),aes(x=x,y=y/10000,col="observed"))


