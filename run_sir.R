###simulate data

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
times      <- seq(0, 29, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)

plot(out[,3],type='l')

library(rjags)
library(R2jags)
model.loc = ("rw_intercept.txt")
jagsscript =  cat("
                    model {
                        S[1] <- pop - infected
                        I[1] <- infected
                        R[1] <- 0
                        beta ~ dgamma(1,1)
                        gamma ~ dgamma(1,1)
                                      
                                  
                                   
                     for(i in 2:N) {
                              Sk1[i] <- -beta*S[i-1]*I[i-1]
                              Ik1[i] <- beta * S[i-1] * I[i-1] - gamma * I[i-1]

                              Sk2[i] <- -beta*(S[i-1] +  .5* Sk1[i])*(I[i-1] + .5*Ik1[i])
                              Ik2[i] <- beta * (S[i-1] +  .5* Sk1[i]) * (I[i-1] + .5*Ik1[i]) - gamma * (I[i-1] + .5*Ik1[i])

                              

                              Sk3[i] <-  -beta*(S[i-1] + .5* Sk2[i])*(I[i-1] + .5*Ik2[i])
                              Ik3[i] <- beta*(S[i-1] +.5*Sk2[i])*(I[i-1] + .5*Ik2[i]) - gamma*(I[i-1] + .5*Ik2[i])
                  
                              Sk4[i] <- -beta*(S[i-1] + Sk3[i])*(I[i-1] + Ik3[i])
                              Ik4[i] <- beta*(S[i-1] + Sk3[i])*(I[i-1] + Ik3[i])-gamma*(I[i-1] + Ik3[i])

                              S[i] <- S[i-1] +  .16666* (Sk1[i] + 2 * Sk2[i] + 2*Sk3[i] +Sk4[i] )
                              I[i]  <- I[i-1] +  .166666* (Ik1[i] + 2 * Ik2[i] + 2*Ik3[i] +Ik4[i] )
  
                              R[i] <- 1-S[i]-I[i]
                      }
                      
                      for(i in 1:N){
                            I_rand[i] ~ dnorm(I[i],1000)
                            Y[i] ~dnorm(I_rand[i],1000000)
                      }
  

}
                                 ", 
                 file = model.loc)

I_sim <- out[,3]
jags.data = list(Y = I_sim,N=30,infected=.00000001,pop=1)
jags.params = c("Y","I","beta","gamma")
mod_rw_intercept = jags(jags.data, parameters.to.save = jags.params, 
                              model.file = model.loc, n.chains = 3, n.burnin = 5000, n.thin = 1, 
                              n.iter = 10000, DIC = TRUE)


library(ggplot2)
ggplot()
df_plot <- data.frame(y=c(t(mod_rw_intercept$BUGSoutput$sims.matrix[1:1000,paste0("I[",2:30,"]")])),x=rep(2:30,1000))
ggplot(df_plot,aes(x=x,y=y))+ geom_line() + ylim(0,1)#+ geom_line(data=data.frame(x=1:30,y=I_sim),aes(x=x,y=y,col='truth')) + ylim(0,1)

mean((mod_rw_intercept$BUGSoutput$sims.matrix[,"beta"]))
mean((mod_rw_intercept$BUGSoutput$sims.matrix[,"gamma"]))

