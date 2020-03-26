nyc_data <- read.delim("/Users/gcgibson/covid/A_EXPORT_QL_crosstab.csv", fileEncoding="UCS-2LE")
nyc_data_city_wide <- nyc_data[nyc_data$Dim1Name == "Citywide",]
nyc_data_city_wide_all_age <- nyc_data_city_wide[nyc_data_city_wide$Dim2Value=="All age groups",]
nyc_data_city_wide_all_age$X <- as.numeric(gsub(",", "", nyc_data_city_wide_all_age$X))

nyc_data_ili <- read.delim("/Users/gcgibson/covid/A_EXPORT_QL_crosstab_ili.csv", fileEncoding="UCS-2LE")
nyc_data_ili_city_wide <- nyc_data_ili[nyc_data_ili$Dim1Name == "Citywide",]
nyc_data_ili_city_wide_all_age <- nyc_data_ili_city_wide[nyc_data_ili_city_wide$Dim2Value=="All age groups",]
nyc_data_ili_city_wide_all_age$X <- as.numeric(gsub(",", "", nyc_data_ili_city_wide_all_age$X))
nyc_data_ili_city_wide_all_age$date <- as.Date(as.character(nyc_data_ili_city_wide_all_age$Date),format ="%m/%d/%y")


covid_ts <- read.csv("COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv")
covid_ts_us <- covid_ts[covid_ts$Province.State == "New York",]
covid_ts_us <- (colSums(covid_ts_us[,5:ncol(covid_ts_us)]))
covid_ts_us <- head(covid_ts_us,length(covid_ts_us)-2)
hist_t <- 100
ili_vs_covid_cases <- data.frame(t=tail(nyc_data_ili_city_wide_all_age$date,hist_t),ili_cases=tail(nyc_data_city_wide_all_age$X,hist_t),
                                 covid=c(rep(NA,hist_t-length(covid_ts_us)),tail(covid_ts_us,hist_t)),ili_percent=tail(nyc_data_ili_city_wide_all_age$X,hist_t))
write.csv(ili_vs_covid_cases$ili_percent,"/Users/gcgibson/covid/ili.csv")
library(ggplot2)
library(scales)

ggplot(ili_vs_covid_cases,aes(x=t,y=ili_cases,col='ili cases',group=1))+ geom_line() + geom_line(aes(x=t,y=covid,col='covid',group=1))+ theme_bw() + scale_x_date(labels = date_format("%d-%m"))




cat(
  '
  functions {
  
  // This largely follows the deSolve package, but also includes the x_r and x_i variables.
  // These variables are used in the background.
  
  real[] SI(real t,
  real[] y,
  real[] params,
  real[] x_r,
  int[] x_i) {
  
  real dydt[3];
  
  dydt[1] = - params[1] * y[1] * y[2];
  dydt[2] = params[1] * y[1] * y[2] - params[2] * y[2];
  dydt[3] = params[2] * y[2];
  
  return dydt;
  }
  
  }
  
  data {
    int<lower = 1> n_obs; // Number of days sampled
    int<lower = 1> n_params; // Number of model parameters
    int<lower = 1> n_difeq; // Number of differential equations in the system
    int<lower = 1> n_sample; // Number of hosts sampled at each time point.
    int<lower = 1> n_fake; // This is to generate "predicted"/"unsampled" data
    
    real y[n_obs]; // The binomially distributed data
    real t0; // Initial time point (zero)
    real ts[n_obs]; // Time points that were sampled
    
    real fake_ts[n_fake]; // Time points for "predicted"/"unsampled" data
  }
  
  transformed data {
  real x_r[0];
  int x_i[0];
  }
  
  parameters {
  real<lower = 0> params[n_params]; // Model parameters
  real<lower = 0, upper = 1> S0; // Initial fraction of hosts susceptible
  }
  
  transformed parameters{
    real y_hat[n_obs, n_difeq]; // Output from the ODE solver
    real y0[n_difeq]; // Initial conditions for both S and I
    
    y0[1] = S0;
    y0[2] = 1 - S0;
    y0[3] = 0;
    
    y_hat = integrate_ode_bdf(SI, y0, t0, ts, params, x_r, x_i);
  
  }
  
  model {
  params ~ normal(0, 2); //constrained to be positive
  S0 ~ normal(0.5, 0.5); //constrained to be 0-1.
  
  y ~ normal(y_hat[, 2],.00000001); //y_hat[,2] are the fractions infected from the ODE solver
  
  }
  
  generated quantities {
  // Generate predicted data over the whole time series:
  real fake_I[n_fake, n_difeq];
  
  fake_I = integrate_ode_bdf(SI, y0, t0, fake_ts, params, x_r, x_i);
  
  }
  
  ', 
  file = "SI_fit.stan", sep="", fill=T)

# FITTING

# For stan model we need the following variables:
params <- list(beta = beta,
               gamma = gamma)
I0 = 0.02    # initial fraction infected
S0 = 1 - I0 # initial fraction susceptible
R0 = 0
inits <- c(S0, I0, R0)

# Assign transmission and pathogen-induced death rates:
beta = 0.60
gamma = 0.10

# We will use the package deSolve to integrate, which requires certain data structures.
# Store parameters and initial values
# Parameters must be stored in a named list.
params <- list(beta = beta,
               gamma = gamma)

# Initial conditions are
# Initial conditions are stored in a vector
t_min = 0
t_max = 50
times = t_min:t_max

# We must create a function for the system of ODEs.
# See the 'ode' function documentation for further insights.
SIR <- function(t, y, params) {
  with(as.list(c(params, y)), {
    
    dS = - beta * y[1] * y[2]
    
    dI = beta * y[1] * y[2] - gamma * y[2]
    
    dR = gamma * y[2]
    
    res <- c(dS,dI,dR)
    list(res)
  })
}

# Run the integration:
out <- ode(inits, times, SIR, params, method="ode45")

sample_days = 20 # number of days sampled throughout the epidemic
sample_n = 25 # number of host individuals sampled per day

# Choose which days the samples were taken. 
# Ideally this would be daily, but we all know that is difficult.
sample_time = 1:sample_days
sample_propinf = out[, 3]

sample_y = rbinom(sample_days, sample_n, sample_propinf)
sample_prop = sample_y / sample_n

stan_d = list(n_obs = sample_days,
              n_params = length(params),
              n_difeq = length(inits),
              n_sample = sample_n,
              n_fake = length(1:t_max),
              y = sample_y/1000,
              t0 = 0,
              ts = sample_time,
              fake_ts = c(1:t_max))


# Which parameters to monitor in the model:
params_monitor = c("y_hat", "y0", "params", "fake_I")



