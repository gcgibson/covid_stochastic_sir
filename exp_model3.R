

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

    # Season intecepts 
    for (i in 1:8){
      X_season_intercept[i] ~ dnorm(0,1)
    }
  
    
    for (r in 1:53){
      I[1,r] ~ dbeta(100*.02,100*(1-.98))
      S[1,r] ~ dbeta(100*.98,100*(.02))
      for(i in 2:33) {
          Ik1[i,r] <- beta * S[i-1,r] * I[i-1,r]- gamma * I[i-1,r]
          Sk1[i,r] <- -beta*S[i-1,r]*I[i-1,r]
  
          Ik2[i,r] <- beta * (S[i-1,r] +  .5*Sk1[i,r]) * (I[i-1,r] + .5*Ik1[i,r]) - gamma * (I[i-1,r] + .5*Ik1[i,r])
          Sk2[i,r] <- -beta*(S[i-1,r] +  .5*Sk1[i,r])*(I[i-1,r] + .5*Ik1[i,r])
  
          Ik3[i,r] <- beta*(S[i-1,r] +.5*Sk2[i,r])*(I[i-1,r] + .5*Ik2[i,r]) - gamma*(I[i-1,r] + .5*Ik2[i,r])
          Sk3[i,r] <-  -beta*(S[i-1,r] + .5*Sk2[i,r])*(I[i-1,r] + .5*Ik2[i,r])
  
          Ik4[i,r] <- beta*(S[i-1,r] + Sk3[i,r])*(I[i-1,r] + Ik3[i,r])-gamma*(I[i-1,r] + Ik3[i,r])
          Sk4[i,r] <- -beta*(S[i-1,r] + Sk3[i,r])*(I[i-1,r] + Ik3[i,r])
  
  
          mean_I[i,r] <- I[i-1,r] +  .166666* (Ik1[i,r] + 2 * Ik2[i,r] + 2*Ik3[i,r] +Ik4[i,r] )
          mean_S[i,r] <- S[i-1,r] +  .16666* (Sk1[i,r] + 2 * Sk2[i,r] + 2*Sk3[i,r] +Sk4[i,r] )
          S[i,r]  <- mean_S[i,r]
          I[i,r]  <- mean_I[i,r]
      }
      
  
      I_corona[1,r] ~ dbeta(100*.02,100*(1-.98))
      S_corona[1,r] ~ dbeta(100*.98,100*(.02))
      for(i in 2:33) {
        Ik1_corona[i,r] <- beta_corona * S_corona[i-1,r] * I_corona[i-1,r]- gamma_corona * I_corona[i-1,r]
        Sk1_corona[i,r] <- -beta_corona*S_corona[i-1,r]*I_corona[i-1,r]
        
        Ik2_corona[i,r] <- beta_corona * (S_corona[i-1,r] +  .5*Sk1_corona[i,r]) * (I_corona[i-1,r] + .5*Ik1_corona[i,r]) - gamma_corona * (I_corona[i-1,r] + .5*Ik1_corona[i,r])
        Sk2_corona[i,r] <- -beta_corona*(S_corona[i-1,r] +  .5*Sk1_corona[i,r])*(I_corona[i-1,r] + .5*Ik1_corona[i,r])
        
        Ik3_corona[i,r] <- beta_corona*(S_corona[i-1,r] +.5*Sk2_corona[i,r])*(I_corona[i-1,r] + .5*Ik2_corona[i,r]) - gamma_corona*(I_corona[i-1,r] + .5*Ik2_corona[i,r])
        Sk3_corona[i,r] <-  -beta_corona*(S_corona[i-1,r] + .5*Sk2_corona[i,r])*(I_corona[i-1,r] + .5*Ik2_corona[i,r])
        
        Ik4_corona[i,r] <- beta_corona*(S_corona[i-1,r] + Sk3_corona[i,r])*(I_corona[i-1,r] + Ik3_corona[i,r])-gamma_corona*(I_corona[i-1,r] + Ik3_corona[i,r])
        Sk4_corona[i,r] <- -beta_corona*(S_corona[i-1,r] + Sk3_corona[i,r])*(I_corona[i-1,r] + Ik3_corona[i,r])
        
        
        mean_I_corona[i,r] <- I_corona[i-1,r] +  .166666* (Ik1_corona[i,r] + 2 * Ik2_corona[i,r] + 2*Ik3_corona[i,r] +Ik4_corona[i,r] )
        mean_S_corona[i,r] <- S_corona[i-1,r] +  .16666* (Sk1_corona[i,r] + 2 * Sk2_corona[i,r] + 2*Sk3_corona[i,r] +Sk4_corona[i,r] )
        S_corona[i,r]  <- mean_S_corona[i,r]
        I_corona[i,r]  <- mean_I_corona[i,r]
      }

    }
    
    X_interaction[1] ~ dnorm(0,1)

    for (i in 2:N){
        X_interaction[i] ~ dnorm(X_interaction[i-1],8000)
        I_dis[i] <- ifelse(d_idx[i] == 1,I[s_idx[i],r_idx[i]],I_corona[s_idx[i],r_idx[i]])
        y_mean[i] <- I_dis[i]  + X_interaction[i] # + X_season[s_idx[i]] + X_region[r_idx[i]]
        y[i] ~ dnorm(y_mean[i],100000)
    }


}
"

library(cdcForecastUtils)
#epi parameters
flu_beta <- 2
flu_gamma <-  1.4


state_data <- download_and_preprocess_state_flu_data(latest_year = 2020)

# get unique regions
unique_regions <- unique(state_data$region)
unique_regions <- unique_regions[unique_regions != "Florida"]
unique_regions <- unique_regions[unique_regions != "Commonwealth of the Northern Mariana Islands"]
#unique_regions <- c(sample(unique_regions,10),"New York","New York City")

data_frame_for_fit <- data.frame()

for (location_itr in 1:length(unique_regions)){
    location <- unique_regions[location_itr]
    state_data_ny <- state_data[state_data$region == location ,]
    state_data_ny <- state_data_ny[state_data_ny$season != "2019/2020",]
    state_data_ny <- state_data_ny[state_data_ny$week %in% c(40:52,1:20),]
    state_data_ny <- state_data_ny[ state_data_ny$season >= "2012/2013",]
    
    current_season <- state_data[state_data$region == location & state_data$year == "2020" & state_data$week >= 10,]$unweighted_ili
    s_idx <- c(rep(1:33,length.out=length(state_data_ny$unweighted_ili)),1:7)
    y <- c(state_data_ny$unweighted_ili/100,current_season/100,rep(NA,4))
    d_idx <- c(rep(1,length(state_data_ny$unweighted_ili)),rep(2,7))
    int_idx <- c(rep(1:length(unique(state_data_ny$season)),each=33),rep(8,7))
    data_frame_for_fit <- rbind(data_frame_for_fit,data.frame(s_idx=s_idx,y=y,r_idx=location_itr,d_idx=d_idx,int_idx=int_idx))
}


dat <- list(beta=flu_beta,gamma=flu_gamma,N=nrow(data_frame_for_fit),beta_corona=3*.2,gamma_corona=3*1/14,
            y=data_frame_for_fit$y,s_idx=data_frame_for_fit$s_idx,r_idx=data_frame_for_fit$r_idx,d_idx=data_frame_for_fit$d_idx,int_idx=data_frame_for_fit$int_idx)
            

## add to dataframe
jgs <- jags.model(file = textConnection(model), data = dat, n.adapt = 1000)
update(jgs, 1000)
out_jags <- jags.samples(jgs, c('y','y_mean'), 3000, 3)
y_samp <- matrix(out_jags$y_mean,ncol=1000)


## set up submission
library(cdcForecastUtils)

season_start_ew <- "2020-ew10"
season_end_ew <- "2020-ew35"
cdc_report_ew <- get_current_date_from_flu_data(state_data)

get_trajectories_one_location <- function(
  nsim,
  location,
  flu_data,
  target_variable,
  season_start_ew,
  season_end_ew,
  cdc_report_ew,
  targets) {
  # subset to location data
  location_data <- flu_data[flu_data$region == location,]
  
  # prepend observed data
  time_from_start_of_season <- get_time_from_start_of_season(season_start_ew, cdc_report_ew)
  local_ridx <- which(unique_regions==location)
  trajectory_matrix_idx <- tail(which(data_frame_for_fit$r_idx == local_ridx),time_from_start_of_season+4)

  trajectory_matrix <- y_samp[trajectory_matrix_idx,]
  trajectory_matrix[trajectory_matrix<0] <-0 
  return(t(trajectory_matrix)*100)
}

  ##
library(dplyr)
trajectories_by_location <- tibble(
  location = unique_regions
) %>%
  mutate(
    trajectories = purrr::map(
      location,
      get_trajectories_one_location,
      nsim = 1000,
      flu_data = state_data,
      target_variable = "weighted_ili",
      season_start_ew = season_start_ew,
      season_end_ew = season_end_ew,
      cdc_report_ew = cdc_report_ew,
      targets = targets)
  )

distributional_submission_df <- multi_trajectories_to_binned_distributions(
  multi_trajectories = trajectories_by_location,
  targets = "wk ahead",
  h_max = 4,
  bins = c(seq(0, 25, by = .1), 100),
  season_start_ew = season_start_ew,
  season_end_ew = season_end_ew,
  cdc_report_ew = cdc_report_ew)

generate_csv_from_submission_df(distributional_submission_df,"covid-19-ili-forecasting-models/submissions/state-UMassCoE-SIR/2020-ew12-SIR.csv")


### plotting
plot_locations <- c("New York","New Jersey","California","Massachusetts","New Mexico", "New York City")

distributional_submission_df_to_plot <- distributional_submission_df[distributional_submission_df$location %in% plot_locations,]
trajectories_by_location_to_plot <- trajectories_by_location[trajectories_by_location$location %in% plot_locations,]
plot_to_save <- plot_trajectories_and_intervals(
  flu_data = state_data,
  target_variable = "unweighted_ili",
  trajectories_by_location = trajectories_by_location_to_plot,
  submission = distributional_submission_df_to_plot,
  season_start_ew = season_start_ew,
  season_end_ew = season_end_ew,
  cdc_report_ew = cdc_report_ew
)
pdf("states.pdf")
print (plot_to_save)
dev.off()
