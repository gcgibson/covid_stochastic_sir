

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
      
      beta_corona_t[1,r] <- beta_corona
      I_corona[1,r] ~ dbeta(100*.02,100*(1-.98))
      S_corona[1,r] ~ dbeta(100*.98,100*(.02))
      for(i in 2:33) {
        beta_corona_t[i,r] ~ dnorm(.75*beta_corona_t[i-1,r] ,1000)
        Ik1_corona[i,r] <- beta_corona_t[i,r]  * S_corona[i-1,r] * I_corona[i-1,r]- gamma_corona * I_corona[i-1,r]
        Sk1_corona[i,r] <- -beta_corona_t[i,r] *S_corona[i-1,r]*I_corona[i-1,r]
        
        Ik2_corona[i,r] <- beta_corona_t[i,r]  * (S_corona[i-1,r] +  .5*Sk1_corona[i,r]) * (I_corona[i-1,r] + .5*Ik1_corona[i,r]) - gamma_corona * (I_corona[i-1,r] + .5*Ik1_corona[i,r])
        Sk2_corona[i,r] <- -beta_corona_t[i,r] *(S_corona[i-1,r] +  .5*Sk1_corona[i,r])*(I_corona[i-1,r] + .5*Ik1_corona[i,r])
        
        Ik3_corona[i,r] <- beta_corona_t[i,r] *(S_corona[i-1,r] +.5*Sk2_corona[i,r])*(I_corona[i-1,r] + .5*Ik2_corona[i,r]) - gamma_corona*(I_corona[i-1,r] + .5*Ik2_corona[i,r])
        Sk3_corona[i,r] <-  -beta_corona_t[i,r] *(S_corona[i-1,r] + .5*Sk2_corona[i,r])*(I_corona[i-1,r] + .5*Ik2_corona[i,r])
        
        Ik4_corona[i,r] <- beta_corona_t[i,r] *(S_corona[i-1,r] + Sk3_corona[i,r])*(I_corona[i-1,r] + Ik3_corona[i,r])-gamma_corona*(I_corona[i-1,r] + Ik3_corona[i,r])
        Sk4_corona[i,r] <- -beta_corona_t[i,r] *(S_corona[i-1,r] + Sk3_corona[i,r])*(I_corona[i-1,r] + Ik3_corona[i,r])
        
        
        mean_I_corona[i,r] <- I_corona[i-1,r] +  .166666* (Ik1_corona[i,r] + 2 * Ik2_corona[i,r] + 2*Ik3_corona[i,r] +Ik4_corona[i,r] )
        mean_S_corona[i,r] <- S_corona[i-1,r] +  .16666* (Sk1_corona[i,r] + 2 * Sk2_corona[i,r] + 2*Sk3_corona[i,r] +Sk4_corona[i,r] )
        S_corona[i,r]  <- mean_S_corona[i,r]
        I_corona[i,r]  <- mean_I_corona[i,r]
      }

    }
    
    #Covid percent case likelihood
  
    for (covid_idx in 1:covid_length){
      state_percent[covid_idx] ~ dnorm(I_corona[state_week[covid_idx],state_idx[covid_idx]],100)
    }



    #ILI likelihood
    X_interaction[1] ~ dnorm(0,1)

    for (i in 2:N){
        X_interaction[i] ~ dnorm(X_interaction[i-1],8000)
        I_dis[i] <- ifelse(d_idx[i] == 1,I[s_idx[i],r_idx[i]],I_corona[s_idx[i],r_idx[i]])
        y_mean[i] <- I_dis[i]  + X_interaction[i]  + X_season[s_idx[i]] + X_region[r_idx[i]]
        y[i] ~ dnorm(y_mean[i],10000)
    }


}
"

library(cdcForecastUtils)
library(openintro)
#epi parameters
flu_beta <- 2
flu_gamma <-  1.4
state_testing_data <- read.csv("/Users/gcgibson/Downloads/states-daily.csv")

#state_data <- download_and_preprocess_state_flu_data(latest_year = 2020)

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

state_testing_data$date_formatted <- as.Date(state_testing_data$dateChecked)

state_testing_data$date_formatted <- unlist(lapply(state_testing_data$date,function(x){
  year <- substr(x,1,4)
  month <- substr(x,5,6)
  day <- substr(x,7,8)
  return (paste0(year,"/",month,"/",day))
}))
state_testing_data$date_formatted <- as.Date(state_testing_data$date_formatted)
state_testing_data$week <- lubridate::week(state_testing_data$date_formatted)

state_testing_data_by_week <- state_testing_data %>% group_by(week,state) %>% summarise(week_positive = sum(positive,na.rm=T))

popoulation<- state.x77[,1]
population_df <- data.frame(names = names(popoulation),pop=popoulation*100)
unique_regions_df <- data.frame(reg =unique_regions, reg_idx = 1:length(unique_regions))
  


state_testing_data_by_week$state_name <- abbr2state(state_testing_data_by_week$state)

state_testing_data_by_week$pop <-population_df[state_testing_data_by_week$state_name,]$pop
state_testing_data_by_week$percent <- state_testing_data_by_week$week_positive/state_testing_data_by_week$pop

state_testing_data_by_week$state_idx <- unique_regions_df[as.factor(state_testing_data_by_week$state_name)  ,]$reg_idx
state_testing_data_by_week <- state_testing_data_by_week[complete.cases(state_testing_data_by_week),]

dat <- list(beta=flu_beta,gamma=flu_gamma,N=nrow(data_frame_for_fit),beta_corona=3*.2,gamma_corona=3*1/14,
            y=data_frame_for_fit$y + 1e-10,s_idx=data_frame_for_fit$s_idx,r_idx=data_frame_for_fit$r_idx,d_idx=data_frame_for_fit$d_idx,int_idx=data_frame_for_fit$int_idx,
            state_percent=state_testing_data_by_week$percent,state_idx=state_testing_data_by_week$state_idx,
            state_week=state_testing_data_by_week$week-9,covid_length = length(state_testing_data_by_week$week))
            

## add to dataframe
jgs <- jags.model(file = textConnection(model), data = dat, n.adapt = 1000)
update(jgs, 1000)
out_jags <- jags.samples(jgs, c('y','y_mean'), 3000, 3)
y_samp <- matrix(out_jags$y,ncol=1000)


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

state_point <- generate_point_forecasts(distributional_submission_df)
state_distributional_submission_df <- sanitize_entry(rbind(distributional_submission_df,state_point))
verify_entry(state_distributional_submission_df,challenge = "state_ili")
generate_csv_from_submission_df(state_distributional_submission_df,"covid-19-ili-forecasting-models/submissions/state-UMassCoE-SIR/2020-ew12-SIR")


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


regions <- unique(distributional_submission_df$location)
census_2010 <- read.csv("census_2010.csv")
census_2010 <- census_2010 %>% mutate(weight = year_2010 / sum(year_2010))
census_2010$location <- census_2010$state
census_2010$location <- unlist(lapply( census_2010$state,function(x){
  if (grepl("_",x)){
    print (x)
    return (paste0(strsplit(as.character(x),"_")[[1]],collapse=" "))
  } else{
    return (as.character(x))
  }
}))
census_2010$location <- as.factor(census_2010$location)

distributional_submission_df$location <- as.factor(distributional_submission_df$location)

levels(distributional_submission_df$location) <- levels(census_2010$location)

regional_weight_df <- data.frame()

for (r in c(as.vector(census_2010$hhs_region),"National")){
  if (r != "National"){
    local_census <- census_2010[census_2010$hhs_region == r, ]
    local_census$weight <- local_census$weight#/sum(local_census$weight)
     
    local_census$region <- r
    regional_weight_df <-rbind(regional_weight_df,data.frame(local_census))
  } else{
    census_2010$region <- r
    
    regional_weight_df <-rbind(regional_weight_df,data.frame(census_2010))
    
  }
}

merged_df <- merge(data.frame(distributional_submission_df),regional_weight_df,all=TRUE)


regional_df <- merged_df %>% group_by(region,bin,target,type) %>% summarize(value=sum(value*weight/(sum(weight))))

regional_df$location <- recode(regional_df$region,hhs1="HHS Region 1",hhs2="HHS Region 2",
                                hhs3="HHS Region 3",hhs4="HHS Region 4",hhs5="HHS Region 5",
                               hhs6 = "HHS Region 6",hhs7="HHS Region 7",hhs8="HHS Region 8",
                               hhs9="HHS Region 9",hhs10="HHS Region 10",National="US National")







regional_df$bin <- as.character(regional_df$bin)
regional_df$type
regional_df <- regional_df[complete.cases(regional_df),]


regional_df_sanitized <- sanitize_entry(regional_df)
point_forecasts <- generate_point_forecasts(regional_df_sanitized)
regional_df_sanitized <-rbind(point_forecasts,regional_df_sanitized[,2:ncol(regional_df_sanitized)])
sum(regional_df_sanitized$type == "point")

verify_entry(regional_df_sanitized,challenge = "ilinet")
generate_csv_from_submission_df(regional_df_sanitized,"covid-19-ili-forecasting-models/submissions/region-UMassCoE-SIR/2020-ew12-SIR")

trajectories_by_location_nat <- tibble(
  location = unique(regional_df_sanitized$location) 
) %>%
  mutate(
    trajectories = purrr::map(
      location,
      function(l) {
        matrix(0, nrow = 1, ncol = 26)
      }
    )
  )

regional_data <- download_and_preprocess_flu_data()
regional_data$region <- recode(regional_data$region, "Region 1" = "HHS Region 1","Region 2" = "HHS Region 2",
                               "Region 3" = "HHS Region 3","Region 4" = "HHS Region 4","Region 5" = "HHS Region 5",
                               "Region 6" = "HHS Region 6","Region 7" = "HHS Region 7","Region 8" = "HHS Region 8",
                               "Region 9" = "HHS Region 9","Region 10" = "HHS Region 10","National"="US National")

regional_df_sanitized_tmp <- regional_df_sanitized[regional_df_sanitized$location %in% c("HHS Region 1","US National"),]

trajectories_by_location_nat <- tibble(
  location = unique(regional_df_sanitized_tmp$location) 
) %>%
  mutate(
    trajectories = purrr::map(
      location,
      function(l) {
        matrix(0, nrow = 1, ncol = 26)
      }
    )
  )

plot_trajectories_and_intervals(
  flu_data = regional_data,
  target_variable = "weighted_ili",
  trajectories_by_location = trajectories_by_location_nat,
  submission = regional_df_sanitized_tmp,
  season_start_ew = season_start_ew,
  season_end_ew = season_end_ew,
  cdc_report_ew = cdc_report_ew
)
