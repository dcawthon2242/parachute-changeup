#===============================================================================
# Packages
#==============================================================================
install.packages("tidyverse")
install.packages("jsonlite")

library(tidyverse)
library(jsonlite)
#===============================================================================
# Convert to Lists (JSON method)
#==============================================================================
# Convert to lists using fromJSon and initialize pitch_trajectory_distance and weighted_trajectory_distance columns
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(trajectory = lapply(trajectory, fromJSON)) %>%
  group_by(game_pk, at_bat_number) %>%
  mutate(pitch_trajectory_distance = NA)

#===============================================================================
# Convert to Lists (String split method)
#==============================================================================
convert_trajectory_to_lists <- function(trajectory_string) {
  # Split the string into separate components based on pattern "), c("
  components <- strsplit(trajectory_string, "\\), c\\(")[[1]]
  
  # Remove any leading "c(" and trailing ")" from each component
  components <- gsub("^c\\(|\\)$", "", components)
  
  # Convert each component to a numeric vector
  trajectory_lists <- lapply(components, function(comp) {
    as.numeric(unlist(strsplit(comp, ", ")))
  })
  
  return(trajectory_lists)
}

# Initialize a list to store results
trajectory_lists_all <- vector("list", length(pbp_2022_2024$trajectory))

# Loop through each trajectory string, convert it, and store the result
for (i in seq_along(pbp_2022_2024$trajectory)) {
  trajectory_lists_all[[i]] <- convert_trajectory_to_lists(pbp_2022_2024$trajectory[i])
  
  # Print progress every 1000 rows
  if (i %% 1000 == 0) {
    print(paste("Processed", i, "rows"))
  }
}

# Add the results back to the dataset
pbp_2022_2024$trajectory_lists <- trajectory_lists_all


#===============================================================================
# Calculate Pitch Trajectory Distance
#===============================================================================
# Load libraries
library(Rcpp)
library(dplyr)
library(jsonlite)

# Rcpp function for fast distance calculation
Rcpp::cppFunction('
  double calculateDistance(NumericVector x1, NumericVector y1, NumericVector z1,
                           NumericVector x2, NumericVector y2, NumericVector z2) {
    int n = x1.size();
    double distance = 0.0;
    
    for (int i = 0; i < n; i++) {
      double dx = x1[i] - x2[i];
      double dy = y1[i] - y2[i];
      double dz = z1[i] - z2[i];
      distance += sqrt(dx * dx + dy * dy + dz * dz);
    }
    return distance;
  }
')

# Wrapper function
add_pitch_distance <- function(df, trajectory_column = "trajectory_lists") {
  
  message("Parsing trajectory JSON...")
  df$trajectory_list <- lapply(df[[trajectory_column]], fromJSON)
  
  message("Starting distance calculations...")
  df <- df %>%
    arrange(game_pk, at_bat_number, pitch_number) %>%
    group_by(game_pk, at_bat_number) %>%
    mutate(pitch_trajectory_distance = NA_real_)
  
  n_rows <- nrow(df)
  
  for (i in 2:n_rows) {
    if (df$game_pk[i] == df$game_pk[i-1] &&
        df$at_bat_number[i] == df$at_bat_number[i-1]) {
      
      traj1 <- df$trajectory_list[[i]]
      traj2 <- df$trajectory_list[[i-1]]
      
      # Extract and unlist
      x1 <- unlist(traj1[[1]])
      y1 <- unlist(traj1[[2]])
      z1 <- unlist(traj1[[3]])
      
      x2 <- unlist(traj2[[1]])
      y2 <- unlist(traj2[[2]])
      z2 <- unlist(traj2[[3]])
      
      len <- min(length(x1), length(x2))  # Ensure same length
      
      dist <- calculateDistance(x1[1:len], y1[1:len], z1[1:len],
                                x2[1:len], y2[1:len], z2[1:len])
      
      df$pitch_trajectory_distance[i] <- dist
    }
    
    # Print progress every 1000 rows
    if (i %% 1000 == 0) {
      message(paste("Processed", i, "rows of", n_rows))
    }
  }
  
  message("All rows processed!")
  return(df)
}

pbp_2022_2024 <- add_pitch_distance(pbp_2022_2024)

write.csv(pbp_2022_2024_with_distance, "pbp_2022_2024_with_distance.csv", row.names = FALSE)


#===============================================================================
# Trajectory Weights Model
#==============================================================================

# Create 100,000 row random subset for use in the model 
set.seed(22)
pbp_2022_2024_subset <- pbp_2022_2024 %>%
  sample_n(100000)

# Function to find optimal weights using correlation to swing decisions
calculate_r_squared <- function(data, x_weight, y_weight, z_weight) {
  
  data_grouped <- data %>% group_by(game_pk, at_bat_number)
  
  data_grouped$weighted_pitch_trajectory_distance <- NA
  
  for (i in 1:nrow(data_grouped)) {
    at_bat_data <- data_grouped[i, ]
    pitch_data <- at_bat_data$trajectory
    
    pitch_df <- data.frame(t = seq(0, at_bat_data$time_of_flight, length.out = 50),
                           x = pitch_data[[1]]$x,
                           y = pitch_data[[1]]$y,
                           z = pitch_data[[1]]$z)
    
    if (i > 1) {
      prev_pitch_data <- data_grouped[i - 1, ]$trajectory
      prev_pitch_df <- data.frame(t = seq(0, data_grouped[i - 1, ]$time_of_flight, length.out = 50),
                                  x = prev_pitch_data[[1]]$x,
                                  y = prev_pitch_data[[1]]$y,
                                  z = prev_pitch_data[[1]]$z)
      
      t1 <- seq(0, max(c(pitch_df$t, prev_pitch_df$t)), length.out = 200)
      x1 <- rep(0, 200)
      x2 <- rep(0, 200)
      y1 <- rep(0, 200)
      y2 <- rep(0, 200)
      z1 <- rep(0, 200)
      z2 <- rep(0, 200)
      h1 <- 0.01
      h2 <- 0.01
      
      for (j in 1:200) {
        x1[j] <- sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2) * pitch_df$x) /
          sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2))
        y1[j] <- sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2) * pitch_df$y) /
          sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2))
        z1[j] <- sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2) * pitch_df$z) /
          sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2))
        x2[j] <- sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2) * prev_pitch_df$x) /
          sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2))
        y2[j] <- sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2) * prev_pitch_df$y) /
          sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2))
        z2[j] <- sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2) * prev_pitch_df$z) /
          sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2))
      }
      
      distance_df <- data.frame(t1, d1 = rep(0, 200), x1, y1, z1, x2, y2, z2)
      distance_df$d1 <- sqrt((x_weight * (distance_df$x1 - distance_df$x2))^2 + (y_weight * (distance_df$y1 - distance_df$y2))^2 + (z_weight * (distance_df$z1 - distance_df$z2))^2) # Weighted Distance
      delta <- distance_df$t1[2] - distance_df$t1[1]
      distance <- sum(delta * distance_df$d1)
      
      data_grouped$weighted_pitch_trajectory_distance[i] <- distance
    }
  }
  
  data_grouped$weighted_pitch_trajectory_distance[data_grouped$pitch_number == 1] <- NA
  data_grouped <- data_grouped %>%
    group_by(game_pk, at_bat_number) %>%
    mutate(
      skip_occurred = pitch_number - lag(pitch_number, default = 0) > 1,
      weighted_pitch_trajectory_distance = ifelse(skip_occurred, NA, weighted_pitch_trajectory_distance)
    ) %>%
    ungroup()
  data_grouped <- data_grouped %>% select(-skip_occurred)
  
  # Calculate R-squared
  r_squared <- cor(data_grouped$weighted_pitch_trajectory_distance, data_grouped$swing_decision_rv, use = "complete.obs")^2
  return(r_squared)
}

# Optimization using optim (more robust than brute force)
optimize_weights <- function(data) {
  optimization_result <- optim(par = c(1, 1, 1), # Initial weights
                               fn = function(weights) {
                                 -calculate_r_squared(data, weights[1], weights[2], weights[3]) # Negative because optim minimizes
                               },
                               lower = c(0, 0, 0), # Weights must be non-negative
                               upper = c(10,10,10),
                               method = "L-BFGS-B") # Use bounded optimization
  
  best_weights <- optimization_result$par
  best_r_squared <- -optimization_result$value
  return(list(weights = best_weights, r_squared = best_r_squared))
}

# Run optimization
optimized_weights <- optimize_weights(pbp_2022_2024_subset)

# Print results
print(paste("Optimized X weight:", optimized_weights$weights[1]))
print(paste("Optimized Y weight:", optimized_weights$weights[2]))
print(paste("Optimized Z weight:", optimized_weights$weights[3]))
print(paste("Best R-squared:", optimized_weights$r_squared))



#===============================================================================
# Weighted Pitch Trajectory Calculation
#==============================================================================

# Print command for data updates
for (i in 1:nrow(pbp_2022_2024_grouped)) {
  if (i %% 1000 == 0) {
    message(paste("Processing row", i, "of", nrow(pbp_2022_2024_grouped)))
  }
  # Creating individual df for pitch
  at_bat_data <- pbp_2022_2024_grouped[i, ]
  pitch_data <- at_bat_data$trajectory
  pitch_df <- data.frame(
    t = seq(0, at_bat_data$reaction_zone, length.out = 50),
    x = pitch_data[[1]]$x,
    y = pitch_data[[1]]$y,
    z = pitch_data[[1]]$z
  )
  # Creating individual df for previous pitch
  if (i > 1) {
    prev_pitch_data <- pbp_2022_2024_grouped[i - 1, ]$trajectory
    prev_pitch_df <- data.frame(
      t = seq(0, pbp_2022_2024_grouped[i - 1, ]$reaction_zone, length.out = 50), # Using time of flight inside hitter reaction zone
      x = prev_pitch_data[[1]]$x,
      y = prev_pitch_data[[1]]$y,
      z = prev_pitch_data[[1]]$z
    )
    #Output 200 points for each pitch to be compared
    t1 <- seq(0, max(c(pitch_df$t, prev_pitch_df$t)), length.out = 200)
    x1 <- rep(0, 200)
    x2 <- rep(0, 200)
    y1 <- rep(0, 200)
    y2 <- rep(0, 200)
    z1 <- rep(0, 200)
    z2 <- rep(0, 200)
    h1 <- 0.01
    h2 <- 0.01
    #Calculate gaussian smoothing to create 200 points on the pitches trajectory
    for (j in 1:200) {
      x1[j] <- sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2) * pitch_df$x) /
        sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2))
      y1[j] <- sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2) * pitch_df$y) /
        sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2))
      z1[j] <- sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2) * pitch_df$z) /
        sum(exp(-0.5 * ((t1[j] - pitch_df$t) / h1)^2))
      x2[j] <- sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2) * prev_pitch_df$x) /
        sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2))
      y2[j] <- sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2) * prev_pitch_df$y) /
        sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2))
      z2[j] <- sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2) * prev_pitch_df$z) /
        sum(exp(-0.5 * ((t1[j] - prev_pitch_df$t) / h2)^2))
    }
    
    # Calculate euclidean distances between two points (with weights)
    distance_df <- data.frame(t1, d1 = rep(0, 200), x1, y1, z1, x2, y2, z2)
    distance_df$d1 <- sqrt(
      optimized_weights[1] * (distance_df$x1 - distance_df$x2)^2 +
        optimized_weights[2] * (distance_df$y1 - distance_df$y2)^2 +
        optimized_weights[3] * (distance_df$z1 - distance_df$z2)^2
    )
    delta <- distance_df$t1[2] - distance_df$t1[1]
    distance <- sum(delta * distance_df$d1)
    
    pbp_2022_2024_grouped$weighted_trajectory_distance[i] <- distance
  }
}

pbp_2022_2024 <- pbp_2022_2024_grouped  

# Update pitch_trajectory to NA where pitch_number is 1 
pbp_2022_2024$weighted_trajectory_distance[pbp_2022_2024$pitch_number == 1] <- NA 

# Identify skips in pitch_number within each plate appearance 
pbp_2022_2024 <- pbp_2022_2024 %>%
  arrange(game_pk, at_bat_number) %>%
  group_by(game_pk, at_bat_number) %>%  
  mutate(     skip_occurred = pitch_number - lag(pitch_number, default = 0) > 1,     
              weighted_trajectory_distance = ifelse(skip_occurred, NA, weighted_trajectory_distance)   ) %>%   
  ungroup()  

# Drop the temporary column
pbp_2022_2024 <- pbp_2022_2024 %>% select(-skip_occurred)                       
                                
                                
                                
