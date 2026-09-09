#===============================================================================
# Libraries
#===============================================================================
#download necessary libraries
install.packages("tidyverse")
install.packages("lightgbm")
install.packages("xgboost")
install.packages("rBayesianOptimization")
install.packages("caret")
install.packages("pROC")
install.packages("caTools")

# Load necessary libraries
library(tidyverse)
library(caTools)
library(lightgbm)
library(xgboost)
library(pROC)
library(rBayesianOptimization)
library(xgboost)
library(caret)

#===============================================================================
# Data Prep
#===============================================================================

pbp_2022_2024 <- pbp_2022_2024 %>%
  arrange(game_pk, inning, inning_topbot, at_bat_number, pitch_number) %>%
  mutate(row_index = row_number())

pbp_2022_2024 <- pbp_2022_2024 %>% drop_na(pitch_type)


# One-hot encode pitch_type
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(pitch_type = as.factor(pitch_type))

# Create one-hot encoding matrix for pitch_type
pitch_type_onehot <- model.matrix(~ pitch_type - 1, data = pbp_2022_2024)

calculate_VAA <- function(data, vy0_col, ay_col, vz0_col, az_col) {
  # Constants
  y0 <- 50                 # Initial y-position in feet
  yf <- 17 / 12            # Final y-position (home plate) in feet
  
  # Ensure required columns exist in the dataset
  required_cols <- c(vy0_col, ay_col, vz0_col, az_col)
  if (!all(required_cols %in% names(data))) {
    stop("The dataset must contain all specified columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Perform calculations
  data$vy_f <- -sqrt(data[[vy0_col]]^2 - (2 * data[[ay_col]] * (y0 - yf)))
  data$t <- (data$vy_f - data[[vy0_col]]) / data[[ay_col]]
  data$vz_f <- data[[vz0_col]] + (data[[az_col]] * data$t)
  data$VAA <- -atan(data$vz_f / data$vy_f) * (180 / pi)
  
  # Return updated dataset
  return(data)
}

calculate_HAA <- function(data, vy0_col, ay_col, vx0_col, ax_col) {
  # Constants
  y0 <- 50                 # Initial y-position in feet
  yf <- 17 / 12            # Final y-position (home plate) in feet
  
  # Ensure required columns exist in the dataset
  required_cols <- c(vy0_col, ay_col, vx0_col, ax_col)
  if (!all(required_cols %in% names(data))) {
    stop("The dataset must contain all specified columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Perform calculations
  data$vy_f <- -sqrt(data[[vy0_col]]^2 - (2 * data[[ay_col]] * (y0 - yf)))
  data$t <- (data$vy_f - data[[vy0_col]]) / data[[ay_col]]
  data$vx_f <- data[[vx0_col]] + (data[[ax_col]] * data$t)
  data$HAA <- -atan(data$vx_f / data$vy_f) * (180 / pi)
  
  # Return updated dataset
  return(data)
}


# Apply the function to the dataset `pbp_2022_2024`
pbp_2022_2024 <- calculate_VAA(
  data = pbp_2022_2024,
  vy0_col = "vy0",
  ay_col = "ay",
  vz0_col = "vz0",
  az_col = "az"
)



# Combine one-hot encoded pitch_type with plate_z and VAA
model_data <- cbind(pbp_2022_2024 %>% select(plate_z, VAA), pitch_type_onehot)

# Fit the linear model
VAA_model <- lm(VAA ~ plate_z + ., data = model_data)

# Predict expected VAA using the linear model
pbp_2022_2024$xVAA <- predict(VAA_model, newdata = model_data)

# One-hot encode pitch_type
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(pitch_type = as.factor(pitch_type))

# Apply the calculate_HAA function to the dataset
pbp_2022_2024 <- calculate_HAA(
  data = pbp_2022_2024,
  vy0_col = "vy0",
  ay_col = "ay",
  vx0_col = "vx0",
  ax_col = "ax"
)

# Create one-hot encoding matrix for pitch_type
pitch_type_onehot <- model.matrix(~ pitch_type, data = pbp_2022_2024)

# Combine one-hot encoded pitch_type with plate_z and VAA
model_data <- cbind(pbp_2022_2024 %>% select(plate_x, HAA), pitch_type_onehot)

# Fit the linear model
HAA_model <- lm(HAA ~ plate_x + ., data = model_data)

# Predict expected VAA using the linear model
pbp_2022_2024$xHAA <- predict(HAA_model, newdata = model_data)

# Add LIVAA and LIHAA
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(LIVAA = xVAA - VAA) %>%
  mutate(LIHAA = xHAA - HAA)

# Perform a left join to add the height column
pbp_2022_2024 <- merge(pbp_2022_2024, 
                       player_heights, 
                       by.x = "pitcher",  # Column in pbp_2022_2024 to match
                       by.y = "id",       # Column in player_heights to match
                       all.x = TRUE)      # Keeps all rows in pbp_2022_2024
colnames(pbp_2022_2024)[colnames(pbp_2022_2024) == "height.x"] <- "height"
# Convert the height column to numeric
pbp_2022_2024$height <- as.numeric(pbp_2022_2024$height)

# Update the calculate_arm_angle function
calculate_arm_angle <- function(data, release_pos_x, release_pos_z, player_height) {
  
  # Ensure required columns exist in the dataset
  required_cols <- c(release_pos_x, release_pos_z, player_height)
  
  if (!all(required_cols %in% names(data))) {
    stop("The dataset must contain all specified columns: ", paste(required_cols, collapse = ", "))
  }
  
  # Extract relevant columns from data
  x <- data[[release_pos_x]]
  z <- data[[release_pos_z]]
  height <- data[[player_height]]
  
  # Calculate components for arm angle
  adjacent <- z - (height * 0.7)  # Adjusted to proper scaling
  opposite <- abs(x)
  hypotenuse <- sqrt(opposite^2 + adjacent^2)
  
  # Calculate arm angle using cosine rule
  data$arm_angle <- acos(adjacent / hypotenuse) * (180 / pi)  # Convert to degrees
  
  # Return updated dataset
  return(data)
}


pbp_2022_2024 <- calculate_arm_angle(
  data = pbp_2022_2024,
  release_pos_x = "release_pos_x",
  release_pos_z = "release_pos_z",
  player_height = "height"
)

# Fit the models
model_az <- lm(arm_angle ~ az, data = pbp_2022_2024)
model_ax <- lm(arm_angle ~ ax, data = pbp_2022_2024)

# Add predicted values as new columns
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    x_az = predict(model_az, newdata = pbp_2022_2024),
    x_ax = predict(model_ax, newdata = pbp_2022_2024)
  )

# Add LIVAA and LIHAA
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(az_AA = x_az - az) %>%
  mutate(ax_AA = x_ax - ax)

pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(batted_ball_type = case_when(
    description == "hit_into_play" & launch_angle < 10 ~ "GroundBall",
    description == "hit_into_play" & launch_angle >= 10 & launch_angle < 25 ~ "LineDrive",
    description == "hit_into_play" & launch_angle >= 25 & launch_angle < 50 ~ "FlyBall",
    description == "hit_into_play" & launch_angle >= 50 ~ "PopUp",
    TRUE ~ NA_character_  # Keep other rows as NA without removing them
  ))

# Create one-hot encoding columns for each batted ball type
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    GroundBall = ifelse(batted_ball_type == "GroundBall", 1, 0),
    LineDrive = ifelse(batted_ball_type == "LineDrive", 1, 0),
    FlyBall = ifelse(batted_ball_type == "FlyBall", 1, 0),
    PopUp = ifelse(batted_ball_type == "PopUp", 1, 0)
  )


#===============================================================================
# PITCHES CALLED MODEL
#===============================================================================
# Define the vector of pitch calls
pitch_call <- c("ball", "called_strike", "hit_by_pitch", "ball_in_dirt")

# Filter and preprocess the dataset
pbp_2022_2024_pitch_call <- filter(pbp_2022_2024, description %in% pitch_call) %>%
  mutate(description = case_when(
    description == "ball" ~ 0,
    description == "called_strike" ~ 1,
    description == "blocked_ball" ~ 2,
    description == "blocked_ball" ~ 0
  )) %>%
  filter(!is.na(description))

# Prepare features and labels
pbp_2022_2024_pitch_call <- pbp_2022_2024_pitch_call %>%
  select(plate_x, plate_z, row_index, description)

# Split data into training and testing sets
set.seed(123)
train_indices <- sample(1:nrow(pbp_2022_2024_pitch_call), 
                        size = floor(0.75 * nrow(pbp_2022_2024_pitch_call)))
train_set <- pbp_2022_2024_pitch_call[train_indices, ]
test_set <- pbp_2022_2024_pitch_call[-train_indices, ]

# Extract labels and features
y_train <- train_set$description
y_test <- test_set$description
X_train <- train_set %>% select(-description, -row_index)
X_test <- test_set %>% select(-description, -row_index)

# Convert to DMatrix for XGBoost
dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
dtest <- xgb.DMatrix(data = as.matrix(X_test), label = y_test)

# Define XGBoost parameters
params <- list(
  objective = "multi:softprob",
  eval_metric = "mlogloss",
  num_class = 3,
  eta = 0.1,
  max_depth = 6,
  gamma = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8
)

# Train the model with early stopping
xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 100,
  watchlist = list(train = dtrain, eval = dtest),
  early_stopping_rounds = 10,
  verbose = 1
)

# Apply the model to the full pbp_2022_2024 dataset
pbp_2022_2024_features <- pbp_2022_2024 %>%
  select(plate_x, plate_z, row_index) %>%
  filter(!is.na(plate_x) & !is.na(plate_z)) # Ensure no missing values in features

dpbp <- xgb.DMatrix(data = as.matrix(pbp_2022_2024_features %>% select(-row_index)))

# Predict probabilities
pred_probs <- predict(xgb_model, dpbp)
pred_probs <- matrix(pred_probs, ncol = 3, byrow = TRUE)

# Add predicted probabilities to pbp_2022_2024
# Create a data frame with predicted probabilities and their corresponding indices
predicted_probs <- data.frame(
  row_index = pbp_2022_2024_features$row_index,
  xBall = pred_probs[, 1],
  xStrike = pred_probs[, 2],
  xHBP = pred_probs[, 3]
)

# Join the predictions back to pbp_2022_2024 by row_index
pbp_2022_2024 <- pbp_2022_2024 %>%
  left_join(predicted_probs, by = "row_index")


rm(pbp_2022_2024_pitch_call)
rm(pbp_2022_2024_features)
rm(dpbp)
rm(predictions_by_type)
rm(predicted_probs)
rm(train_set)
rm(test_set)
rm(VAA_model)
rm(X_test)
rm(X_train)
rm(xgb_model)
rm(pred_probs)
rm(evaluation_summary)
rm(pitch_type_onehot)
rm(roc_ball)
rm(roc_strike)
rm(pred_probs_test)

#===============================================================================
# PITCHES SWUNG AT MODEL
#===============================================================================
# Define the vector of swing types
swing_types <- c("swinging_strike", "foul_tip", "foul", "hit_into_play", "swinging_strike_blocked")

# Filter the data for the specified pitch calls to create training dataset
filtered_pbp_2022_2024_swing <- filter(pbp_2022_2024, description %in% swing_types)

# Create the 'contact' column based on the description for training data
pbp_2022_2024_swing <- filtered_pbp_2022_2024_swing %>%
  mutate(whiff = if_else(description %in% c("foul", "hit_into_play"), 0, 1))

# Prepare features for training
one_hot_pbp_2022_2024_swing <- pbp_2022_2024_swing %>%
  select(pitch_type, p_throws, balls, strikes, p_throws, stand, release_spin_rate, spin_axis, release_speed, az, ax, 
         release_pos_x, release_pos_z, release_extension, 
         plate_x, plate_z, row_index, LIHAA, LIVAA, arm_angle, whiff)

# Transform categorical variables
one_hot_pbp_2022_2024_swing <- one_hot_pbp_2022_2024_swing %>%
  mutate(
    stand = ifelse(stand == "Right", 1, 0),
    p_throws = ifelse(p_throws == "Right", 1, 0)
  )

# One-hot encode pitch_type
pitch_type_onehot_matrix <- model.matrix(~ pitch_type - 1, one_hot_pbp_2022_2024_swing)
one_hot_pbp_2022_2024_swing <- cbind(one_hot_pbp_2022_2024_swing, pitch_type_onehot_matrix) %>%
  select(-pitch_type)

# Prepare features for full dataset (excluding whiff column)
one_hot_pbp_2022_2024_full <- pbp_2022_2024 %>%
  select(pitch_type, p_throws, balls, strikes, p_throws, stand, release_spin_rate, spin_axis, release_speed, az, ax, 
         release_pos_x, release_pos_z, release_extension, 
         plate_x, plate_z, row_index, LIHAA, LIVAA, arm_angle)

# Transform categorical variables for full dataset
one_hot_pbp_2022_2024_full <- one_hot_pbp_2022_2024_full %>%
  mutate(
    stand = ifelse(stand == "Right", 1, 0),
    p_throws = ifelse(p_throws == "Right", 1, 0)
  )

# One-hot encode pitch_type for full dataset
pitch_type_onehot_matrix_full <- model.matrix(~ pitch_type - 1, one_hot_pbp_2022_2024_full)
one_hot_pbp_2022_2024_full <- cbind(one_hot_pbp_2022_2024_full, pitch_type_onehot_matrix_full) %>%
  select(-pitch_type)

# Ensure full dataset has the same columns as training dataset
# Identify columns present in training dataset but not in full dataset
training_columns <- colnames(one_hot_pbp_2022_2024_swing %>% select(-whiff))
missing_columns <- setdiff(training_columns, colnames(one_hot_pbp_2022_2024_full))

# Add missing columns to full dataset with 0 values if necessary
for (col in missing_columns) {
  one_hot_pbp_2022_2024_full[[col]] <- 0
}

# Ensure columns are in the same order
one_hot_pbp_2022_2024_full <- one_hot_pbp_2022_2024_full[, training_columns]

# Prepare training data
X_train <- one_hot_pbp_2022_2024_swing %>% select(-whiff)
y_train <- one_hot_pbp_2022_2024_swing$whiff

# Convert to matrix format
X_train_matrix <- as.matrix(X_train)
X_full_matrix <- as.matrix(one_hot_pbp_2022_2024_full)

# Create LightGBM dataset for training
lgb_train <- lgb.Dataset(data = X_train_matrix, label = y_train)

# Modified Bayesian Optimization Function
whiff_objective <- function(num_leaves, max_depth, learning_rate) {
  # Prepare parameters
  params <- list(
    objective = "binary",
    metric = "binary_logloss",
    num_leaves = round(num_leaves),
    max_depth = round(max_depth),
    learning_rate = learning_rate,
    min_data_in_leaf = 30,  # Fixed value instead of dynamic
    feature_pre_filter = FALSE,  # Add this to allow dynamic parameter changes
    feature_fraction = 0.8,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    scale_pos_weight = (length(y_train) - sum(y_train)) / sum(y_train)
  )
  
  # Perform cross-validation
  cv_results <- lgb.cv(
    params = params,
    data = lgb_train,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 50,
    verbose = -1
  )
  
  # Return best score (minimize binary logloss)
  best_score <- min(cv_results$best_score)
  return(list(Score = -best_score))  # Negative because optimizer maximizes
}

# Updated parameter search bounds
bounds <- list(
  num_leaves = c(20L, 100L),
  max_depth = c(3L, 15L),
  learning_rate = c(0.01, 0.3)
)

# Perform Bayesian Optimization
opt_results <- BayesianOptimization(
  FUN = whiff_objective,
  bounds = bounds,
  init_points = 10,
  n_iter = 50,
  acq = "ucb"
)

# Extract best parameters
best_params <- list(
  objective = "binary",
  metric = "binary_logloss",
  num_leaves = round(opt_results$Best_Par["num_leaves"]),
  max_depth = round(opt_results$Best_Par["max_depth"]),
  learning_rate = opt_results$Best_Par["learning_rate"],
  min_data_in_leaf = 30,  # Fixed value
  feature_pre_filter = FALSE,
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq = 1,
  scale_pos_weight = (length(y_train) - sum(y_train)) / sum(y_train)
)

# Train final model with optimized parameters
final_model <- lgb.train(
  params = best_params,
  data = lgb_train,
  nrounds = 1000,
  verbose = 0
)

# Make predictions on full dataset
pred_probs <- predict(final_model, X_full_matrix)

# Add predictions to the original dataset
pbp_2022_2024$xWhiff <- pred_probs

# Optional: ROC Curve for Training Data
# First, create a test set from the training data
set.seed(123)
sample_split <- sample.split(Y = y_train, SplitRatio = 0.75)
train_subset <- subset(one_hot_pbp_2022_2024_swing, sample_split == TRUE)
test_subset <- subset(one_hot_pbp_2022_2024_swing, sample_split == FALSE)

X_test_subset <- test_subset %>% select(-whiff)
y_test_subset <- test_subset$whiff
X_test_matrix_subset <- as.matrix(X_test_subset)

# Predictions on test subset
test_pred_probs <- predict(final_model, as.matrix(X_test_subset))

# ROC Curve
roc_curve <- roc(y_test_subset, test_pred_probs)
auc_value <- auc(roc_curve)

# Plot ROC Curve
roc_data <- data.frame(
  tpr = rev(roc_curve$sensitivities),
  fpr = rev(1 - roc_curve$specificities)
)

ggplot(roc_data, aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve for Whiff Prediction Model",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  annotate("text", x = 0.7, y = 0.2, 
           label = paste("AUC =", round(auc_value, 4)), 
           color = "red") +
  theme_minimal()

# Feature Importance Plot
importance_matrix <- lgb.importance(final_model, percentage = TRUE)
importance_matrix <- importance_matrix %>%
  rename(Feature = Feature, Importance = Gain)

ggplot(data = importance_matrix, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(
    title = "Feature Importance",
    x = "Features",
    y = "Importance (Gain)"
  ) +
  theme_minimal()

rm(pbp_2022_2024_swing)
rm(filtered_pbp_2022_2024_swing)
rm(one_hot_pbp_2022_2024_swing)
rm(pitch_type_onehot_matrix)
rm(one_hot_pbp_2022_2024_full)
rm(training_columns)
rm(missing_columns)
rm(X_train)
rm(y_train)
rm(X_train_matrix)
rm(X_full_matrix)
rm(train_subset)
rm(test_subset)
rm(X_test_subset)
rm(y_test_subset)
rm(X_test_matrix_subset)
rm(test_pred_probs)
rm(pred_probs)
rm(best_params)
rm(bounds)
rm(final_model)
rm(importance_matrix)
rm(lgb_train)
rm(opt_results)
rm(params)
rm(pitch_type_onehot_matrix_full)
rm(sample_split)
rm(col)
rm(dtest)
rm(dtrain)
rm(train_indices)
rm(y_test)
rm(best_params)




#===============================================================================
# PITCHES MADE CONTACT WITH MODEL
#===============================================================================
# Define the vector of pitch calls
contact_types <- c("foul", "hit_into_play")

# Filter the data for the specified pitch calls to create training dataset
filtered_pbp_2022_2024_contact <- filter(pbp_2022_2024, description %in% contact_types)

# Create the 'contact' column based on the description for training data
pbp_2022_2024_contact <- filtered_pbp_2022_2024_contact %>%
  mutate(in_play = if_else(description == "foul", 0, 1))

# Prepare features for training
one_hot_pbp_2022_2024_contact <- pbp_2022_2024_contact %>%
  select(balls, strikes, p_throws, stand, release_spin_rate, spin_axis, release_speed, az, ax, 
         release_pos_x, release_pos_z, , 
         plate_x, plate_z, row_index, arm_angle, LIHAA, LIVAA,  in_play)

# Transform categorical variables
one_hot_pbp_2022_2024_contact <- one_hot_pbp_2022_2024_contact %>%
  mutate(
    stand = ifelse(stand == "Right", 1, 0),
    p_throws = ifelse(p_throws == "Right", 1, 0)
  )

# Prepare features for full dataset (excluding in_play column)
one_hot_pbp_2022_2024_full <- pbp_2022_2024 %>%
  select(balls, strikes, p_throws, stand, release_spin_rate, spin_axis, release_speed, az, ax, 
         release_pos_x, release_pos_z, release_extension, 
         plate_x, plate_z, row_index, arm_angle, LIHAA, LIVAA)

# Transform categorical variables for full dataset
one_hot_pbp_2022_2024_full <- one_hot_pbp_2022_2024_full %>%
  mutate(
    stand = ifelse(stand == "Right", 1, 0),
    p_throws = ifelse(p_throws == "Right", 1, 0)
  )

# Ensure full dataset has the same columns as training dataset
# Identify columns present in training dataset but not in full dataset
training_columns <- colnames(one_hot_pbp_2022_2024_contact %>% select(-in_play))
missing_columns <- setdiff(training_columns, colnames(one_hot_pbp_2022_2024_full))

# Add missing columns to full dataset with 0 values if necessary
for (col in missing_columns) {
  one_hot_pbp_2022_2024_full[[col]] <- 0
}

# Ensure columns are in the same order
one_hot_pbp_2022_2024_full <- one_hot_pbp_2022_2024_full[, training_columns]

# Prepare training data
X_train <- one_hot_pbp_2022_2024_contact %>% select(-in_play)
y_train <- one_hot_pbp_2022_2024_contact$in_play

# Convert to matrix format
X_train_matrix <- as.matrix(X_train)
X_full_matrix <- as.matrix(one_hot_pbp_2022_2024_full)

# Create LightGBM dataset for training
lgb_train <- lgb.Dataset(data = X_train_matrix, label = y_train)

# Define the Objective Function for Bayesian Optimization
lgb_opt_function <- function(num_leaves, max_depth, learning_rate, feature_fraction, bagging_fraction) {
  
  params <- list(
    objective = "binary",
    metric = "binary_logloss",
    num_leaves = round(num_leaves),
    max_depth = round(max_depth),
    learning_rate = learning_rate,
    feature_fraction = feature_fraction,
    bagging_fraction = bagging_fraction,
    bagging_freq = 1,
    min_data_in_leaf = 20,
    scale_pos_weight = (length(y_train) - sum(y_train)) / sum(y_train)
  )
  
  lgb_cv <- lgb.cv(
    params = params,
    data = lgb_train,
    nrounds = 100,
    nfold = 5,
    verbose = -1,
    stratified = TRUE,
    early_stopping_rounds = 10
  )
  
  # Extract the minimum binary_logloss value
  binary_logloss_vals <- unlist(lgb_cv$record_evals$valid$binary_logloss$eval)
  min_logloss <- min(binary_logloss_vals, na.rm = TRUE)
  
  list(Score = -min_logloss, Pred = NA)
}

# Perform Bayesian Optimization
opt_results <- BayesianOptimization(
  FUN = lgb_opt_function,
  bounds = list(
    num_leaves = c(20L, 50L),
    max_depth = c(5L, 15L),
    learning_rate = c(0.01, 0.2),
    feature_fraction = c(0.6, 1.0),
    bagging_fraction = c(0.6, 1.0)
  ),
  init_points = 10,
  n_iter = 40,
  acq = "ucb",
  kappa = 2.576,
  verbose = TRUE
)

# Get the best parameters from the Bayesian Optimization result
best_params <- opt_results$Best_Par

# Extract the parameters as a flat list
best_params <- unlist(best_params)

# Ensure the parameters are in the correct format for LightGBM
params <- list(
  objective = "binary",
  metric = "binary_logloss",
  num_leaves = round(best_params['num_leaves']),
  max_depth = round(best_params['max_depth']),
  learning_rate = best_params['learning_rate'],
  feature_fraction = best_params['feature_fraction'],
  bagging_fraction = best_params['bagging_fraction'],
  bagging_freq = 1,
  min_data_in_leaf = 20,
  scale_pos_weight = (length(y_train) - sum(y_train)) / sum(y_train)
)

# Train the final model using the optimized parameters
final_model <- lgb.train(
  params = params,
  data = lgb_train,
  nrounds = 1000,
  verbose = 0
)

# Make predictions on the full dataset
pred_probs <- predict(final_model, X_full_matrix)

# Add predictions to the original dataset
pbp_2022_2024$xIn_Play <- pred_probs


# Optional: ROC Curve for Training Data
# First, create a test set from the training data
set.seed(123)
sample_split <- sample.split(Y = y_train, SplitRatio = 0.75)
train_subset <- subset(one_hot_pbp_2022_2024_contact, sample_split == TRUE)
test_subset <- subset(one_hot_pbp_2022_2024_contact, sample_split == FALSE)

X_test_subset <- test_subset %>% select(-in_play)
y_test_subset <- test_subset$in_play
X_test_matrix_subset <- as.matrix(X_test_subset)

# Predictions on test subset
test_pred_probs <- predict(final_model, X_test_matrix_subset)

# ROC Curve
roc_curve <- roc(y_test_subset, test_pred_probs)
auc_value <- auc(roc_curve)

# Plot ROC Curve
roc_data <- data.frame(
  tpr = rev(roc_curve$sensitivities),
  fpr = rev(1 - roc_curve$specificities)
)

ggplot(roc_data, aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve for Contact Type Prediction Model",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  annotate("text", x = 0.7, y = 0.2, 
           label = paste("AUC =", round(auc_value, 4)), 
           color = "red") +
  theme_minimal()

# Feature Importance Plot
importance_matrix <- lgb.importance(final_model, percentage = TRUE)
importance_matrix <- importance_matrix %>%
  rename(Feature = Feature, Importance = Gain)

ggplot(data = importance_matrix, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(
    title = "Feature Importance",
    x = "Features",
    y = "Importance (Gain)"
  ) +
  theme_minimal()

rm(pbp_2022_2024_contact)
rm(filtered_pbp_2022_2024_contact)
rm(one_hot_pbp_2022_2024_contact)
rm(pitch_type_onehot_matrix)
rm(one_hot_pbp_2022_2024_full)
rm(training_columns)
rm(missing_columns)
rm(X_train)
rm(y_train)
rm(X_train_matrix)
rm(X_full_matrix)
rm(train_subset)
rm(test_subset)
rm(X_test_subset)
rm(y_test_subset)
rm(X_test_matrix_subset)
rm(test_pred_probs)
rm(pred_probs)
rm(roc_data)
rm(roc_curve)
rm(lgb_train)

#===============================================================================
# HIT TYPE PREDICTION MODEL
#===============================================================================

# Define target columns
target_columns <- c("GroundBall", "LineDrive", "FlyBall", "PopUp")

# Remove rows with NA in GroundBall column and filter in-play hits for training dataset
pbp_2022_2024_in_play <- pbp_2022_2024 %>%
  filter(description == "hit_into_play") %>%
  filter(!is.na(GroundBall))

# Data Preparation
# Remove rows with NA in target columns and select relevant features
pbp_2022_2024_in_play <- pbp_2022_2024_in_play %>%
  filter(!is.na(GroundBall) & !is.na(LineDrive) & !is.na(FlyBall) & !is.na(PopUp))
relevant_features <- c("pitch_type", "p_throws", "balls", "strikes", "stand", "release_speed", "az", "ax", 
                       "release_pos_x", "release_pos_z", "release_extension", 
                       "plate_x", "plate_z", "release_spin_rate", "spin_axis", 
                       "LIHAA", "LIVAA", "arm_angle", "row_index")

# Prepare dataset for training
data_prep <- pbp_2022_2024_in_play %>%
  select(all_of(relevant_features), GroundBall, LineDrive, FlyBall, PopUp) %>%
  mutate(
    p_throws = ifelse(p_throws == "Right", 1, 0),
    stand = ifelse(stand == "Right", 1, 0)  # Encode categorical variables
  )

# Prepare dataset for prediction (entire dataset)
prediction_data_prep <- pbp_2022_2024 %>%
  select(all_of(relevant_features)) %>%
  mutate(
    p_throws = ifelse(p_throws == "Right", 1, 0),
    stand = ifelse(stand == "Right", 1, 0)  # Encode categorical variables
  )

# One-hot encode pitch_type for training data
pitch_type_onehot_matrix_train <- model.matrix(~ pitch_type - 1, data_prep)
data_prep <- cbind(data_prep, pitch_type_onehot_matrix_train) %>%
  select(-pitch_type)

# One-hot encode pitch_type for prediction data
pitch_type_onehot_matrix_prediction <- model.matrix(~ pitch_type - 1, prediction_data_prep)
prediction_data_prep <- cbind(prediction_data_prep, pitch_type_onehot_matrix_prediction) %>%
  select(-pitch_type)

# Ensure prediction data has the same columns as training data
missing_cols <- setdiff(names(data_prep)[!names(data_prep) %in% target_columns], 
                        names(prediction_data_prep))
for (col in missing_cols) {
  prediction_data_prep[[col]] <- 0
}

# Reorder columns to match training data
prediction_data_prep <- prediction_data_prep[, names(data_prep)[!names(data_prep) %in% target_columns]]

# Split training data into train/test sets
set.seed(123)
sample_split <- sample.split(Y = data_prep$GroundBall, SplitRatio = 0.75)
train_set <- subset(data_prep, sample_split == TRUE)
test_set <- subset(data_prep, sample_split == FALSE)

# Prepare for model training
X_train <- train_set %>% select(-all_of(target_columns))
X_test <- test_set %>% select(-all_of(target_columns))

# Updated train_model function with extensive error checking
train_model <- function(target, X_train, X_test, train_set, test_set) {
  # Check for valid data
  if (nrow(X_train) == 0 || nrow(train_set) == 0) {
    stop(paste("No training data for", target))
  }
  
  y_train <- train_set[[target]]
  y_test <- test_set[[target]]
  
  # Validate target variable
  if (length(unique(y_train)) < 2) {
    stop(paste("Insufficient variation in target variable", target))
  }
  
  # LightGBM Dataset preparation
  lgb_train <- lgb.Dataset(data = as.matrix(X_train), label = y_train)
  lgb_test <- lgb.Dataset(data = as.matrix(X_test), label = y_test, reference = lgb_train)
  
  # Objective function for Bayesian Optimization
  objective_function <- function(num_leaves, max_depth, learning_rate, min_data_in_leaf) {
    params <- list(
      objective = "binary",
      metric = "binary_logloss",
      num_leaves = round(num_leaves),
      learning_rate = learning_rate,
      max_depth = round(max_depth),
      min_data_in_leaf = round(min_data_in_leaf),
      feature_fraction = 0.8,
      bagging_fraction = 0.8,
      bagging_freq = 1,
      scale_pos_weight = (length(y_train) - sum(y_train)) / sum(y_train),
      feature_pre_filter = FALSE
    )
    
    tryCatch({
      model <- lgb.train(
        params = params,
        data = lgb_train,
        nrounds = 1000,
        verbose = -1
      )
      
      pred_probs <- predict(model, as.matrix(X_test))
      roc_curve <- roc(y_test, pred_probs)
      auc_value <- auc(roc_curve)
      
      return(list(Score = auc_value))
    }, error = function(e) {
      cat("Error in optimization for", target, ":", e$message, "\n")
      return(list(Score = 0))
    })
  }
  
  # Bayesian Optimization
  bounds <- list(
    num_leaves = c(20L, 50L),
    max_depth = c(5L, 15L),
    learning_rate = c(0.01, 0.2),
    min_data_in_leaf = c(10L, 30L)
  )
  
  opt_results <- tryCatch({
    BayesianOptimization(
      FUN = objective_function,
      bounds = bounds,
      init_points = 10,
      n_iter = 40,
      acq = "ucb"
    )
  }, error = function(e) {
    cat("Bayesian Optimization failed for", target, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(opt_results)) {
    stop(paste("Optimization failed for", target))
  }
  
  # Train final model using best parameters
  best_params <- list(
    objective = "binary",
    metric = "binary_logloss",
    num_leaves = opt_results$Best_Par[["num_leaves"]],
    learning_rate = opt_results$Best_Par[["learning_rate"]],
    max_depth = opt_results$Best_Par[["max_depth"]],
    min_data_in_leaf = opt_results$Best_Par[["min_data_in_leaf"]],
    feature_fraction = 0.8,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    scale_pos_weight = (length(y_train) - sum(y_train)) / sum(y_train),
    feature_pre_filter = FALSE
  )
  
  final_model <- tryCatch({
    lgb.train(
      params = best_params,
      data = lgb_train,
      nrounds = 1000,
      verbose = -1
    )
  }, error = function(e) {
    cat("Model training failed for", target, ":", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(final_model)) {
    stop(paste("Failed to train model for", target))
  }
  
  # Predictions and Evaluation
  pred_probs <- predict(final_model, as.matrix(X_test))
  pred_classes <- ifelse(pred_probs > 0.5, 1, 0)
  
  auc <- auc(roc(y_test, pred_probs))
  accuracy <- sum(pred_classes == y_test) / length(y_test)
  
  return(list(
    model = final_model,
    auc = auc,
    accuracy = accuracy,
    predictions = list(test_probs = pred_probs, test_classes = pred_classes)
  ))
}

# Train models for all target variables and store results
results <- list()
for (target in target_columns) {
  cat("Training model for:", target, "\n")
  results[[target]] <- train_model(target, X_train, X_test, train_set, test_set)
}

# Predict probabilities for the entire dataset
predictions <- list()
for (target in target_columns) {
  cat("Predicting for:", target, "\n")
  predictions[[target]] <- predict(results[[target]]$model, as.matrix(prediction_data_prep))
}

# Add predictions to the original pbp_2022_2024 dataset
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    xGB = predictions$GroundBall,
    xLD = predictions$LineDrive,
    xFB = predictions$FlyBall,
    xPU = predictions$PopUp
  )

# Normalize predictions so they sum to 1
pbp_2022_2024 <- pbp_2022_2024 %>%
  rowwise() %>%
  mutate(
    sum_pred = xGB + xLD + xFB + xPU,
    xGB = xGB / sum_pred,
    xLD = xLD / sum_pred,
    xFB = xFB / sum_pred,
    xPU = xPU / sum_pred
  ) %>%
  select(-sum_pred)  # Remove the temporary sum column


# Evaluation Summary
evaluation_summary <- lapply(results, function(res) {
  list(auc = res$auc, accuracy = res$accuracy)
})
print(evaluation_summary)

# Visualization: ROC Curve for GroundBall
roc_curve <- roc(test_set$GroundBall, results$GroundBall$predictions$test_probs)
ggplot(data.frame(tpr = rev(roc_curve$sensitivities), fpr = rev(roc_curve$specificities)), 
       aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve for Hit Type Prediction",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  annotate("text", x = 0.7, y = 0.2, label = paste("AUC =", round(roc_curve$auc, 4)), color = "red") +
  theme_minimal()

# Run Values by hit type
#Flyball      0.586
#Line Drive   0.528
#Ground Ball  0.164
#PopUp        0.0186

rm(pbp_2022_2024_in_play)
rm(data_prep)
rm(pitch_type_onehot_matrix)
rm(prediction_data_prep)
rm(training_columns)
rm(missing_columns)
rm(X_train)
rm(y_train)
rm(X_train_matrix)
rm(X_full_matrix)
rm(train_subset)
rm(test_subset)
rm(X_test_subset)
rm(y_test_subset)
rm(X_test_matrix_subset)
rm(test_pred_probs)
rm(pred_probs)
rm(lgb_train)
rm(lgb_test)
rm(roc_data)
rm(roc_curve)

#===============================================================================
# SWING PREDICTION MODEL
#===============================================================================
# Create a new column `swing` based on the `description` column
pbp_2022_2024$swing <- ifelse(
  pbp_2022_2024$description %in% c("hit_into_play", "foul", 
                                 "swinging_strike", "swinging_strike_blocked", "foul_tip"),
  1,
  0
)


# Prepare features for training
one_hot_pbp_2022_2024 <- pbp_2022_2024 %>%
  select(pitch_type, p_throws, balls, strikes, p_throws, stand, release_spin_rate, spin_axis, release_speed, az, ax, 
         release_pos_x, release_pos_z, release_extension, 
         plate_x, plate_z, row_index, arm_angle, LIHAA, LIVAA, swing)

# Transform categorical variables
one_hot_pbp_2022_2024 <- one_hot_pbp_2022_2024 %>%
  mutate(
    stand = ifelse(stand == "Right", 1, 0),
    p_throws = ifelse(p_throws == "Right", 1, 0)
  )

# One-hot encode pitch_type
pitch_type_onehot_matrix <- model.matrix(~ pitch_type - 1, one_hot_pbp_2022_2024)
one_hot_pbp_2022_2024 <- cbind(one_hot_pbp_2022_2024, pitch_type_onehot_matrix) %>%
  select(-pitch_type)

# Prepare features for full dataset (excluding whiff column)
one_hot_pbp_2022_2024_full <- pbp_2022_2024 %>%
  select(pitch_type, p_throws, balls, strikes, p_throws, stand, release_spin_rate, spin_axis, release_speed, az, ax, 
         release_pos_x, release_pos_z, release_extension, 
         plate_x, plate_z, arm_angle, row_index, LIHAA, LIVAA)

# Transform categorical variables for full dataset
one_hot_pbp_2022_2024_full <- one_hot_pbp_2022_2024_full %>%
  mutate(
    stand = ifelse(stand == "Right", 1, 0),
    p_throws = ifelse(p_throws == "Right", 1, 0)
  )

# One-hot encode pitch_type for full dataset
pitch_type_onehot_matrix_full <- model.matrix(~ pitch_type - 1, one_hot_pbp_2022_2024_full)
one_hot_pbp_2022_2024_full <- cbind(one_hot_pbp_2022_2024_full, pitch_type_onehot_matrix_full) %>%
  select(-pitch_type)

# Ensure full dataset has the same columns as training dataset
# Identify columns present in training dataset but not in full dataset
training_columns <- colnames(one_hot_pbp_2022_2024 %>% select(-swing))
missing_columns <- setdiff(training_columns, colnames(one_hot_pbp_2022_2024_full))

# Add missing columns to full dataset with 0 values if necessary
for (col in missing_columns) {
  one_hot_pbp_2022_2024_full[[col]] <- 0
}

# Ensure columns are in the same order
one_hot_pbp_2022_2024_full <- one_hot_pbp_2022_2024_full[, training_columns]

# Prepare training data
X_train <- one_hot_pbp_2022_2024 %>% select(-swing)
y_train <- one_hot_pbp_2022_2024$swing

# Convert to matrix format
X_train_matrix <- as.matrix(X_train)
X_full_matrix <- as.matrix(one_hot_pbp_2022_2024_full)

# Create LightGBM dataset for training
lgb_train <- lgb.Dataset(data = X_train_matrix, label = y_train)

# Modified Bayesian Optimization Function
whiff_objective <- function(num_leaves, max_depth, learning_rate) {
  # Prepare parameters
  params <- list(
    objective = "binary",
    metric = "binary_logloss",
    num_leaves = round(num_leaves),
    max_depth = round(max_depth),
    learning_rate = learning_rate,
    min_data_in_leaf = 30,  # Fixed value instead of dynamic
    feature_pre_filter = FALSE,  # Add this to allow dynamic parameter changes
    feature_fraction = 0.8,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    scale_pos_weight = (length(y_train) - sum(y_train)) / sum(y_train)
  )
  
  # Perform cross-validation
  cv_results <- lgb.cv(
    params = params,
    data = lgb_train,
    nrounds = 1000,
    nfold = 5,
    early_stopping_rounds = 50,
    verbose = -1
  )
  
  # Return best score (minimize binary logloss)
  best_score <- min(cv_results$best_score)
  return(list(Score = -best_score))  # Negative because optimizer maximizes
}

# Updated parameter search bounds
bounds <- list(
  num_leaves = c(20L, 100L),
  max_depth = c(3L, 15L),
  learning_rate = c(0.01, 0.3)
)

# Perform Bayesian Optimization
opt_results <- BayesianOptimization(
  FUN = whiff_objective,
  bounds = bounds,
  init_points = 10,
  n_iter = 40,
  acq = "ucb"
)

# Extract best parameters
best_params <- list(
  objective = "binary",
  metric = "binary_logloss",
  num_leaves = round(opt_results$Best_Par["num_leaves"]),
  max_depth = round(opt_results$Best_Par["max_depth"]),
  learning_rate = opt_results$Best_Par["learning_rate"],
  min_data_in_leaf = 30,  # Fixed value
  feature_pre_filter = FALSE,
  feature_fraction = 0.8,
  bagging_fraction = 0.8,
  bagging_freq = 1,
  scale_pos_weight = (length(y_train) - sum(y_train)) / sum(y_train)
)

# Train final model with optimized parameters
final_model <- lgb.train(
  params = best_params,
  data = lgb_train,
  nrounds = 1000,
  verbose = 0
)

# Make predictions on full dataset
pred_probs <- predict(final_model, X_full_matrix)

# Add predictions to the original dataset
pbp_2022_2024$xSwing <- pred_probs

# Optional: ROC Curve for Training Data
# First, create a test set from the training data
set.seed(123)
sample_split <- sample.split(Y = y_train, SplitRatio = 0.75)
train_subset <- subset(one_hot_pbp_2022_2024, sample_split == TRUE)
test_subset <- subset(one_hot_pbp_2022_2024, sample_split == FALSE)

X_test_subset <- test_subset %>% select(-swing)
y_test_subset <- test_subset$swing
X_test_matrix_subset <- as.matrix(X_test_subset)

# Predictions on test subset
test_pred_probs <- predict(final_model, as.matrix(X_test_subset))

# ROC Curve
roc_curve <- roc(y_test_subset, test_pred_probs)
auc_value <- auc(roc_curve)

# Plot ROC Curve
roc_data <- data.frame(
  tpr = rev(roc_curve$sensitivities),
  fpr = rev(1 - roc_curve$specificities)
)

ggplot(roc_data, aes(x = fpr, y = tpr)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  labs(
    title = "ROC Curve for Swing Prediction Model",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)"
  ) +
  annotate("text", x = 0.7, y = 0.2, 
           label = paste("AUC =", round(auc_value, 4)), 
           color = "red") +
  theme_minimal()

# Feature Importance Plot
importance_matrix <- lgb.importance(final_model, percentage = TRUE)
importance_matrix <- importance_matrix %>%
  rename(Feature = Feature, Importance = Gain)

ggplot(data = importance_matrix, aes(x = reorder(Feature, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "blue") +
  coord_flip() +
  labs(
    title = "Feature Importance",
    x = "Features",
    y = "Importance (Gain)"
  ) +
  theme_minimal()

#===============================================================================
# Find most likely expected event
#===============================================================================

pbp_2022_2024$xContact <- 1 - pbp_2022_2024$xWhiff
pbp_2022_2024$xFoul <- 1 - pbp_2022_2024$xIn_Play

# Ensure pbp_2022_2024 is a data frame
pbp_2022_2024 <- as.data.frame(pbp_2022_2024)

# Calculate the probabilities at each leaf of the tree
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    # Takes
    pct_CalledStrike = (1 - xSwing) * xStrike * 100,
    pct_Ball = (1 - xSwing) * xBall * 100,
    pct_HBP = (1 - xSwing) * xHBP * 100,
    
    # Swings
    pct_Whiff = xSwing * xWhiff * 100,
    pct_Contact_Foul = xSwing * xContact * xFoul * 100,
    pct_Contact_in_play_GB = xSwing * xContact * xIn_Play * xGB * 100,
    pct_Contact_in_play_LD = xSwing * xContact * xIn_Play * xLD * 100,
    pct_Contact_in_play_FB = xSwing * xContact * xIn_Play * xFB * 100,
    pct_Contact_in_play_PU = xSwing * xContact * xIn_Play * xPU * 100
  )

# Verify the percentages sum to 100 for each row
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    total_pct = pct_CalledStrike + pct_Ball + pct_HBP +
      pct_Whiff + pct_Contact_Foul + 
      pct_Contact_in_play_GB + pct_Contact_in_play_LD + 
      pct_Contact_in_play_FB + pct_Contact_in_play_PU
  )


# Normalize probabilities for hit types so they sum to 1
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    total = xGB + xLD + xFB + xPU,
    xGB = xGB / total,
    xLD = xLD / total,
    xFB = xFB / total,
    xPU = xPU / total
  )

# Identify the most likely event
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    expected_event = case_when(
      pct_Contact_in_play_GB == pmax(pct_Contact_in_play_GB, pct_Contact_in_play_LD, pct_Contact_in_play_FB, pct_Contact_in_play_PU, pct_CalledStrike, pct_Ball, pct_HBP, pct_Whiff, pct_Contact_Foul) ~ "GroundBall",
      pct_Contact_in_play_LD == pmax(pct_Contact_in_play_GB, pct_Contact_in_play_LD, pct_Contact_in_play_FB, pct_Contact_in_play_PU, pct_CalledStrike, pct_Ball, pct_HBP, pct_Whiff, pct_Contact_Foul) ~ "LineDrive",
      pct_Contact_in_play_FB == pmax(pct_Contact_in_play_GB, pct_Contact_in_play_LD, pct_Contact_in_play_FB, pct_Contact_in_play_PU, pct_CalledStrike, pct_Ball, pct_HBP, pct_Whiff, pct_Contact_Foul) ~ "FlyBall",
      pct_Contact_in_play_PU == pmax(pct_Contact_in_play_GB, pct_Contact_in_play_LD, pct_Contact_in_play_FB, pct_Contact_in_play_PU, pct_CalledStrike, pct_Ball, pct_HBP, pct_Whiff, pct_Contact_Foul) ~ "PopUp",
      pct_CalledStrike == pmax(pct_Contact_in_play_GB, pct_Contact_in_play_LD, pct_Contact_in_play_FB, pct_Contact_in_play_PU, pct_CalledStrike, pct_Ball, pct_HBP, pct_Whiff, pct_Contact_Foul) ~ "CalledStrike",
      pct_Ball == pmax(pct_Contact_in_play_GB, pct_Contact_in_play_LD, pct_Contact_in_play_FB, pct_Contact_in_play_PU, pct_CalledStrike, pct_Ball, pct_HBP, pct_Whiff, pct_Contact_Foul) ~ "CalledBall",
      pct_HBP == pmax(pct_Contact_in_play_GB, pct_Contact_in_play_LD, pct_Contact_in_play_FB, pct_Contact_in_play_PU, pct_CalledStrike, pct_Ball, pct_HBP, pct_Whiff, pct_Contact_Foul) ~ "hit_by_pitch",
      pct_Whiff == pmax(pct_Contact_in_play_GB, pct_Contact_in_play_LD, pct_Contact_in_play_FB, pct_Contact_in_play_PU, pct_CalledStrike, pct_Ball, pct_HBP, pct_Whiff, pct_Contact_Foul) ~ "Whiff",
      pct_Contact_Foul == pmax(pct_Contact_in_play_GB, pct_Contact_in_play_LD, pct_Contact_in_play_FB, pct_Contact_in_play_PU, pct_CalledStrike, pct_Ball, pct_HBP, pct_Whiff, pct_Contact_Foul) ~ "Foul",
      TRUE ~ NA_character_
    )
  )


# Identify the most likely event
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    actual_event = case_when(
      description %in% c("ball", "blocked_ball") ~ "CalledBall",
      description == "called_strike" ~ "CalledStrike",
      description == "foul" ~ "Foul",
      description == "swinging_strike" ~ "Whiff",
      description == "hit_by_pitch" ~ "hit_by_pitch",
      batted_ball_type == "GroundBall" ~ "GroundBall",
      batted_ball_type == "LineDrive" ~ "LineDrive",
      batted_ball_type == "FlyBall" ~ "FlyBall",
      batted_ball_type == "PopUp" ~ "PopUp",
      TRUE ~ NA_character_
    )
  )

# Calculate accuracy
accuracy <- pbp_2022_2024 %>%
  summarise(
    total_events = n(),
    correct_predictions = sum(actual_event == expected_event, na.rm = TRUE),
    accuracy = correct_predictions / total_events
  )

# View the result
accuracy


#===============================================================================
# Get run values for events
#===============================================================================

# Define the run values for each count
run_values <- c(
  "0-0" = 0.001695997, "1-0" = 0.039248042, "0-1" = -0.043581338,
  "2-0" = -0.043581338, "1-1" = -0.015277684, "0-2" = -0.103242476,
  "3-0" = 0.200960731, "2-1" = 0.034545018, "1-2" = -0.080485991,
  "3-1" = 0.138254876, "2-2" = -0.039716495, "3-2" = 0.048505049,
  "Walk" = 0.325, "Strikeout" = -0.284,
  "FlyBall" = 0.586, "LineDrive" = 0.528,
  "GroundBall" = 0.164, "PopUp" = 0.0186
)

# Calculate expected run value change
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    # Current count as a string
    current_count = paste(balls, strikes, sep = "-"),
    current_run_value = run_values[current_count],
    
    # New count based on expected event
    new_count = case_when(
      expected_event == "CalledBall" & balls < 3 ~ paste(balls + 1, strikes, sep = "-"),
      expected_event %in% c("CalledStrike", "Whiff", "Foul") & strikes < 2 ~ paste(balls, strikes + 1, sep = "-"),
      expected_event == "Walk" ~ "Walk",
      expected_event == "Strikeout" ~ "Strikeout",
      expected_event %in% c("FlyBall", "LineDrive", "GroundBall", "PopUp") ~ expected_event,
      TRUE ~ current_count # Default to current count if no change
    ),
    new_run_value = run_values[new_count],
    
    # Calculate the run value change
    run_value_change = new_run_value - current_run_value
  )

# Calculate actual run value change
pbp_2022_2024 <- pbp_2022_2024 %>%
  mutate(
    # Current count as a string
    current_count = paste(balls, strikes, sep = "-"),
    current_run_value = run_values[current_count],
    
    # New count based on actual event
    new_count_actual = case_when(
      actual_event == "CalledBall" & balls < 3 ~ paste(balls + 1, strikes, sep = "-"),
      actual_event %in% c("CalledStrike", "Whiff", "Foul") & strikes < 2 ~ paste(balls, strikes + 1, sep = "-"),
      actual_event == "Walk" ~ "Walk",
      actual_event == "Strikeout" ~ "Strikeout",
      actual_event %in% c("FlyBall", "LineDrive", "GroundBall", "PopUp") ~ actual_event,
      TRUE ~ current_count # Default to current count if no change
    ),
    new_run_value_actual = run_values[new_count_actual],
    
    # Calculate the run value change
    run_value_change_actual = new_run_value_actual - current_run_value
  )

pbp_2022_2024$swing_decision_rv <- pbp_2022_2024$run_value_change_actual - pbp_2022_2024$run_value_change

pbp_2022_2024 <- pbp_2022_2024 %>%
  select(-run_value_change_actual, -new_run_value_actual, -new_count_actual, -run_value_change, -new_run_value, -new_count, -current_run_value, -current_count, -total, - total_pct, -xFoul, -xContact, -xSwing, -swing, -xPU, -xFB, -xLD, -xGB, -xIn_Play, -xWhiff, -xHBP, -xStrike, -xBall)

rm(accuracy)
rm(best_params)
rm(bounds)
rm(evaluation_summary)
rm(final_model)
rm(HAA_model)
rm(importance_matrix)
rm(lgb_train)
rm(model_data)
rm(one_hot_pbp_2022_2024)
rm(one_hot_pbp_2022_2024_full)
rm(opt_results)
rm(params)
rm(predictions)
rm(results)
rm(roc_data)
rm(roc_curve)
rm(pitch_type_onehot_matrix)
rm(pitch_type_onehot_matrix_full)
rm(pitch_type_onehot_matrix_prediction)
rm(pitch_type_onehot_matrix_train)
rm(test_set)
rm(test_subset)
rm(train_set)
rm(train_subset)
rm(X_full_matrix)
rm(X_test)
rm(X_test_matrix_subset)
rm(X_test_subset)
rm(X_train)
rm(X_train_matrix)
rm(auc_value)
rm(col)
rm(contact_types)
rm(dtest)
rm(dtrain)
rm(missing_cols)
rm(missing_columns)
rm(pitch_call)
rm(pred_probs)
rm(relevant_features)
rm(run_values)
rm(sample_split)
rm(swing_types)
rm(target)
rm(target_columns)
rm(test_pred_probs)
rm(train_indices)
rm(training_columns)
rm(y_test)
rm(y_test_subset)
rm(y_train)
rm(lgb_opt_function)
rm(train_model)
rm(whiff_objective)

