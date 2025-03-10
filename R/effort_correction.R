#' A helper function to generate effort (detectability) corrected counts for each checklist
#'
#' This function runs XGboost to predict count for
#'    original conditions (checlist duration, distance, hour of day, number of observer, mode of checklist,
#'    cci), as well for standarndized condition (medians). Then, it estimates the 'correction factor, cf'
#'    which is used to standardize (correct for effort and detection) counts.
#'
#' @param sp_data zero-filled 'sp_data' data.table object, output from \code{import_from_erd} function
#' @return A \code{data.table} object with zero-filled abundance (effort corrected) and presence-only reported columns
#'    as well as additional columns for latitude, longitude, year, month, day_of_year, hours_of_day,
#'    protocol_id, is_stationary, is_traveling, effort_hrs, effort_distance_km, cci
#' @export
effort_correction <-
  function(sp_data) {

    # take variables of interest and then filter only occurence data
    sp_data1 <- sp_data %>% select(obs_count, day_of_year, effort_distance_km, # should exclude latlong
                                   hours_of_day, num_observers, effort_hrs, is_stationary,
                                   cci) %>%
      mutate(day_of_year = (day_of_year - min(day_of_year)) / (max(day_of_year) - min(day_of_year)))  # re-scaling day_of_year as 0-1
    # %>%
 #     filter(obs_count > 0)

   # track indices
   sp_data1_indices <- which(!is.na(sp_data$obs_count))

    # Categorical values need to be converted to numeric
    sp_data1$is_stationary  <- as.integer(sp_data1$is_stationary)

    ### Each cell-season should have at least 10 checklists!
    if(nrow(sp_data1) < 10) {
      message("less than 10 data points detected. No model for this cell and season!")
      return(data.frame()) # returning an empty df to keep consistent str for concatenating the list of df later on!
    }

    ### check for NA values (XGBoost does not take NA values)
    if(sum(is.na(sp_data1)) > 0) {
      message("NA exists, XGboost does not take NA values!")
      return(data.frame()) # returning an empty df to keep consistent str for concatenating the list of df later on!
    }

    ### data preparation for XGboost

    # Convert to XGBoost D.matrix and then divide into training set (to train data), validation set (to tune hyperparameters),
    # and test set (model evaluation).

    ### Set seed
    set.seed(123)

    ### randomize the dataset
    df_shuffled <- sample(nrow(sp_data1))
    sp_data1 <- sp_data1[df_shuffled, ]

    ### Get the number of each set (in case its not an even break like 1000)
    n_train <- round(0.6 * nrow(sp_data1))
    n_val <- round(0.2 * nrow(sp_data1))
    n_test<- nrow(sp_data1) - n_train - n_val

    ### Create indices for the training, validation, and test sets (so you know which rows you sampled from)
    train_indices <- df_shuffled [1:n_train]
    valid_indices <- df_shuffled [(n_train + 1):(n_train + n_val)]
    test_indices  <- df_shuffled [(n_train + n_val + 1):nrow(sp_data1)]

    ### Divide the dataset into groups:
    train_set <- sp_data1[train_indices, ]
    valid_set <- sp_data1[valid_indices, ]
    test_set  <- sp_data1[test_indices, ]

    ### Prep training set and convert to data matrix
    # remove label variable from the training data
    train_x = data.matrix(train_set[, -1]) #removing the column we wish to predict
    train_y = data.matrix(train_set[,1]) ### this is the label

    ### Prep validation set
    valid_x = data.matrix(valid_set[, -1])
    valid_y = data.matrix(valid_set[, 1])

    ### Prep test set
    test_x = data.matrix(test_set[, -1])
    test_y = data.matrix(test_set[, 1])

    # Ensure train_y is a numeric vector
    train_y <- as.vector(train_y)
    valid_y <- as.vector(valid_y)

    # Convert training and validation sets to DMatrix
    xgb_train <- xgb.DMatrix(data = train_x, label = train_y)
    xgb_val <- xgb.DMatrix(data = valid_x, label = valid_y)

    # # should you tune hyperparameters? Define a grid of hyperparameters
    # tune_grid <- expand.grid(
    #   nrounds = seq(100, 1000, by = 100),            # Number of boosting iterations
    #   max_depth = seq(3, 10, by = 2),                # Maximum depth of a tree
    #   eta = c(0.01, 0.05, 0.1, 0.3),                 # Learning rate
    #   gamma = c(0, 0.1, 0.2, 0.5),                   # Minimum loss reduction required to make a split
    #   colsample_bytree = c(0.5, 0.7, 1),             # Subsample ratio of columns when constructing each tree
    #   min_child_weight = c(1, 3, 5),                 # Minimum sum of instance weight needed in a child
    #   subsample = c(0.5, 0.7, 1)                     # Subsample ratio of the training instances
    # )
    #
    # # run hyperparameter tuning in parallel
    # cores <- 7
    # cl <- makeCluster(cores)
    # registerDoParallel(cl)
    #
    # # Set up training control
    # train_control <- trainControl(
    #   method = "cv",                   # cross-validation
    #   number = 3,                      # number of folds
    #   allowParallel = TRUE,
    #   verboseIter = TRUE             # print training log
    # )
    #
    # tictoc::tic()
    # # Train the model using caret
    # xgb_train_caret <- train(
    #   x = train_x,
    #   y = train_y,
    #   method = "xgbTree",
    #   trControl = train_control,
    #   tuneGrid = tune_grid,
    #   metric = "RMSE"
    # ). # this proccess took 3.5 hours!
    # tictoc::toc()
    # stopCluster(cl)
    # registerDoSEQ()  # Revert to sequential processing
    #
    # # Extract the best hyperparameters
    # best_params <- xgb_train_caret$bestTune

    # best hyperparameters from above tuning!
    # best_params <- structure(
    #   list(
    #     nrounds = 400,
    #     max_depth = 7,
    #     eta = 0.01,
    #     gamma = 0.5,
    #     colsample_bytree = 0.7,
    #     min_child_weight = 5,
    #     subsample = 0.7
    #   ),
    #   row.names = 3134L,
    #   class = "data.frame"
    # )
    #
    # model <- xgb.train(
    #   params = list(
    #     objective = "reg:tweedie",
    #     eval_metric = "poisson-nloglik",
    #     booster = "gbtree",
    #     eta = best_params$eta,
    #     max_depth = best_params$max_depth,
    #     gamma = best_params$gamma,
    #     colsample_bytree = best_params$colsample_bytree,
    #     min_child_weight = best_params$min_child_weight,
    #     subsample = best_params$subsample,
    #     nthread = 8
    #   ),
    #   data = xgb_train,
    #   nrounds = best_params$nrounds,
    #   watchlist = list(train = xgb_train, eval = xgb_val),
    #   early_stopping_rounds = 10,
    #   verbose = 1
    # )


    # Train the final model with the default parameters
    model <- xgb.train(
      params = list(
        objective = "reg:tweedie",
        eval_metric = "poisson-nloglik",
        booster = "gbtree",
        eta = 0.3,                                  # default: 0.3
        max_depth = 6,
        nthread = 8
        ),
      data = xgb_train,
      nrounds = 500,
      watchlist = list(train = xgb_train, eval = xgb_val),
      early_stopping_rounds = 10,
      verbose = 1
    )

    # Evaluate the model performance on training and test sets
    train_set$pred <- predict(model, xgb.DMatrix(data = train_x))
    train_rmse <- sqrt(mean((train_y - train_set$pred)^2, na.rm = TRUE))

    test_set$pred <- predict(model, xgb.DMatrix(data = test_x))
    test_rmse <- sqrt(mean((test_y - test_set$pred)^2, na.rm = TRUE))

    # Print diagnostic information
    cat("Train RMSE:", train_rmse, "\n")
    cat("Test RMSE:", test_rmse, "\n")
    cat("Train SD:", sd(train_set$obs_count, na.rm = TRUE), "\n")
    cat("Test SD:", sd(test_set$obs_count, na.rm = TRUE), "\n")

    # # Check for NA values in RMSE and standard deviations
    # if (is.na(train_rmse) | is.na(test_rmse) | is.na(sd(train_set$obs_count, na.rm = TRUE)) | is.na(sd(test_set$obs_count, na.rm = TRUE))) {
    #   message("NA values detected in RMSE or standard deviation calculations. Please check your data.")
    #   return(data.frame())
    # }
    #
    # # Run predictions if RMSE are lower than SD
    # if (!(train_rmse < sd(train_set$obs_count, na.rm = TRUE)) & !(test_rmse < sd(test_set$obs_count, na.rm = TRUE))) {
    #   message("rmse is larger than sd !")
    #   return(data.frame())
    # }


    all_data_x <- rbind(train_x, valid_x, test_x)  # re-combine training, validation test data
    pred_org <- predict(model, xgb.DMatrix(data = all_data_x))

    sim_data_x <- all_data_x %>% as.data.frame() %>%
      mutate(effort_distance_km = 1,
             hours_of_day = 8, #  8 am
             day_of_year = median(day_of_year),
             num_observers = 1,
             effort_hrs = 1,
             is_stationary = 0,       # travelling checklist
             cci = 3 # 75% ?
      ) %>% as.matrix()

    pred_sim <- predict(model, xgb.DMatrix(data = sim_data_x))


    # calculation of CF
    cf <- pred_sim / pred_org

    # re-order the predictions to match the original order of sp_data1
    cf <- cf[order(df_shuffled)]
    # pred_org <- pred_org[order(df_shuffled)]
    # pred_sim <- pred_sim[order(df_shuffled)]

    # corrected count
    sp_data1 <- sp_data1[order(df_shuffled), ] %>% mutate(corr_count = obs_count * cf, cf = cf)

    # merge back to the orginal sp_data, creating sp_data2
    sp_data2 <- sp_data
    sp_data2$corr_count <- 0 # initialize corr_count column with 0, store both!
    sp_data2$corr_count[sp_data1_indices] <- sp_data1$corr_count
    sp_data2$correction_factor <- NA
    sp_data2$correction_factor[sp_data1_indices] <- sp_data1$cf
    return(sp_data2)

  }



#' A helper function to generate effort (detectability) corrected counts for each checklist
#'
#' This function runs XGboost to predict count for
#'    original conditions (checlist duration, distance, hour of day, number of observer, mode of checklist,
#'    cci), as well for standarndized condition (medians). Then, it estimates the 'correction factor, cf'
#'    which is used to standardize (correct for effort and detection) counts.
#'
#' @param sp_data zero-filled 'sp_data' data.table object, output from \code{import_from_erd} function
#' @return A \code{data.table} object with zero-filled abundance (effort corrected) and presence-only reported columns
#'    as well as additional columns for latitude, longitude, year, month, day_of_year, hours_of_day,
#'    protocol_id, is_stationary, is_traveling, effort_hrs, effort_distance_km, cci
#' @export
effort_residual_correction <-
  function(sp_data) {

    message("Original row count: ", nrow(sp_data))

    # take variables of interest and then filter only occurence data
    sp_data1 <- sp_data %>% select(obs_count, day_of_year, effort_distance_km, # should exclude latlong
                                   hours_of_day, num_observers, effort_hrs, is_stationary,
                                   cci) %>%
      mutate(day_of_year = (day_of_year - min(day_of_year)) / (max(day_of_year) - min(day_of_year)))  # re-scaling day_of_year as 0-1


    # track indices
    sp_data1_indices <- which(!is.na(sp_data1$obs_count))

    # Categorical values need to be converted to numeric
    sp_data1$is_stationary  <- as.integer(sp_data1$is_stationary)

    ### Each cell-season should have at least 10 checklists!
    if(nrow(sp_data1) < 10) {
      message("less than 10 data points detected. No model for this cell and season!")
      return(data.frame()) # returning an empty df to keep consistent str for concatenating the list of df later on!
    }

    ### check for NA values (XGBoost does not take NA values)
    if(sum(is.na(sp_data1)) > 0) {
      message("NA exists, XGboost does not take NA values!")
      return(data.frame()) # returning an empty df to keep consistent str for concatenating the list of df later on!
    }

    ### data preparation for XGboost

    # Convert to XGBoost D.matrix and then divide into training set (to train data), validation set (to tune hyperparameters),
    # and test set (model evaluation).

    ### Set seed
    set.seed(123)

    ### randomize the dataset
    df_shuffled <- sample(nrow(sp_data1))
    sp_data1 <- sp_data1[df_shuffled, ]

    message("sp_data1 row count after shuffling: ", nrow(sp_data1))

    ### Get the number of each set (in case its not an even break like 1000)
    n_train <- round(0.6 * nrow(sp_data1))
    n_val <- round(0.2 * nrow(sp_data1))
    n_test<- nrow(sp_data1) - n_train - n_val

    ### Create indices for the training, validation, and test sets (so you know which rows you sampled from)
    train_indices <- df_shuffled [1:n_train]
    valid_indices <- df_shuffled [(n_train + 1):(n_train + n_val)]
    test_indices <- df_shuffled [(n_train + n_val + 1):nrow(sp_data1)]

    ### Divide the dataset into groups:
    train_set <- sp_data1[train_indices, ]
    message("Train set row count: ", nrow(train_set))

    valid_set <- sp_data1[valid_indices, ]
    message("valid set row count: ", nrow(valid_set))

    test_set  <- sp_data1[test_indices, ]
    message("Test set row count: ", nrow(test_set))

    ### Prep training set and convert to data matrix
    # remove label variable from the training data
    train_x = data.matrix(train_set[, -1]) #removing the column we wish to predict
    train_y = data.matrix(train_set[,1]) ### this is the label

    ### Prep validation set
    valid_x = data.matrix(valid_set[, -1])
    valid_y = data.matrix(valid_set[, 1])

    ### Prep test set
    test_x = data.matrix(test_set[, -1])
    test_y = data.matrix(test_set[, 1])

    if (nrow(train_x) + nrow(valid_x) + nrow(test_x) != nrow(sp_data1)) {
      message("Warning: The number of rows in train, validation, and test sets does not add up!")
    }

    # Ensure train_y is a numeric vector
    train_y <- as.vector(train_y)
    valid_y <- as.vector(valid_y)

    # Convert training and validation sets to DMatrix
    xgb_train <- xgb.DMatrix(data = train_x, label = train_y)
    xgb_val <- xgb.DMatrix(data = valid_x, label = valid_y)

    # Train the final model with the default parameters
    model <- xgb.train(
      params = list(
        objective = "reg:tweedie",
        eval_metric = "poisson-nloglik",
        booster = "gbtree",
        eta = 0.3,                                  # default: 0.3
        max_depth = 6,
        nthread = 8
      ),
      data = xgb_train,
      nrounds = 500,
      watchlist = list(train = xgb_train, eval = xgb_val),
      early_stopping_rounds = 10,
      verbose = 1
    )

    # Evaluate the model performance on training and test sets
    train_set$pred <- predict(model, xgb.DMatrix(data = train_x))
    train_rmse <- sqrt(mean((train_y - train_set$pred)^2, na.rm = TRUE))

    test_set$pred <- predict(model, xgb.DMatrix(data = test_x))
    test_rmse <- sqrt(mean((test_y - test_set$pred)^2, na.rm = TRUE))

    # Print diagnostic information
    cat("Train RMSE:", train_rmse, "\n")
    cat("Test RMSE:", test_rmse, "\n")
    cat("Train SD:", sd(train_set$obs_count, na.rm = TRUE), "\n")
    cat("Test SD:", sd(test_set$obs_count, na.rm = TRUE), "\n")

    # Check for NA values in RMSE and standard deviations
    if (is.na(train_rmse) | is.na(test_rmse) | is.na(sd(train_set$obs_count)) | is.na(sd(test_set$obs_count))) {
      message("NA values detected in RMSE or standard deviation calculations. Please check your data.")
      return(data.frame())
    }

    # Run predictions if RMSE are lower than SD
    if (!(train_rmse < sd(train_set$obs_count, na.rm = TRUE)) & !(test_rmse < sd(test_set$obs_count, na.rm = TRUE))) {
      message("rmse is larger than sd !")
      return(data.frame())
    }


    all_data_x <- rbind(train_x, valid_x, test_x)  # re-combine training, validation test test data
    message("all_data_x row count: ", nrow(all_data_x))

    pred_org <- predict(model, xgb.DMatrix(data = all_data_x))

    sim_data_x <- all_data_x %>% as.data.frame() %>%
      mutate(effort_distance_km = 1,
             hours_of_day = 8, #  8 am
             day_of_year = median(day_of_year),
             num_observers = 1,
             effort_hrs = 1,
             is_stationary = 0,       # travelling checklist
             cci = 3 # 75% ?
      ) %>% as.matrix()

    message("sim_data_x set row count: ", nrow(sim_data_x))

    pred_sim <- predict(model, xgb.DMatrix(data = sim_data_x))

    # re-order the predictions to match the original order of sp_data1
    pred_org <- pred_org[order(df_shuffled)]
    pred_sim <- pred_sim[order(df_shuffled)]

    # corrected count
    sp_data1 <- sp_data1[order(df_shuffled), ] %>% mutate(corr_count = pred_sim - obs_count + pred_org)

    # merge back to the orginal sp_data, creating sp_data2
    sp_data2 <- sp_data
    sp_data2$corr_count <- 0 # initialize corr_count column with 0, store both!
    sp_data2$corr_count[sp_data1_indices] <- sp_data1$corr_count

    message("Final sp_data2 row count: ", nrow(sp_data2))

    return(sp_data2)


  }
