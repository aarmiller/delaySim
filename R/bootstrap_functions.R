
#' Run multiple simulations of missed visits
#'
#' @param sim_data a dataset containing the time_map, miss bins and other parameters used for the simulation.
#'                 This dataset should be created using the `prep_sim_data()` function
#' @param trials number of simulationed trials to run (default is 50)
#' @param no_bootstrapping Specifies whether you want to run the simulations without bootstrapping the original dataset
#' @param num_cores The number of worker cores to use. If not specified will determined the number of cores based on the which ever
#' is the smallest value between number of trials or detected number of cores - 1
#' @param sim_algorithm a string defining the algorithm to use. Option include "simple", 
#'                      "simple_correlated" and "general"
#' @param sim_ctrl the sim_ctrl function that defines additional parameter values passed to the 
#'                 simulation. If no value is supplied (i.e., sim_ctrl = NULL) then default 
#'                 values for the sim_ctrl() function will be used
#' @export
#'
run_sim_miss_visits <- function (sim_data, trials = 50, no_bootstrapping = FALSE, sim_algorithm="simple", sim_ctrl = NULL,
                                 num_cores = NULL) {
  
  if (no_bootstrapping == FALSE){
    tmp <- tibble::tibble(trial = 1:trials) %>%
      dplyr::mutate(data = purrr::map(trial,
                                    ~sim_miss_visits(sim_data = sim_data, sim_algorithm = sim_algorithm,
                                                     sim_ctrl = sim_ctrl)))
  } else {
    simulation_data <- sim_data

    if (is.null(num_cores)) {
      num_cores <- min(trials, parallel::detectCores() - 1)
    } else {
      num_cores <- num_cores
    }

    cluster <- parallel::makeCluster(num_cores)

    parallel::clusterCall(cluster, function() library(tidyverse))
    parallel::clusterCall(cluster, function() library(delayDX))

    test <- parallel::parLapply(cl = cluster,
                               1:trials,
                               function(x){sim_miss_visits(sim_data = simulation_data,  sim_algorithm = sim_algorithm,
                                                           sim_ctrl = sim_ctrl)})
    parallel::stopCluster(cluster)
    gc()

    tmp <- tibble()
    for (i in 1:length(test)){
      tmp1 <- tibble(trial = i,
                     data = test[i])
      tmp <- bind_rows(tmp, tmp1)

    }
  }
  return(tmp)
}

#' Bootstrap Estimation of Changepoint and simulated miss visits
#'
#' @param sim_data a dataset containing the time_map, miss bins and other parameters used for the simulation.
#'                 This dataset should be created using the `prep_sim_data()` function
#' @param n_sim_trials Number of trials to run in simulation of miss visits or miss patients
#' @param sim_algorithm a string defining the algorithm to use. Option include "simple", 
#'                      "simple_correlated" and "general"
#' @param sim_ctrl the sim_ctrl function that defines additional parameter values passed to the 
#'                 simulation. If no value is supplied (i.e., sim_ctrl = NULL) then default 
#'                 values for the sim_ctrl() function will be used
#' @param no_bootstrapping Specifies whether you want to run the simulations without bootstrapping the original dataset
#' @param num_cores The number of worker cores to use. If not specified will detect cores and use 1 less than the number of cores
#' @export
#'
boot_change_point <- function (sim_data, n_sim_trials = 100L,
                               sim_algorithm="simple", sim_ctrl = NULL,
                               eval_criteria="AIC", week_period=FALSE, num_cores = NULL,
                               no_bootstrapping = FALSE, auto_reg = FALSE) {
  if (no_bootstrapping == FALSE){
    # draw bootstrapped samples
    draw_time_map <- sim_data$time_map %>%
      dplyr::distinct(enrolid) %>%
      dplyr::sample_frac(1, replace = TRUE) %>%                 # sample enrolids with replacement
      dplyr::mutate(enrolid_new = row_number()) %>%             # generate new unique enrolids
      dplyr::inner_join(select(sim_data$time_map,-enrolid_new), # remove old enrolid_new
                        by = "enrolid") %>%                     # merge back into time map
      dplyr::mutate(enrolid_old=enrolid,                        # old enrolid for tracking
                    enrolid=enrolid_new,
                    patient_id=enrolid_new)                        # change to new enrolid, change name to patient id  (NOTE NEEDS TO OCCUR TO COUNT PATIENTS CORRECTLY)
  } else {
    draw_time_map <- sim_data$time_map
  }
  # simulation change-point data
  sim_cp <- prep_cp_data(draw_time_map,
                         event_name = "miss_ind",   # note the count value has already been selected in the sim data
                         start_day = 0L,  # to account for adjusted days_since_dx value
                         by_days = 1L)

    if (sim_data$cp_method != "set_cp"){
      sim_cp <- sim_cp %>%
        find_change_point(var_name = "n_miss_visits",
                            method = sim_data$cp_method,
                            eval_criteria = sim_data$eval_criteria,
                            week_period = sim_data$week_period,
                            auto_reg = sim_data$auto_reg,
                            specify_cp = sim_data$specify_cp)
    } else {
      sim_cp <- sim_cp %>%
        set_change_point(var_name = "n_miss_visits",
                         method = sim_data$set_cp_method,
                         return_miss_only = FALSE,
                         compare_all_methods = FALSE,
                         week_period = sim_data$week_period,
                         specify_cp = sim_data$specify_cp)
    }
      # pull out the miss bins
      miss_bins <- sim_cp$miss_bins

      # update miss bins to reflect predicted value or upper bound prediction interval
      if (sim_data$prediction_bound_for_sim == FALSE){

        miss_bins <- miss_bins %>% select(period, Y,
                                          pred1, pred,
                                          num_miss)

      } else {
        miss_bins <- miss_bins %>% select(period, Y,
                                           pred1 = upper_int_pred1,
                                           pred,
                                           num_miss = num_miss_upper_int)
      }

      # estimated number of missed visits and observed missed visits
      miss_stats <- miss_bins %>%
        dplyr::summarise(miss_visits_est = sum(num_miss),
                         miss_visits_obs = sum(Y - pred1))

      # update sim data, prepare for run_sim_miss_visits
      new_sim_data <- sim_data 
      new_sim_data$time_map <- draw_time_map
      new_sim_data$miss_bins_visits <- miss_bins
      new_sim_data$change_point <- sim_cp$change_point$period

      # run simulation
      sim_res_results <- run_sim_miss_visits(sim_data = new_sim_data,
                                              trials = n_sim_trials,
                                              sim_algorithm = sim_algorithm,
                                              sim_ctrl = sim_ctrl,
                                              no_bootstrapping = no_bootstrapping,
                                              num_cores = num_cores)

      
      # aggregate results
      results <- list(change_point = sim_cp$change_point,
                      pred = sim_cp$pred,
                      miss_counts = miss_stats,
                      sim_visit_results = sim_res_results)
      return(results)
    }

  


#' Run multiple bootstrapped change_point simulations (in parallel)
#'
#' @param sim_data A dataset containing the time_map, miss bins and other parameters used for the simulation.
#'                 This dataset should be created using the `prep_sim_data()` function
#' @param new_draw_weight A weighing parameter used to assign preference to drawing previously "missed" patients
#'                        at each time step. A value of 0 applies strict preference to drawing patients who
#'                        have been assigned to miss in prior time steps, while a value 0.5 applies equal weight
#'                        to patients who have and have not been previously selected
#' @param boot_trials The number of bootstrapped trials to run (default is 100)
#' @param n_sim_trials The number of trials to run in simulation of miss visits (default is 50)
#' @param num_cores The number of worker cores to use. If not specified will determined the number of cores based on the which ever
#' is the smallest value between number of boot_trials or detected number of cores - 1
#' @param no_bootstrapping Specifies whether you want to run the simulations without bootstrapping the original dataset
#' 
#' @examples
#'
#' ### Run simulations with bootstrapping and allow change point to vary with each bootstrap sample ###
#'
#' #load example final_time_map dataset
#' load("/Shared/Statepi_Diagnosis/grant_projects/hsv_enceph/scripts/validation/enrolled_ge_365/report_data.RData")
#'
#' # rename ED column
#' final_time_map <- final_time_map %>% rename(ed = ED)
#'
#' # run prep sim function
#' tmp_sim_data <- prep_sim_data(final_time_map, event_name = "any_ssd", cp_method = "cusum", start_day = 1L, by_days = 1L,
#'                               week_period = TRUE)
#'
#' #run simulations on number of visits
#' simulation_results <- run_cp_bootstrap(tmp_sim_data,
#'                                        boot_trials = 500,
#'                                        n_sim_trials = 50,
#'                                        sim_algorithm = sim_algorithm,
#'                                        sim_ctrl = sim_ctrl,
#'                                        num_cores = NULL,
#'                                        no_bootstrapping = FALSE)
#'
#'
#' ### Run simulations with bootstrapping and specify a change point applied to each bootstrap sample ###
#'
#' # set a change point
#' cp <- 20L
#'
#' # run prep sim function
#' tmp_sim_data <- prep_sim_data(final_time_map, event_name = "any_ssd", cp_method = "cusum", start_day = 1L, by_days = 1L,
#'                               week_period = TRUE, specify_cp = cp)
#'
#' #run simulations on number of visits
#' simulation_results <- run_cp_bootstrap(tmp_sim_data,
#'                                        boot_trials = 500,
#'                                        n_sim_trials = 50,
#'                                        sim_algorithm = sim_algorithm,
#'                                        sim_ctrl = sim_ctrl,
#'                                        num_cores = NULL,
#'                                        no_bootstrapping = FALSE)
#'
#'
#' ### Run simulations without bootstrapping and allow function to find the optimal change point for inputted data ###
#'
#' # run prep sim function
#' tmp_sim_data <- prep_sim_data(final_time_map, event_name = "any_ssd", cp_method = "cusum", start_day = 1L, by_days = 1L,
#'                               week_period = TRUE)
#'
#' #run simulations on number of visits
#' simulation_results <- run_cp_bootstrap(tmp_sim_data,
#'                                        boot_trials = 0,
#'                                        n_sim_trials = 25000,
#'                                        sim_algorithm = sim_algorithm,
#'                                        sim_ctrl = sim_ctrl,
#'                                        num_cores = NULL,
#'                                        no_bootstrapping = TRUE)
#'
#'
#' ### Run simulations without bootstrapping and specify a change point instead of allowing function to find the optimal change point to apply to the data ###
#'
#' # set a change point
#' cp <- 20L
#'
#' # run prep sim function
#' tmp_sim_data <- prep_sim_data(final_time_map, event_name = "any_ssd", cp_method = "cusum", start_day = 1L, by_days = 1L,
#'                               week_period = TRUE, specify_cp = cp)
#'
#' #run simulations on number of visits
#' simulation_results <- run_cp_bootstrap(tmp_sim_data,
#'                                        boot_trials = 0,
#'                                        n_sim_trials = 25000,
#'                                        sim_algorithm = sim_algorithm,
#'                                        sim_ctrl = sim_ctrl,
#'                                        num_cores = NULL,
#'                                        no_bootstrapping = TRUE)
#' @export
#'
#'
run_cp_bootstrap <-   function (sim_data, boot_trials = 100, n_sim_trials = 50,
                                sim_algorithm="simple", sim_ctrl = NULL, num_cores = NULL,
                                no_bootstrapping = FALSE)   {
  simulation_data <- sim_data


  # Add an warning if you specify no boostrapping but set boot_trials >0
  if ( boot_trials > 0 & no_bootstrapping == TRUE)
    stop("If you specify 'no_boostrapping' == TRUE, you have to set 'boot_trial' to 0")

  if (boot_trials>0 & no_bootstrapping == FALSE){
    # set up clusters
    if (is.null(num_cores)) {
      num_cores <- min(boot_trials, parallel::detectCores() - 1)
    } else {
      num_cores <- num_cores
    }

    cluster <- parallel::makeCluster(num_cores)
    
    #Pass things to the clusters that we need
    parallel::clusterCall(cluster, function() library(tidyverse))
    parallel::clusterCall(cluster, function() library(delaySim))
    parallel::clusterExport(cluster, varlist=c("boot_trials","simulation_data",
                                               "n_sim_trials"),
                            envir=environment())

      tmp <- parallel::parLapply(cl = cluster,
                                 1:boot_trials,
                                 function(x){boot_change_point(sim_data = simulation_data,
                                                               n_sim_trials = n_sim_trials)})

    # pull out change points
    change_point <- map2(tmp,1:boot_trials,~.x$change_point %>%
                           mutate(boot_trial=.y)) %>%
      bind_rows() %>%
      select(boot_trial,dplyr::everything())

    # pull out miss visit counts
    miss_counts <- map2(tmp,1:boot_trials,~.$miss_counts %>%
                                mutate(boot_trial=.y)) %>%
      bind_rows() %>%
      select(boot_trial,dplyr::everything())

    # pull out predicted visits
    preds <- map2(tmp,1:boot_trials,~.$pred %>%
                    mutate(boot_trial=.y)) %>%
      bind_rows() %>%
      select(boot_trial,dplyr::everything())

    # pull out simulation results
    sim_visit_results <- map2(tmp,1:boot_trials,~.$sim_visit_results %>%
                                mutate(boot_trial=.y))  %>%
      bind_rows() %>%
      select(boot_trial,dplyr::everything())

    parallel::stopCluster(cluster)
    gc()

  } else {
    
      tmp <- boot_change_point(sim_data = simulation_data,
                               n_sim_trials = n_sim_trials,
                               no_bootstrapping = TRUE,
                               num_cores = num_cores)
    

    # pull out change points
    change_point <- tmp$change_point %>%
                           mutate(boot_trial=0) %>%
      bind_rows() %>%
      select(boot_trial,dplyr::everything())

    # pull out miss visit counts
    miss_counts <- tmp$miss_counts %>%
                          mutate(boot_trial=0) %>%
      bind_rows() %>%
      select(boot_trial,dplyr::everything())

    # pull out predicted visits
    preds <- tmp$pred %>%
                    mutate(boot_trial=0) %>%
      bind_rows() %>%
      select(boot_trial,dplyr::everything())

    # pull out simulation results
    sim_visit_results <- tmp$sim_visit_results %>%
                                mutate(boot_trial=0)  %>%
      bind_rows() %>%
      select(boot_trial,dplyr::everything())

  }
  return(list(change_point = change_point,
              miss_counts = miss_counts,
              preds = preds,
              sim_visit_results = sim_visit_results))
}

