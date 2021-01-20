#' Simulated miss visits from an estimated set of miss bins
#'
#' @param sim_data a dataset containing the time_map, miss bins and other parameters used for the simulation.
#'                 This dataset should be created using the `prep_sim_data()` function
#' @export
#'
sim_miss_visits <- function (sim_data, sim_algorithm="simple", sim_ctrl = NULL) {

  if (is.null(sim_ctrl)){
    sim_ctrl <- sim_ctrl()
  }

  #### Draw missed visits based on selected algorithm ####
  if (sim_algorithm=="simple"){

    # simple algorithm with uncorrelated draw
    miss_draw <- run_sim_draw_simple(sim_data=sim_data)

  } else if (sim_algorithm=="simple_correlated") {

    # simple algorithm with correlated draw
    miss_draw <- run_sim_draw_simple_cor(sim_data=sim_data,sim_ctrl = sim_ctrl)

  } else if (sim_algorithm=="general") {
    
    miss_draw <- run_sim_draw_general(sim_data=sim_data,sim_ctrl = sim_ctrl)
    
  }

  ## Compute number missed and miss duration by patient ##
  sim_miss_num <- compute_miss_num(miss_draw)

  # for computing mean and median duration with the 0 visits
  tmp_num_not_drawn <- sim_data$total_patients-nrow(sim_miss_num)

  # stats when 0 miss patients are included
  w0_stats <- sim_miss_num %>%
    dplyr::full_join(tibble::tibble(patient_id=rep(-99,tmp_num_not_drawn)),by = "patient_id") %>%
    dplyr::mutate(dplyr::across(.cols = dplyr::everything(),~tidyr::replace_na(.,0))) %>%
    dplyr::summarise_at(dplyr::vars(-patient_id),list(mean_w0=mean,median_w0=median))

  ## build table of number of missed visits
  sim_trial_n_visit_table <- compute_sim_trial_n_visit_table(sim_miss_num_data = sim_miss_num,
                                                             total_patients = sim_data$total_patients)


  ## build table of duration of missed visits
  sim_trial_duration_table <- compute_sim_trial_duration_table(sim_miss_num_data = sim_miss_num,
                                                                upper_bound = sim_data$change_point)
  #sim_trial_duration_table <- NULL

  ## Compute summary statistics across all patients ###
  miss_summary <- compute_miss_summary(sim_miss_num_data=sim_miss_num) %>%
    dplyr::bind_cols(w0_stats)

  return(list(sum_n_miss=sim_trial_n_visit_table,
              sum_duration=sim_trial_duration_table,
              miss_summary=miss_summary))
}



#' Control function for specifying simulation parameters
#'
#' @param alpha draw parameter for simple_correlated simulation
#' @param weight_function a weighting function to use in the third simulation algorithm
#'
#' @export
#'
sim_ctrl <- function(alpha=0.5,weight_function=simple_weight_sum){
  list(alpha = alpha,
       weight_function = weight_function)
}

