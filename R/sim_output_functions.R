#### This file contains functions to create output of the simulation ####

#' Compute simulation trial miss summary
#'
#'
#' @export
compute_miss_summary <- function(sim_miss_num_data){

  out <- sim_miss_num_data %>%
    dplyr::summarise(n_pat = dplyr::n(),                          # number of patients missed
                     mean_n_vis = mean(n_vis),             # mean number of misses per patient
                     median_n_vis =median(n_vis),          # median number of misses per patient
                     mean_n_vis_out = mean(n_vis_out),     # mean number of ed misses per patient
                     median_n_vis_out =median(n_vis_out),  # median number of ed misses per patient
                     mean_n_vis_ed = mean(n_vis_ed),       # mean number of ed misses per patient
                     median_n_vis_ed =median(n_vis_ed),    # median number of ed misses per patient
                     mean_n_vis_inpatient = mean(n_vis_inpatient),     # mean number of ed misses per patient
                     median_n_vis_inpatient = median(n_vis_inpatient),  # median number of ed misses per patient
                     min_dur = min(dur),           # min duration
                     mean_dur = mean(dur),         # mean miss duration
                     median_dur = median(dur),     # median miss duration
                     max_dur = max(dur),
                     n_vis = sum(n_vis),
                     n_vis_out = sum(n_vis_out),
                     n_vis_ed = sum(n_vis_ed),
                     n_vis_inpatient = sum(n_vis_inpatient))

  return(out)
}

#' Compute number of missed visits per patient
#'
#'
#' @export
compute_miss_num <- function(miss_draw_data){
  miss_draw_data %>%
    dplyr::group_by(patient_id) %>%
    dplyr::summarise(n_vis = dplyr::n(),
                     n_vis_out = sum(ed==0 & inpatient==0),
                     n_vis_ed = sum(ed==1),
                     n_vis_inpatient = sum(inpatient==1),
                     dur = max(-(days_since_dx)))
}


#' Compute simulation trial visit table
#'
#'
#' @export
compute_sim_trial_n_visit_table <- function(sim_miss_num_data,total_patients){
  ## Compute bins for number of missed visits ##
  visit <- tibble::tibble(n_vis_groups = 1:5) %>%
    dplyr::mutate(n=purrr::map(n_vis_groups, ~ dplyr::filter(sim_miss_num_data, n_vis >= .) %>% nrow())) %>%
    tidyr::unnest(n) %>%
    dplyr::mutate(percent = n/total_patients*100)

  ## Aggregate number of missed visits table ##
  n_miss <- visit %>%
    rbind(visit %>%
            dplyr::filter(n_vis_groups==1) %>%
            dplyr::mutate(n_zero=total_patients-n, percent_zero = (n_zero/total_patients)*100) %>%
            dplyr::mutate(n_vis_groups = 0)  %>%
            dplyr::select(n_vis_groups, n = n_zero, percent = percent_zero)) %>%
    dplyr::arrange(n_vis_groups) %>%
    dplyr::mutate(n_vis_groups=ifelse(n_vis_groups !=0, paste0(">= ", n_vis_groups), as.character(n_vis_groups))) %>%
    rbind(sim_miss_num_data %>%
            dplyr::mutate(n_vis_groups=ifelse(n_vis>=mean(sim_miss_num_data$n_vis),
                                              ">= Mean", "Greater that or equal to mean")) %>%
            dplyr::group_by(n_vis_groups) %>%
            dplyr::count() %>%
            dplyr::ungroup() %>%
            dplyr::mutate(percent=n/total_patients*100) %>%
            dplyr::filter(n_vis_groups==">= Mean")) %>%
    rbind(sim_miss_num_data %>%
            dplyr::mutate(n_vis_groups=ifelse(n_vis>=median(sim_miss_num_data$n_vis),
                                              ">= Median", "Greater that or equal to median")) %>%
            dplyr::group_by(n_vis_groups) %>%
            dplyr::count() %>%
            dplyr::ungroup() %>%
            dplyr::mutate(percent=n/total_patients*100) %>%
            dplyr::filter(n_vis_groups==">= Median"))

  return(n_miss)
}

#' Compute simulation trial duration table
#'
#'
#' @export
compute_sim_trial_duration_table <- function(sim_miss_num_data, upper_bound){
  ## compute range of miss durations ##
  # use change point to define upper bound of range
  x <-  cut(0:upper_bound, 10)

  range <- c(0)
  for (i in levels(x)){
    y <- stringr::str_split(i, ",")
    y <- y[[1]][2]
    z <- floor(as.numeric(stringr::str_remove(y, "]")))
    range <- c(range, z)
  }

  ## compute table of miss durations ##
  duration_miss <- tibble::tibble(bins= range) %>%
    dplyr::mutate(n=purrr::map(bins, ~dplyr::filter(sim_miss_num_data, dur>= .) %>%
                                 nrow())) %>%
    tidyr::unnest(n) %>%
    dplyr::mutate(percent=n/nrow(sim_miss_num_data)*100) %>%
    dplyr::mutate(bins = paste0(">= ", bins)) %>%
    rbind(tibble::tibble(bins = ">= Mean",
                         n =dplyr::filter(sim_miss_num_data,dur >= mean(sim_miss_num_data$dur)) %>%
                           nrow(),
                         percent = n/nrow(sim_miss_num_data)*100)) %>%
    rbind(tibble::tibble(bins = ">= Median",
                         n =dplyr::filter(sim_miss_num_data, dur >= median(sim_miss_num_data$dur)) %>%
                           nrow(),
                         percent = n/nrow(sim_miss_num_data)*100))

  return(duration_miss)
}


