
#### This file contains the helper functions and algorithms for running the simulations ###

#' Draw Missed Visits Using the simple algorithm
#'
#'
#' @export
run_sim_draw_simple <- function(sim_data){
  out <- sim_data$time_map %>%
    dplyr::filter(miss_ind==1) %>%  # filter to miss (i.e., SSD) visits
    dplyr::inner_join(dplyr::select(sim_data$miss_bins_visits,period, num_miss), by = "period") %>% # join in number missed by bin
    dplyr::mutate(rand = runif(n = dplyr::n())) %>% # draw random number
    dplyr::arrange(period, rand) %>% # arrange by random draw
    dplyr::group_by(period) %>%
    dplyr::filter(dplyr::row_number() <= num_miss) %>% # filter to the number in each bin corresponding to number missed
    dplyr::ungroup() %>%
    dplyr::select(-num_miss,-rand)

  return(out)
}

#' Draw Missed Visits Using the simple correlated draw algorithm
#'
#' @param sim_ctrl the control function specifying the simulation parameters. This simulation uses
#'                 the parameter "alpha"
#'
#' @export
run_sim_draw_simple_cor <- function(sim_data,sim_ctrl){

  draw_set <- c()
  new_set <- list()
  prior_set <- list()


  for (i in sim_data$change_point:0){

    # set of visits at given period from patients previously drawn
    prior_draw <- sim_data$time_map %>%
      dplyr::filter(period==i, patient_id %in% draw_set) %>%
      dplyr::filter(miss_ind==1)   # filter to miss (i.e., SSD) visits

    # set of visits at given period from patients not previously drawn
    new_draw <- sim_data$time_map %>%
      dplyr::filter(period==i, !(patient_id %in% draw_set)) %>%
      dplyr::filter(miss_ind==1)  # filter to miss (i.e., SSD) visits

    # number of misses to draw
    num_miss_visits <- sim_data$miss_bins_visits[sim_data$miss_bins_visits$period==i,"num_miss"]

    # number of misses to draw from each group
    miss_draw_prior <- num_miss_visits*(1-sim_ctrl$alpha)
    miss_draw_new <- num_miss_visits*sim_ctrl$alpha

    # correct the draw counts if there are not enough cases to draw from
    if (miss_draw_prior>nrow(prior_draw)){
      # if not enough prior visits apply correction to miss_draw_new
      correction <- miss_draw_prior-nrow(prior_draw)
      miss_draw_new <- miss_draw_new + correction

    } else if (miss_draw_new>nrow(new_draw)){
      # if not enough new visits apply correction to miss_draw_prior
      correction <- miss_draw_new-nrow(new_draw)
      miss_draw_prior <- miss_draw_prior + correction

    }

    # cases drawn from previously drawn patients
    out1 <- prior_draw %>%
      dplyr::mutate(rand = runif(n = dplyr::n())) %>% # draw random number
      dplyr::arrange(rand) %>% # arrange by random draw
      dplyr::filter(dplyr::row_number() <= miss_draw_prior)  # filter to the number in each bin corresponding to number missed

    # cases drawn from new patients
    out2 <- new_draw %>%
      dplyr::mutate(rand = runif(n = dplyr::n())) %>% # draw random number
      dplyr::arrange(rand) %>% # arrange by random draw
      dplyr::filter(dplyr::row_number() <= miss_draw_new)  # filter to the number in each bin corresponding to number missed

    # add in patients previously drawn
    draw_set <- c(draw_set,unique(out1$patient_id),unique(out2$patient_id)) %>% unique()

    # update sets with patients drawn
    new_set[[i+1]] <- out2
    prior_set[[i+1]] <- out1

  }

  out <- dplyr::bind_rows(new_set,prior_set) %>%
    dplyr::select(-rand)

  return(out)
}
