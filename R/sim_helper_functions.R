
#### This file contains the helper functions and algorithms for running the simulations ###

#' Draw Missed Visits Using the simple algorithm
#'
#' @param sim_data a simulation dataset containing time_map, miss_bins, change_point and number of patients
#'
#' @export
#'
#' @examples
#' # Basic usage of run_sim_draw_simple
#' run_sim_draw_simple(sim_data = example_sim_data)
#'
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
#' @param sim_data a simulation dataset containing time_map, miss_bins, change_point and number of patients
#' @param sim_ctrl the control function specifying the simulation parameters. This simulation uses
#'                 the parameter "alpha." The default value is set to
#'
#' @export
#'
#' @examples
#' # Usage of run_sim_draw_simple_cor() while specifying alpha=0
#' run_sim_draw_simple_cor(sim_data = example_sim_data,
#'                         sim_ctrl = sim_ctrl(alpha = 0))
#'
run_sim_draw_simple_cor <- function(sim_data,sim_ctrl){

  draw_set <- c()
  new_set <- list()
  prior_set <- list()


  for (i in sim_data$change_point:1){

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
      miss_draw_prior <- nrow(prior_draw)  # update the over draw

    } else if (miss_draw_new>nrow(new_draw)){
      # if not enough new visits apply correction to miss_draw_prior
      correction <- miss_draw_new-nrow(new_draw)
      miss_draw_prior <- miss_draw_prior + correction
      miss_draw_new <- nrow(new_draw) # update the over draw

    }

    # cases drawn from previously drawn patients
    out1 <- prior_draw %>%
      dplyr::mutate(rand = runif(n = dplyr::n())) %>% # draw random number
      dplyr::arrange(rand) %>% # arrange by random draw
      dplyr::filter(dplyr::row_number() <= miss_draw_prior[[1]]) # filter to the number in each bin corresponding to number missed

    # cases drawn from new patients
    out2 <- new_draw %>%
      dplyr::mutate(rand = runif(n = dplyr::n())) %>% # draw random number
      dplyr::arrange(rand) %>% # arrange by random draw
      dplyr::filter(dplyr::row_number() <= miss_draw_new[[1]])  # filter to the number in each bin corresponding to number missed

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

#' Draw Missed Visits Using the generalized algorithm
#'
#' @param sim_data a simulation dataset containing time_map, miss_bins, change_point and number of patients
#' @param sim_ctrl the control function specifying the simulation parameters. This simulation uses
#'                 as a parameter a weighting function, weight_function. The default weight_function is
#'                 the additive function simple_weight_sum() which simply adds together the number of distinct
#'                 visits, number of distinct SSDs and the number of times previously selected.
#'
#' @export
#'
#' @examples
#' # Usage of run_sim_draw_simple_cor() while specifying weight_function=simple_weight_sum
#' run_sim_draw_general(sim_data = example_sim_data,
#'                      sim_ctrl = sim_ctrl(weight_function=simple_weight_sum))
#'
run_sim_draw_general <- function(sim_data,sim_ctrl){

  # create data frame of variables to compute weights
  # rows for 1) SSD Visits, 2) Distinct SSDs, 3) Number of times drawn
  # 1 and 2 are created outside the simulation. 3 is done at each step
  weight_data <- sim_data$time_map %>%
    dplyr::filter(miss_ind==1) %>%
    dplyr::arrange(patient_id,days_since_dx) %>%
    dplyr::group_by(patient_id) %>%
    dplyr::mutate_at(dplyr::vars(contains("ssd")),~cumsum(.)>0) %>%
    dplyr::mutate(distinct_visits=dplyr::row_number()) %>% # note the time-map needs to be for distinct visits
    dplyr::ungroup() %>%
    dplyr::mutate(distinct_ssds=rowSums(.[7:10])) %>%
    dplyr::group_by(period) %>%
    tidyr::nest() %>%
    dplyr::arrange(dplyr::desc(period))

  # update weight data at first step compute weights
  weight_data$data[[1]] <- weight_data$data[[1]] %>%
    dplyr::mutate(prior_draws=0) %>%
    dplyr::mutate(weight=sim_ctrl$weight_function(distinct_visits,distinct_ssds,prior_draws)) %>%
    dplyr::mutate(weight2=weight/sum(weight))

  # draw sample - replace 4 with the number needing to be drawn
  tmp_draw <- sample(1:nrow(weight_data$data[[1]]),
                     size = sim_data$miss_bins_visits$num_miss[1],
                     prob = weight_data$data[[1]]$weight2)

  # the final drawn set
  draw_set <- weight_data$data[[1]][tmp_draw,]

  # loop over remaining periods
  for (i in 2:sim_data$change_point){

    # count the number of times a patient has been previously drawn
    tmp_draw_counts <- draw_set %>%
      dplyr::count(patient_id) %>%
      dplyr::rename(prior_draws=n)

    # update time_map for given period with weights
    weight_data$data[[i]] <- weight_data$data[[i]] %>%
      dplyr::left_join(tmp_draw_counts,by = "patient_id") %>%
      dplyr::mutate(prior_draws=ifelse(is.na(prior_draws),0,prior_draws)) %>%
      dplyr::mutate(weight=sim_ctrl$weight_function(distinct_visits,distinct_ssds,prior_draws)) %>%
      dplyr::mutate(weight2=weight/sum(weight))

    # draw sample - replace 4 with the number needing to be drawn
    tmp_draw <- sample(1:nrow(weight_data$data[[i]]),
                       size = sim_data$miss_bins_visits$num_miss[i],
                       prob = weight_data$data[[i]]$weight2)

    # update the final drawn set
    draw_set <- rbind(draw_set,weight_data$data[[i]][tmp_draw,])
  }

  out <- draw_set
  # out <- select(draw_set,
  #               patient_id,days_since_dx)
  return(out)
}

#' simple weight function function for algorithm 3
#'
#' @param num_distinct_visits number of distinct visits up to the time-point
#' @param num_distinct_ssds number of distinct ssds up to the time-point
#' @param num_prior_draws number of times the individual has been previously drawn
#'
#' @export
#'
simple_weight_sum <- function(num_distinct_visits,num_distinct_ssds,num_prior_draws){
  num_distinct_visits+num_distinct_ssds+num_prior_draws
}
