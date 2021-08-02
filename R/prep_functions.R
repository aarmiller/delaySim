
#' Count prior visits assocatied with a given condition for preparation in finding a changepoint
#'
#' @importFrom rlang .data
#'
#' @param time_map_data A timemap with indicators for specific events (i.e.
#' diagnoses during a given visit)
#' @param event_name The variable name for the even indicator
#' @param start_day When to start counting prior to the index
#' @param by_days The number of days in each bin to count by
#'
#'
#' @export
prep_cp_data <- function(time_map_data,event_name,start_day=1L,by_days=1L){
  
  tmp <- time_map_data %>%
    dplyr::rename(miss_ind=!!event_name) %>%
    dplyr::filter(.data$days_since_dx<=-start_day) %>%
    dplyr::arrange(-.data$days_since_dx) %>%
    dplyr::mutate(period=by_days*((-.data$days_since_dx-start_day) %/% by_days)) %>%
    dplyr::group_by(.data$period)
  
  tmp1 <- tmp %>%
    dplyr::summarise(n_visits=dplyr::n(),
                     n_patients=dplyr::n_distinct(.data$enrolid))
  
  tmp2 <- tmp %>%
    dplyr::filter(miss_ind==1) %>%
    dplyr::summarise(n_miss_visits=dplyr::n(),
                     n_miss_patients=dplyr::n_distinct(.data$enrolid))
  
  dplyr::inner_join(tmp1,tmp2,by="period")
}




#'  Prepare data for simulation
#'
#' @param time_map_data a time_map of visits
#' @param by_days the number of days to aggregate by in counting periods
#' @param start_day When to start counting prior to the index
#' @param event_name The variable name for the event indicator
#' @param cp_method The change-point method to fit visit counts (i.e. "lm","lm_quad","lm_cube", "quad", "cube", "exp", "spline", "cusum",
#' "pettitt", or "set_cp"). "set_cp" is not a change point detection method, rather it is a method to specify the change point and method used
#' to model the data prior to the specified change point
#' @param eval_criteria The evaluation criteria used to find change points, if using a
#' linear regression method
#' @param specify_cp Set a specific change point you want to use instead of searching for optimal change point. Enter a postive integer value
#' repersenting the days before the index on which you you want to specify the change point. (e.g. 100 would be 100 days before the index)
#' This is a required argument if cp_method = "set_cp"
#' @param week_period Logical to incorporate a "day of the week" effect into the linear model, if
#' method is "pettitt" of "cusum". Note this is only sensible for one-day period aggregation
#' @param auto_reg Logical that determines whether expected counts use a time-series framework that incorporates autoregression.
#' Will automatically fit periodicity, automatically setting week_period to TRUE. Only relevant for cusum and pettitt methods
#' @param set_cp_method The method used to model the data prior to a specified change point for the "set_cp" cp_method
#' (i.e. "linear", "cubic", "exponential", etc.)
#' @param prediction_bound_for_sim Logical to specify whether or not to use the estimated predicted value or the upper bound 90%
#' prediction value in the simulations. The defualt is FALSE which uses the estimated predicited value
#' 
#' @export
#' 
#' @examples
#'
#' ## Example to detect change point using cusum ##
#' out <- prep_sim_data_delay(time_map_data, by_days=1, start_day=1, event_name = "any_ssd", cp_method = "cumsum")
#'
#' ## Examples showing different ways to specify a change point ##
#'
#' # Using a change point detection method (e.g. cumsum or lm_cube), but enforce a specific change point
#' cp <- 21
#' out <- prep_sim_data_delaydx(time_map_data, by_days=1, start_day=1, event_name = "any_ssd", cp_method = "cumsum", specify_cp = cp)
#'
#' # Using the "set_cp" in the cp_method argument to apply a cubic model to the data prior to the change point
#' cp <- 21
#' out <- prep_sim_data_delaydx(time_map_data, by_days=1, start_day=1, event_name = "any_ssd", cp_method = "set_cp", specify_cp = cp,
#' set_cp_method = "cubic")
#'
#'
#'
prep_sim_data <- function(time_map_data, by_days=1, start_day=1, time_before=-365, event_name = "any_ssd", cp_method = "lm_quad", specify_cp = NULL,
                                  set_cp_method = NULL, eval_criteria="AIC", week_period=FALSE, prediction_bound_for_sim = FALSE,
                                  auto_reg=FALSE){

  if (cp_method == "set_cp" & (is.null(specify_cp) | is.null(set_cp_method))){
    stop("If using the 'set_cp' method for cp_method, specify_cp and set_cp_method cannot be NULL")
  }

  #Filter time map to be in format that we desire
  sim_time_map <- time_map_data %>% dplyr::mutate(period = ((-days_since_dx - days_before)%/%by_days)) %>%
    dplyr::mutate(enrolid_new = enrolid) %>%
    dplyr::select(enrolid, enrolid_new, period, days_since_dx,
                  miss_ind = (!!sym(event_name)),  inpatient,
                  ed, outpatient) %>% filter(period >= 0)

  #Get miss bins and changepoint

  final_prior_visit_counts <- prep_cp_data(time_map_data,
                                                        event_name = event_name,
                                                        start_day = start_day,
                                                        by_days = by_days) %>%
                                       filter(period<=time_before) %>%
                                       mutate(days_since_dx = -period - by_days)

  if (cp_method != "set_cp"){
    results <- final_prior_visit_counts %>%
      find_change_point(var_name = "n_miss_visits",
                        method = cp_method,
                        return_miss_only = FALSE,
                        eval_criteria = eval_criteria,
                        week_period = week_period,
                        auto_reg = auto_reg,
                        specify_cp = specify_cp)


  } else{
    results <- final_prior_visit_counts %>%
      set_change_point(var_name = "n_miss_visits",
                       method = set_cp_method,
                       return_miss_only = FALSE,
                       compare_all_methods = FALSE,
                       week_period = week_period,
                       specify_cp = specify_cp)

  }

  #Extract what we want from each of these
  change_point <- results$change_point %>% .$period
  miss_bins_visits <- results$pred %>%  filter(period <= change_point) %>% mutate(num_miss = Y - pred1) %>%
    mutate(num_miss = ifelse(num_miss < 0, 0, num_miss)) %>% select(period, Y, pred, pred1, num_miss) %>%
    mutate(num_miss = round(num_miss, 0))



  #Return a list with all data we need
  #Also return parameters with which this was called for future reference
  return(list(time_map = sim_time_map,
              miss_bins_visits = miss_bins_visits,
              total_patients  = time_map_data %>% dplyr::distinct(enrolid) %>% nrow(),
              change_point = change_point,

              cp_method=cp_method,
              event_name=event_name,
              start_day=start_day,
              eval_criteria=eval_criteria,
              week_period=week_period,
              specify_cp = specify_cp,
              auto_reg = auto_reg,
              set_cp_method = set_cp_method,
              prediction_bound_for_sim = prediction_bound_for_sim)
  )


}

