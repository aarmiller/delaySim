require(tidyverse)
require(lubridate)
require(stringr)
require(readr)
#require(codeBuildr)

#' Parse ASCII Files
#'
#' @param path (Character) Path to ASCII file
#' @param file_type (Character) Which file to parse (e.g., core, ahal, severity)
#' @param data_source (Character) The data set source, either SEDD or SID
#' @param state (Character) Two-letter abbervation for the desired state
#' @param year (Character) Typically, a four digit expression of desired year
#' but may be 2015q1q3 or 2015q4 to accomdate the switch from ICD-9 to ICD-10
#' on 2015/10/01
#' @param n_max (Numeric) Number of rows to read in. Defaults to infinity.
#' Useful for testing.
#'
#' @return A tibble of the parsed data from the data located at \code{path}
#' using the parameters detailed by \code{file_type}, \code{data_source},
#' \code{state}, \code{year} and the version guessed by
#' \code{find_data_version()}.
#'
#' @export
#'
parse_data <- function(path, file_type, data_source, state, year, n_max = Inf) {
  f <- tolower(file_type)
  badyear <- FALSE
  if (state == "NC" & file_type == "AHAL") {
    if (year == 2006) {
      year <- 2005
    } else if (year == 2007) {
      year <- 2008
    } else {
      year <- year
    }
    badyear <- TRUE
  }
  if (state == "FL" & year == 2007 & file_type == "AHAL") {
    year <- 2006
    badyear <- TRUE
  }
  
  data_version <- find_data_version(path, file_type, data_source, state, year)
  locations <- find_variable_locations(data_source, state, year)
  
  if (purrr::is_empty(data_version)) {
    data_version <- 0
  }
  
  locations <- locations %>%
    dplyr::filter(version == data_version,
                  file_type == f)
  
  if (badyear) {
    if (state == "NC" & year == 2005) {
      year <- 2005
    } else if (state == "NC" & year == 2008) {
      year <- 2007
    } else {
      year <- 2007
    }
    year <- ifelse(state == "FL", "2007", "2006")
  }
  
  # check if the first two lines are HCUP notices or not
  notice <- readr::read_lines(path, n_max = 1) %>%
    stringr::str_detect("NOTICE")
  
  data <- readr::read_fwf(path,
                          readr::fwf_positions(start = locations$start,
                                               end = locations$end,
                                               col_names = locations$variable),
                          col_types = stringr::str_flatten(ifelse(
                            locations$type == "Char" | tolower(locations$variable) == "key",
                            "c", "n")),
                          skip = ifelse(notice, 2, 0),
                          n_max = n_max)
  
  return(data)
}

#' Determine Data Version Number
#'
#' @param path (Character) Path to CORE ASCII file
#' @param file_type (Character) Which file to parse (e.g., core, ahal, severity)
#' @param data_source (Character) The data set source, either SEDD or SID
#' @param state (Character) Two-letter abbervation for the desired state
#' @param year (Character) Typically, a four digit expression of desired year
#' but may be 2015q1q3 or 2015q4 to accomdate the switch from ICD-9 to ICD-10
#' on 2015/10/01
#'
#' @return a vector of length 1 with the version number of the data located at
#' \code{path} based on the observed line-length and the expected line length
#' for the different versions of data provided for \code{data_source},
#' \code{state} and \code{year}. A version of \code{0} denotes the current
#' version while values of \code{-1}, \code{-2} and so on denote versions
#' 1, 2 and so on.
#'
#' @export
#'
find_data_version <- function(path, file_type, data_source, state, year) {
  s <- state
  y <- year
  f <- tolower(file_type)
  if (data_source == "SID") {
    expected_characters <- sid_meta_data
  } else {
    return(0) # we only have the most current version of the SEDD data sets
  }
  expected_characters <- expected_characters %>%
    dplyr::filter(name == "line_lengths") %>%
    tidyr::unnest() %>%
    dplyr::filter(tolower(file_type) == f,
                  state == s,
                  year == y)
  
  # you need to skip the first few lines because, from time-to-time and likely
  # with some glee, HCUP appends a NOTICE to the start of the data. If we skip
  # the first few lines, then we are certainly in the data in all years.
  observed_characters <- stringr::str_length(readr::read_lines(path,
                                                               skip = 5,
                                                               n_max = 1))
  version <- expected_characters %>%
    dplyr::filter(length == observed_characters)
  version <- version$version
  return(version)
}

#' Find Variables and Locations by Data Souce, State and Year
#'
#' @param data_source character The data set source, either SEDD or SID
#' @param state character Two-letter abbervation for the desired state
#' @param year character Typically, a four digit expression of desired year
#' but may be 2015q1q3 or 2015q4 to accomdate the switch from ICD-9 to ICD-10
#' on 2015/10/01
#'
#' @return a tibble with columns for variable name, encoding type, start and
#' end columns referring to the ASCII file for data source, state and year.
#' @export

find_variable_locations <- function(data_source, state, year) {
  s <- state
  y <- year
  if (data_source == "SID") {
    variables_and_locations <- sid_meta_data
  } else {
    variables_and_locations <- sedd_meta_data
  }
  variables_and_locations <- variables_and_locations %>%
    dplyr::filter(name == "locations") %>%
    tidyr::unnest() %>%
    dplyr::filter(state == s,
                  year == y) %>%
    dplyr::select(variable, file_type, version, start, end, type) %>%
    dplyr::mutate(file_type = tolower(file_type))
  return(variables_and_locations)
}


#' Get SEDD helper functions used in build_HCUP_time_map
#'
#' @param data (tibble) Typically, a tibble containing the pasrsed SID data from the ASCII file
#' @param index_codes (list) Typically, a list length 2 for ICD-9-CM and ICD-10-CM codes for the condition of interest
#' @param ssd_list (tibble) Typically, a tibble containing the SSD names and their respective ICD-9-CM and ICD-10-CM codes
#' @return a tibble for the time map of interest
#'
#' @export
 

get_SEDD <- function(data, index_codes, ssd_list){
  
  temp1 <- data %>% 
    select(KEY, VisitLink, DaysToEvent,
           ATYPE, AWEEKEND,
           AYEAR, AMONTH, AWEEKEND,
           DYEAR = YEAR, DMONTH,
           LOS_X, DIED, 
           AGE, FEMALE,  RACE, 
           matches("^DX[[:digit:]]|^I10_DX[[:digit:]]|^DX_A"),
           matches("^CPT[[:digit:]]")) %>% 
    mutate(daystoend = DaysToEvent + LOS_X,
           SID = 0L)
  
  #find relevant visit keys
  visit_keys_and_KEYs <- temp1 %>% 
    select(KEY, VisitLink, DaysToEvent, AMONTH, contains("DX")) %>% 
    pivot_longer(contains("DX"), values_to = "DX") %>% 
    filter(DX %in% c(index_codes$icd9_codes, index_codes$icd10_codes)) %>% 
    distinct(VisitLink, KEY) %>% 
    mutate(SID = 0L) %>% 
    distinct(VisitLink, KEY, SID) %>% 
    mutate(dx_ind = 1L)
  
  # select visits for specific keys 
  temp2 <-  temp1 %>% inner_join(visit_keys_and_KEYs %>%
                                   filter(SID == 0L) %>% 
                                   distinct(VisitLink), by = "VisitLink") %>% 
    left_join(visit_keys_and_KEYs %>%
                filter(SID == 0L) %>% 
                distinct(VisitLink, KEY, dx_ind))
  
  temp_ssd <-  temp2 %>% 
    select(KEY, VisitLink, contains("DX"), -dx_ind) %>% 
    pivot_longer(contains("DX"), values_to = "DX") 
  
  for(i in unique(ssd_list$ssd_name)){
    temp_codes <- ssd_list %>% filter(ssd_name == i) %>% .$icd_codes
    temp_ssd[[i]] <- ifelse(temp_ssd$DX %in% temp_codes, 1L, 0L)
  }
  
  temp_ssd <- temp_ssd %>% select(-name, -DX) %>% 
    group_by(KEY, VisitLink) %>% 
    summarise_at(vars(unique(ssd_list$ssd_name)), .funs = ~max(.)) %>% 
    ungroup()
  
  out <- temp2 %>% left_join(temp_ssd) %>% 
    mutate_at(vars(c(unique(ssd_list$ssd_name),dx_ind)), .funs = ~ifelse(is.na(.), 0L, .)) %>% 
    select(-contains("DX"), -contains("CPT"), dx_ind)
  
  return(out)
}

#' Get SID helper functions used in build_HCUP_time_map
#'
#' @param data (tibble) Typically, a tibble containing the pasrsed SID data from the ASCII file
#' @param index_codes (list) Typically, a list length 2 for ICD-9-CM and ICD-10-CM codes for the condition of interest
#' @param ssd_list (tibble) Typically, a tibble containing the SSD names and their respective ICD-9-CM and ICD-10-CM codes
#' @return a tibble for the time map of interest
#'
#' @export

get_SID <- function(data, index_codes, ssd_list){
  
  temp1 <- data %>% 
    select(KEY, VisitLink, DaysToEvent,
           ATYPE, AWEEKEND,
           AYEAR, AMONTH, AWEEKEND,
           DYEAR = YEAR, DMONTH,
           LOS_X, DIED, 
           AGE, FEMALE,  RACE, 
           matches("^DX[[:digit:]]|^I10_DX[[:digit:]]|^DX_A"),
           matches("^PR[[:digit:]]|^I10_PR[[:digit:]]")) %>% 
    mutate(daystoend = DaysToEvent + LOS_X,
           SID = 1L)
  
  #find relevant visit keys
  visit_keys_and_KEYs <- temp1 %>% 
    select(KEY, VisitLink, DaysToEvent, AMONTH, contains("DX")) %>% 
    pivot_longer(contains("DX"), values_to = "DX") %>% 
    filter(DX %in% c(index_codes$icd9_codes, index_codes$icd10_codes)) %>% 
    distinct(VisitLink, KEY) %>% 
    mutate(SID = 1L) %>% 
    distinct(VisitLink, KEY, SID) %>% 
    mutate(dx_ind = 1L)
  
  # select visits for specific keys 
  temp2 <-  temp1 %>% inner_join(visit_keys_and_KEYs %>%
                                   filter(SID == 1L) %>% 
                                   distinct(VisitLink), by = "VisitLink") %>% 
    left_join(visit_keys_and_KEYs %>%
                filter(SID == 1L) %>% 
                distinct(VisitLink, KEY, dx_ind))
  
  temp_ssd <-  temp2 %>% 
    select(KEY, VisitLink, contains("DX"), -dx_ind) %>% 
    pivot_longer(contains("DX"), values_to = "DX") 
  
  for(i in unique(ssd_list$ssd_name)){
    temp_codes <- ssd_list %>% filter(ssd_name == i) %>% .$icd_codes
    temp_ssd[[i]] <- ifelse(temp_ssd$DX %in% temp_codes, 1L, 0L)
  }
  
  temp_ssd <- temp_ssd %>% select(-name, -DX) %>% 
    group_by(KEY, VisitLink) %>% 
    summarise_at(vars(unique(ssd_list$ssd_name)), .funs = ~max(.)) %>% 
    ungroup()
  
  out <- temp2 %>% left_join(temp_ssd) %>% 
    mutate_at(vars(c(unique(ssd_list$ssd_name),dx_ind)), .funs = ~ifelse(is.na(.), 0L, .)) %>% 
    select(-contains("DX"), -contains("CPT"), dx_ind)
  
  return(out)
}

#' Build Delay Time Map with ASCII HCUP Data
#'
#' @param path (Character) Path to ASCII files
#' @param state (Character) Two-letter abbervation for the desired state
#' @param year (Character) Typically, a four digit expression of desired year
#' @param index_codes (list) Typically, a list length 2 for ICD-9-CM and ICD-10-CM codes for the condition of interest
#' @param ssd_list (tibble) Typically, a tibble containing the SSD names and their respective ICD-9-CM and ICD-10-CM codes
#' @param n_collect (Numeric) Number of rows to read in. Defaults to infinity.
#' @return a tibble for the time map of interest
#' 
#' @examples
#' index_codes <- codeBuildr::load_disease_codes("ami") # get index codes
#' ssd_list <- read_csv("/Shared/Statepi_Diagnosis/atlan/delay_dx/params/ssd_codes/ami/ssd_codes.csv") #load ssd codes
#' ami_time_map <- build_HCUP_time_map(path = "/Shared/Statepi_Marketscan/raw_data/hcup/",
#'                                     state = "IA",
#'                                     year = "2017",
#'                                     index_codes = index_codes,
#'                                     ssd_list = ssd_list,
#'                                     n_collect = 10000)
#'
#'
#' @export


build_HCUP_time_map <- function(path, state, year, index_codes, ssd_list, n_collect){
  
  if (as.numeric(year) != 2015){
    path_SEDD<- paste0(path, state, "/", state, "_SEDD_", year, "_CORE.asc")
    path_SID<- paste0(path, state, "/", state, "_SID_", year, "_CORE.asc")
    SEDD_exist <- 1L
    SID_exist <- 1L
    
    if (file.exists(path_SEDD)){
      SEDD_temp <- parse_data(path = path_SEDD,
                              file_type = "core",
                              data_source = "SEDD",
                              state = state, 
                              year = year,
                              n_max = n_collect)
    } else {
      SEDD_exist <- 0L
    }
    
    if (file.exists(path_SID)){
      SID_temp <- parse_data(path = path_SID,
                             file_type = "core",
                             data_source = "SID",
                             state = state, 
                             year = year,
                             n_max = n_collect)
    } else {
      SID_exist <- 0L
    }
    
    if((SEDD_exist + SID_exist) <1L){
      stop(paste0("The following files do not exist: \n", path_SEDD, "\n", path_SID))
    } else if (SEDD_exist == 1L & SID_exist == 0L) {
      warning(paste0("The following files do not exist: \n", path_SID))
    } else if (SEDD_exist == 0L & SID_exist == 1L) {
      warning(paste0("The following files do not exist: \n", path_SEDD))
    } 
    
  } else {
    path_SEDD_q1q3<- paste0(path, state, "/", state, "_SEDD_", year, "q1q3_CORE.asc")
    path_SEDD_q4<- paste0(path, state, "/", state, "_SEDD_", year, "q4_CORE.asc")
    path_SID_q1q3<- paste0(path, state, "/", state, "_SID_", year, "q1q3_CORE.asc")
    path_SID_q4<- paste0(path, state, "/", state, "_SID_", year, "q4_CORE.asc")
    SEDD_q1q3_exist <- 1L
    SEDD_q4_exist <- 1L
    SID_q1q3_exist <- 1L
    SID_q4_exist <- 1L
    
    if (file.exists(path_SEDD_q1q3)){
      SEDD_temp_q1q3 <- parse_data(path = path_SEDD_q1q3,
                                   file_type = "core",
                                   data_source = "SEDD",
                                   state = state, 
                                   year = paste0(year,"q1q3"),
                                   n_max = n_collect)
    } else {
      SEDD_q1q3_exist <- 0L
    }
    
    if (file.exists(path_SEDD_q4)){
      SEDD_temp_q4 <- parse_data(path = path_SEDD_q4,
                                 file_type = "core",
                                 data_source = "SEDD",
                                 state = state, 
                                 year = paste0(year,"q4"),
                                 n_max = n_collect)
    } else {
      SEDD_q4_exist <- 0L
    }
    
    if (file.exists(path_SID_q1q3)){
      SID_temp_q1q3 <- parse_data(path = path_SID_q1q3,
                                  file_type = "core",
                                  data_source = "SID",
                                  state = state, 
                                  year = paste0(year,"q1q3"),
                                  n_max = n_collect)
    } else {
      SID_q1q3_exist <- 0L
    }
    
    if (file.exists(path_SID_q4)){
      SID_temp_q4 <- parse_data(path = path_SID_q4,
                                file_type = "core",
                                data_source = "SID",
                                state = state, 
                                year = paste0(year,"q4"),
                                n_max = n_collect)
    } else {
      SID_q4_exist <- 0L
    }
    
    
    if((SEDD_q1q3_exist + SEDD_q4_exist + SID_q1q3_exist + SID_q4_exist) <1L){
      stop(paste0("The following files do not exist: \n", path_SEDD_q1q3, "\n", path_SEDD_q4,
                  "\n", path_SID_q1q3, "\n", path_SID_q4))
    } else if ((SEDD_q1q3_exist == 1L & SEDD_q4_exist == 1L) & (SID_q1q3_exist == 0L & SID_q4_exist == 0L)) {
      warning(paste0("The following files do not exist: \n", path_SID_q1q3, "\n", path_SID_q4))
      SEDD_temp <- SEDD_temp_q1q3 %>% bind_rows(SEDD_temp_q4)
    } else if ((SEDD_q1q3_exist == 0L & SEDD_q4_exist == 0L) & (SID_q1q3_exist == 1L & SID_q4_exist == 1L))  {
      warning(paste0("The following files do not exist: \n", path_SEDD_q1q3, "\n", path_SEDD_q4))
      SID_temp <- SID_temp_q1q3 %>% bind_rows(SID_temp_q4)
    } else {
      SEDD_temp <- SEDD_temp_q1q3 %>% bind_rows(SEDD_temp_q4)
      SID_temp <- SID_temp_q1q3 %>% bind_rows(SID_temp_q4)
    } 
    
  }
    
    # refine ssd list
    temp_ssd_list <- ssd_list
    temp_ssd_list <- temp_ssd_list %>% 
      mutate(ssd_name = str_replace_all(temp_ssd_list$label, " ", "_")) %>% 
      mutate(ssd_name = str_replace_all(temp_ssd_list$ssd_name, ";", "_"))
    
    if(exists("SEDD_temp")){
      # get SEDD
      SEDD_out <- get_SEDD(data = SEDD_temp,
                           index_codes = index_codes,
                           ssd_list = temp_ssd_list)
    }
    
    if(exists("SID_temp")){
     # get SID
      SID_out <- get_SID(data = SID_temp,
                           index_codes = index_codes,
                           ssd_list = temp_ssd_list)
    }
    
    if (exists("SEDD_temp") & exists("SID_temp")){
      
      final_out <- SEDD_out %>% 
        bind_rows(SID_out) %>% 
        arrange(VisitLink, DaysToEvent, SID) %>% 
        group_by(VisitLink) %>%
        mutate(days_since_prev = DaysToEvent - lag(daystoend),
               index_dx = ifelse((dx_ind == 1 & cumsum(dx_ind) ==  1), 1L, 0L),
               days_since_index_dx = DaysToEvent - max(index_dx * DaysToEvent)) %>% 
        ungroup() %>% 
        distinct()
      
    } else if (exists("SEDD_temp") & !exists("SID_temp")){
  
      final_out <- SEDD_out %>% 
        arrange(VisitLink, DaysToEvent, SID) %>% 
        group_by(VisitLink) %>%
        mutate(days_since_prev = DaysToEvent - lag(daystoend),
               index_dx = ifelse((dx_ind == 1 & cumsum(dx_ind) ==  1), 1L, 0L),
               days_since_index_dx = DaysToEvent - max(index_dx * DaysToEvent)) %>% 
        ungroup() %>% 
        distinct()
      
    } else if (!exists("SEDD_temp") & exists("SID_temp")){
      
      final_out <- SID_out %>% 
        arrange(VisitLink, DaysToEvent, SID) %>% 
        group_by(VisitLink) %>%
        mutate(days_since_prev = DaysToEvent - lag(daystoend),
               index_dx = ifelse((dx_ind == 1 & cumsum(dx_ind) ==  1), 1L, 0L),
               days_since_index_dx = DaysToEvent - max(index_dx * DaysToEvent)) %>% 
        ungroup() %>% 
        distinct()
      
    }
    
    return(final_out)
}


