
<!-- README.md is generated from README.Rmd. Please edit that file -->

# delaySim

<!-- badges: start -->
<!-- badges: end -->

The delaySim package performs simulations to estimate the frequency and
duration of diagnostic delays. The package contains three different
algorithms for simulating missed opportunities.

## Installation

You can install the delaySim package using the following:

``` r
# install.packages("devtools")
devtools::install_github("aarmiller/delaySim")
```

Before using the package it is helpful to also load the `tidyverse`
package:

``` r
library(tidyverse)
library(delaySim)
```

## Example

The package comes pre-loaded with an example dataset for performing
example simulations. (To perform simulation using the package, your
dataset should be formatted in a similar manner):

``` r
## View the example dataset
data("example_sim_data")
example_sim_data
#> $time_map
#> # A tibble: 10,875 x 10
#>    patient_id period days_since_dx miss_ind inpatient    ed ssd1  ssd2  ssd3 
#>         <int>  <dbl>         <dbl>    <int>     <int> <int> <lgl> <lgl> <lgl>
#>  1         21     21           -21        1         0     0 FALSE FALSE TRUE 
#>  2         21     20           -20        1         0     0 TRUE  FALSE FALSE
#>  3         21     17           -17        1         0     0 FALSE TRUE  FALSE
#>  4         21     16           -16        1         0     0 FALSE FALSE FALSE
#>  5         21     14           -14        1         0     0 FALSE TRUE  FALSE
#>  6         21     10           -10        1         0     0 TRUE  FALSE FALSE
#>  7         21      9            -9        1         0     0 FALSE FALSE FALSE
#>  8         21      3            -3        1         0     0 FALSE FALSE FALSE
#>  9         21      2            -2        1         0     0 FALSE FALSE FALSE
#> 10         21      1            -1        1         0     0 FALSE TRUE  TRUE 
#> # … with 10,865 more rows, and 1 more variable: ssd4 <lgl>
#> 
#> $miss_bins_visits
#>    period  num_miss
#> 1      21  50.62119
#> 2      20 121.36337
#> 3      19  44.08457
#> 4      18   0.00000
#> 5      17  22.46369
#> 6      16  32.12141
#> 7      15  22.75779
#> 8      14  89.37272
#> 9      13 145.96612
#> 10     12  87.53789
#> 11     11  16.08793
#> 12     10  35.61615
#> 13      9  56.12247
#> 14      8  61.60677
#> 15      7 144.06898
#> 16      6 189.50900
#> 17      5 135.92672
#> 18      4  76.32206
#> 19      3  95.69493
#> 20      2 102.04523
#> 21      1 131.37287
#> 
#> $change_point
#> [1] 21
#> 
#> $total_patients
#> [1] 3371
```

This example dataset is a list with 4 different pieces of information.

1.  `example_sim_data$time_map` - is a time\_map of visits during the
    *diagnostic opportunity window.* This data.frame should contain
    columns with the following names and corresponding information

    -   **patient\_id** - an integer indicating the unique patient id

    -   **period** - a period marker defining each time point in the
        diagnostic opportunity window, with 1 as the time period nearest
        the index diagnosis date and the highest value occurring at the
        change-point period

    -   **days\_since\_dx** - a time variable indicating the number of
        days since diagnosis (these values should be ≤-1

    -   **miss\_ind** - an indicator if the visit contains any SSD
        diagnosis

    -   **inpatient** - an indicator if the visit corresponded to an
        inpatient visits

    -   **ed** - an indicator if the visit corresponded to an inpatient
        visits

    -   **ssd1-ssdx** - indicator variables for the

2.  `example_sim_data$miss_bins_visits` - is a dataset containing the
    number of missed visits that should be drawn at each period during
    the diagnostic opportunity window. This data.frame should contain
    columns with the following names and corresponding information

    -   **period** - a period marker defining each time point in the
        diagnostic opportunity window, with 1 as the time period nearest
        the index diagnosis date and the highest value occurring at the
        change-point period

3.  `example_sim_data$change_point` - is an integer giving the period
    representing the change-point, or furthest period in the diagnostic
    opportunity window

4.  `example_sim_data$total_patients` - is an integer giving the total
    number of patients in the simulation (note: this value must be
    provided since some individuals may not have visits occurring in the
    time\_map provided)

### Simple Uncorrelated Algorithm

The first algorithm draws patients at each time period

``` r
sim_miss_visits(example_sim_data)
#> $sum_n_miss
#> # A tibble: 8 x 3
#>   n_vis_groups     n percent
#>   <chr>        <int>   <dbl>
#> 1 0             2212  65.6  
#> 2 >= 1          1159  34.4  
#> 3 >= 2           357  10.6  
#> 4 >= 3            94   2.79 
#> 5 >= 4            29   0.860
#> 6 >= 5             7   0.208
#> 7 >= Mean        357  10.6  
#> 8 >= Median     1159  34.4  
#> 
#> $sum_duration
#> # A tibble: 13 x 3
#>    bins          n percent
#>    <chr>     <int>   <dbl>
#>  1 >= 0       1159  100   
#>  2 >= 2       1087   93.8 
#>  3 >= 4        984   84.9 
#>  4 >= 6        862   74.4 
#>  5 >= 8        643   55.5 
#>  6 >= 10       556   48.0 
#>  7 >= 12       518   44.7 
#>  8 >= 14       337   29.1 
#>  9 >= 16       257   22.2 
#> 10 >= 18       207   17.9 
#> 11 >= 21        50    4.31
#> 12 >= Mean     531   45.8 
#> 13 >= Median   600   51.8 
#> 
#> $miss_summary
#> # A tibble: 1 x 27
#>   n_pat mean_n_vis median_n_vis mean_n_vis_out median_n_vis_out mean_n_vis_ed
#>   <int>      <dbl>        <int>          <dbl>            <int>         <dbl>
#> 1  1159       1.43            1           1.30                1        0.0328
#> # … with 21 more variables: median_n_vis_ed <int>, mean_n_vis_inpatient <dbl>,
#> #   median_n_vis_inpatient <int>, min_dur <dbl>, mean_dur <dbl>,
#> #   median_dur <dbl>, max_dur <dbl>, n_vis <int>, n_vis_out <int>,
#> #   n_vis_ed <int>, n_vis_inpatient <int>, n_vis_mean_w0 <dbl>,
#> #   n_vis_out_mean_w0 <dbl>, n_vis_ed_mean_w0 <dbl>,
#> #   n_vis_inpatient_mean_w0 <dbl>, dur_mean_w0 <dbl>, n_vis_median_w0 <dbl>,
#> #   n_vis_out_median_w0 <dbl>, n_vis_ed_median_w0 <dbl>,
#> #   n_vis_inpatient_median_w0 <dbl>, dur_median_w0 <dbl>
```

\# the next option is to use the simple\_correlated algorithm

### Simple Correlated Algorithm

The second algorithm draws cases at each period using a preference based
on whether or not the patient was drawn in prior period. The preference
is determined by the parameter `alpha` where `alpha=0` denotes strict
preference to drawing previously selected patients, `alpha=1` denotes
strict preference to patients not previously drawn and `alpha=0.5`
denotes equal preference.

By default alpha is set to 0.5

``` r
sim_miss_visits(example_sim_data,sim_algorithm = "simple_correlated")
#> $sum_n_miss
#> # A tibble: 8 x 3
#>   n_vis_groups     n percent
#>   <chr>        <int>   <dbl>
#> 1 0             2467   73.2 
#> 2 >= 1           904   26.8 
#> 3 >= 2           385   11.4 
#> 4 >= 3           183    5.43
#> 5 >= 4            77    2.28
#> 6 >= 5            42    1.25
#> 7 >= Mean        385   11.4 
#> 8 >= Median      904   26.8 
#> 
#> $sum_duration
#> # A tibble: 13 x 3
#>    bins          n percent
#>    <chr>     <int>   <dbl>
#>  1 >= 0        904  100   
#>  2 >= 2        839   92.8 
#>  3 >= 4        743   82.2 
#>  4 >= 6        638   70.6 
#>  5 >= 8        475   52.5 
#>  6 >= 10       417   46.1 
#>  7 >= 12       392   43.4 
#>  8 >= 14       260   28.8 
#>  9 >= 16       200   22.1 
#> 10 >= 18       173   19.1 
#> 11 >= 21        50    5.53
#> 12 >= Mean     400   44.2 
#> 13 >= Median   475   52.5 
#> 
#> $miss_summary
#> # A tibble: 1 x 27
#>   n_pat mean_n_vis median_n_vis mean_n_vis_out median_n_vis_out mean_n_vis_ed
#>   <int>      <dbl>        <dbl>          <dbl>            <dbl>         <dbl>
#> 1   904       1.82            1           1.68                1        0.0354
#> # … with 21 more variables: median_n_vis_ed <dbl>, mean_n_vis_inpatient <dbl>,
#> #   median_n_vis_inpatient <dbl>, min_dur <dbl>, mean_dur <dbl>,
#> #   median_dur <dbl>, max_dur <dbl>, n_vis <int>, n_vis_out <int>,
#> #   n_vis_ed <int>, n_vis_inpatient <int>, n_vis_mean_w0 <dbl>,
#> #   n_vis_out_mean_w0 <dbl>, n_vis_ed_mean_w0 <dbl>,
#> #   n_vis_inpatient_mean_w0 <dbl>, dur_mean_w0 <dbl>, n_vis_median_w0 <dbl>,
#> #   n_vis_out_median_w0 <dbl>, n_vis_ed_median_w0 <dbl>,
#> #   n_vis_inpatient_median_w0 <dbl>, dur_median_w0 <dbl>
```

You can change the value of alpha using the `sim_ctrl()` function. Here
we compare the number of patients estimated to have been missed using
different values for alpha:

``` r
res0 <-  sim_miss_visits(example_sim_data,sim_algorithm = "simple_correlated",sim_ctrl = sim_ctrl(alpha = 0))
res1 <-  sim_miss_visits(example_sim_data,sim_algorithm = "simple_correlated",sim_ctrl = sim_ctrl(alpha = 1))

res0$miss_summary$n_pat
#> [1] 635
res1$miss_summary$n_pat
#> [1] 1603
```

Notice that:

-   `alpha=0` minimizes the number of individuals with a delay (i.e.,
    under the curve)

-   `alpha=1` maximizes the number of individuals with a delay (i.e.,
    under the curve)

### Generalized Algorithm
