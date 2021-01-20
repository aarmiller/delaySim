
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
#> 1 0             2201  65.3  
#> 2 >= 1          1170  34.7  
#> 3 >= 2           338  10.0  
#> 4 >= 3           108   3.20 
#> 5 >= 4            24   0.712
#> 6 >= 5             9   0.267
#> 7 >= Mean        338  10.0  
#> 8 >= Median     1170  34.7  
#> 
#> $sum_duration
#> # A tibble: 13 x 3
#>    bins          n percent
#>    <chr>     <int>   <dbl>
#>  1 >= 0       1170  100   
#>  2 >= 2       1093   93.4 
#>  3 >= 4        979   83.7 
#>  4 >= 6        862   73.7 
#>  5 >= 8        655   56.0 
#>  6 >= 10       574   49.1 
#>  7 >= 12       535   45.7 
#>  8 >= 14       350   29.9 
#>  9 >= 16       256   21.9 
#> 10 >= 18       208   17.8 
#> 11 >= 21        50    4.27
#> 12 >= Mean     549   46.9 
#> 13 >= Median   616   52.6 
#> 
#> $miss_summary
#> # A tibble: 1 x 27
#>   n_pat mean_n_vis median_n_vis mean_n_vis_out median_n_vis_out mean_n_vis_ed
#>   <int>      <dbl>        <dbl>          <dbl>            <dbl>         <dbl>
#> 1  1170       1.41            1           1.30                1        0.0299
#> # … with 21 more variables: median_n_vis_ed <dbl>, mean_n_vis_inpatient <dbl>,
#> #   median_n_vis_inpatient <dbl>, min_dur <dbl>, mean_dur <dbl>,
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
#> 1 0             2471   73.3 
#> 2 >= 1           900   26.7 
#> 3 >= 2           395   11.7 
#> 4 >= 3           176    5.22
#> 5 >= 4            86    2.55
#> 6 >= 5            43    1.28
#> 7 >= Mean        395   11.7 
#> 8 >= Median      900   26.7 
#> 
#> $sum_duration
#> # A tibble: 13 x 3
#>    bins          n percent
#>    <chr>     <int>   <dbl>
#>  1 >= 0        900  100   
#>  2 >= 2        836   92.9 
#>  3 >= 4        739   82.1 
#>  4 >= 6        636   70.7 
#>  5 >= 8        472   52.4 
#>  6 >= 10       415   46.1 
#>  7 >= 12       390   43.3 
#>  8 >= 14       266   29.6 
#>  9 >= 16       212   23.6 
#> 10 >= 18       185   20.6 
#> 11 >= 21        50    5.56
#> 12 >= Mean     398   44.2 
#> 13 >= Median   472   52.4 
#> 
#> $miss_summary
#> # A tibble: 1 x 27
#>   n_pat mean_n_vis median_n_vis mean_n_vis_out median_n_vis_out mean_n_vis_ed
#>   <int>      <dbl>        <dbl>          <dbl>            <dbl>         <dbl>
#> 1   900       1.83            1           1.67                1        0.0489
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
#> [1] 612
res1$miss_summary$n_pat
#> [1] 1608
```

Notice that:

-   `alpha=0` minimizes the number of individuals with a delay (i.e.,
    under the curve)

-   `alpha=1` maximizes the number of individuals with a delay (i.e.,
    under the curve)

### Generalized Algorithm

The last algorithm is the generalized algorithm which takes a weighting
function and draws visits corresponding to these weights. The weighting
function can account for the number of visits a patient had during the
diagnostic opportunity window, the number of distinct SSDs they
presented with and the number of times they were previously drawn in the
simulation.

Here we apply the general algorithm using the default
`simple_weight_sum()` function; this function simply sums the number of
current and prior visits, distinct ssd and prior draws for each patient
as the weight.

``` r
sim_miss_visits(example_sim_data,
                sim_algorithm = "general",
                sim_ctrl = sim_ctrl(weight_function = simple_weight_sum))
#> $sum_n_miss
#> # A tibble: 8 x 3
#>   n_vis_groups     n percent
#>   <chr>        <int>   <dbl>
#> 1 0             2424   71.9 
#> 2 >= 1           947   28.1 
#> 3 >= 2           365   10.8 
#> 4 >= 3           139    4.12
#> 5 >= 4            75    2.22
#> 6 >= 5            46    1.36
#> 7 >= Mean        365   10.8 
#> 8 >= Median      947   28.1 
#> 
#> $sum_duration
#> # A tibble: 13 x 3
#>    bins          n percent
#>    <chr>     <int>   <dbl>
#>  1 >= 0        947  100   
#>  2 >= 2        889   93.9 
#>  3 >= 4        803   84.8 
#>  4 >= 6        714   75.4 
#>  5 >= 8        543   57.3 
#>  6 >= 10       489   51.6 
#>  7 >= 12       462   48.8 
#>  8 >= 14       310   32.7 
#>  9 >= 16       233   24.6 
#> 10 >= 18       190   20.1 
#> 11 >= 21        50    5.28
#> 12 >= Mean     470   49.6 
#> 13 >= Median   489   51.6 
#> 
#> $miss_summary
#> # A tibble: 1 x 27
#>   n_pat mean_n_vis median_n_vis mean_n_vis_out median_n_vis_out mean_n_vis_ed
#>   <int>      <dbl>        <int>          <dbl>            <int>         <dbl>
#> 1   947       1.74            1           1.57                1        0.0486
#> # … with 21 more variables: median_n_vis_ed <int>, mean_n_vis_inpatient <dbl>,
#> #   median_n_vis_inpatient <int>, min_dur <dbl>, mean_dur <dbl>,
#> #   median_dur <dbl>, max_dur <dbl>, n_vis <int>, n_vis_out <int>,
#> #   n_vis_ed <int>, n_vis_inpatient <int>, n_vis_mean_w0 <dbl>,
#> #   n_vis_out_mean_w0 <dbl>, n_vis_ed_mean_w0 <dbl>,
#> #   n_vis_inpatient_mean_w0 <dbl>, dur_mean_w0 <dbl>, n_vis_median_w0 <dbl>,
#> #   n_vis_out_median_w0 <dbl>, n_vis_ed_median_w0 <dbl>,
#> #   n_vis_inpatient_median_w0 <dbl>, dur_median_w0 <dbl>
```
