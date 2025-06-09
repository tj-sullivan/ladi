# Data cleaning, descriptives, reliability, histograms, correlations
TJ Sullivan

# Load packages and data

``` r
library(tidyverse)
library(sjmisc)
library(haven)
library(easystats)
library(purrr)
library(psych)
library(gridExtra)
library(lme4)
library(lmerTest)
```

Load apps dataset, which contains LADI coding information and relevant
variables

``` r
apps <- read_spss("data/LGBTQ-APPS 2.24.22 KB.sav")
```

Load full dataset to use for reporting demographics, calculating
reliability, and pulling variables not included in prior apps dataset
when initially requested.

``` r
full <- read_spss("data/ESTEEM All Data Merged_updated 3.2.21_with_comorbidity_RB.sav")
```

Reduce full dataset to the n=57 people who had coded data for this
study.

``` r
id <- select(apps, ParticipantID)
full <- id |> left_join(full, by = "ParticipantID")
rm(id)
```

Merge variables not included in the apps dataset (this includes HIV risk
transmission behavior from TLFB, AUDIT, SIDAS, and SIP-AD).

``` r
# isolate variables to merge
merge <- 
  full |> 
  select(ParticipantID, Assessment, Sx_CASriskacts, Sx_totalacts_sum, SIP_sum, AUDIT_sum, SIDAS_sum) |> 
  # change Assessment variable to character to match dataset; add in missing variables for BL 
  mutate(Assessment = case_when(Assessment == 0 ~ "BL",
                                Assessment == 1 ~ "4M",
                                Assessment == 2 ~ "8M",
                                Assessment == 3 ~ "12M")) |> 
  pivot_wider(id_cols = ParticipantID,
              names_from = Assessment,
              values_from = c(-ParticipantID, -Assessment),
              names_sep = ".")
# merge
apps <- left_join(apps, merge, by = "ParticipantID") |> 
  # create a variable for whether pt was engaged in CMHT, general MH services, or medication at BL - should be 0 for everyone
  mutate(CMHT.BL = 0, 
         MHgen.BL = 0,
         Medication.BL = 0)
rm(merge)
```

# Data restructuring

First we need to re-structure the dataset to long format. One dataframe
will contain the outcome/follow-up assessment data, one will contain
alliance data. This is done because there are differing numbers of time
points - 4 assessments (BL, 4MFU, 8MFU, 12MFU) and 10 sessions.

Because we are interested in controlling for the therapy alliance in
analyses, we need to separate out the alliance data and pull the
alliance scores for the sessions where the LGBTQ-APPS was coded. This
code will do that and merge it with the main dataset before
restructuring.

``` r
apps <- mutate(apps, WAI_sum_missing = rowSums(is.na(across(WAI_P_Sum.Session1:WAI_P_Sum.Session10))),
         WAI_sum_missing = recode_factor(if_else(WAI_sum_missing > 0, 1, 0), `0` = "No missing data for any session", `1` = "Missing data for at least 1 session")) # create a variable denoting whether participant has complete WAI data or is missing at least one session

apps_alliance <-
  apps |> 
  select(ParticipantID, Early_Session, Mid_Session, WAI_P_Sum.Session1:WAI_P_Sum.Session10) |> 
  pivot_longer(cols = c(-ParticipantID, -Early_Session, -Mid_Session),
    names_to = c(".value", "Session"), 
    names_sep = "\\.Session") |> 
  # rename alliance variable 
  rename(WAI.Total = WAI_P_Sum)

early <- apps_alliance |> select(ParticipantID, Session, Early_Session, Mid_Session, WAI.Total) |> 
  # match the alliance score for the specific session that was coded 
  filter(if_else(Early_Session == 3, Session == 3, Session == 2)) |> 
  select(ParticipantID, WAI_early = WAI.Total)

mid <- apps_alliance |> select(ParticipantID, Session, Early_Session, Mid_Session, WAI.Total) |> 
  # match the alliance score for the specific session that was coded 
  filter(if_else(Mid_Session == 5, Session == 5, Session == 6)) |> 
  select(ParticipantID, WAI_mid = WAI.Total)

# join the isolated alliance scores with the dataset
apps <- left_join(apps, early, by = "ParticipantID") |> 
  left_join(mid, by = "ParticipantID") 

# create average across the two sessions
apps <- apps |> 
  sjmisc::row_means(WAI_early:WAI_mid, n = 0.5, var = "WAI_avg")

rm(early)
rm(mid)
rm(apps_alliance)

full |> 
  select(contains("WAI"))
```

    # A tibble: 228 × 0

Second, separate out outcome data to separate dataframe & pivot to long
format.

``` r
apps_outcome <- apps |> 
  # select out variables for outcome analyses
  select(ParticipantID, TherapistID, Site, CMHT.BL, CMHT.4M:CMHT.12M, MHgen.BL, MHgen.4M:MHgen.12M, Medication.BL, Medication.4M:Medication.12M, Early_Affirm, Mid_Affirm, Average_Affirm, HAMD_sum.BL:LGBIS_identaffirm_missing, Sx_CASriskacts.BL:SIDAS_sum.12M, WAI_early, WAI_mid, WAI_avg) |> 
  # pivot to long format
    pivot_longer(cols = c(-ParticipantID, -Site, -TherapistID, -Early_Affirm, -Mid_Affirm, -Average_Affirm, -HAMD_missing, -BAI_missing, -RS_missing, -IHS_missing, -SOC_missing, -LGBIS_identaffirm_missing, -WAI_early, -WAI_mid, -WAI_avg),
               names_to = c(".value", "Assessment"),
               names_sep = "\\.") |>
  # clean up assessment variable to be numeric 
  mutate(Assessment = case_when(Assessment == "BL" ~ 0,
                                Assessment == "4M" ~ 1,
                                Assessment == "8M" ~ 2,
                                Assessment == "12M" ~ 3)) |> 
  relocate(Assessment, .after = TherapistID)
```

Thus `apps_outcome` is the data frame for all analyses down below.

# Centering variables

For outcome analyses, we can have assessment be centered around the
post-treatment outcome scores. This would mean that the intercept would
reflect the mean outcome variable score post-treatment/at 4MFU. The LADI
predictor variable would then signify the association b/t affirmative
technique and post-treatment outcome scores while accounting for change
in outcomes over the 12-month follow-up period.

``` r
apps_outcome <- mutate(apps_outcome,
                       Assessment_c = Assessment - 1)
```

For all analyses, we will grand-mean center the LGBTQ-APPS variable.
We’ll also center the average alliance score since we’ll include that as
a covariate.

``` r
apps_outcome <- mutate(apps_outcome,
                       Average_Affirm_c = Average_Affirm - mean(Average_Affirm),
                       WAI_avg_c = WAI_avg - mean(WAI_avg),
                       WAI_early_c = WAI_early - mean(WAI_early),
                       WAI_mid_c = WAI_mid - mean(WAI_mid),
                       Sx_totalacts_sum_c = Sx_totalacts_sum - mean(Sx_totalacts_sum, na.rm = T) 
                       )
```

# Create quadratic polynomial term

This will allow us to estimate a non-linear change trajectory, which is
pretty common in psychotherapy outcomes. Here, we’ll consider just
quadratic because cubic is not very tenable with the small sample size
that we have - it’s likely to overfit to the data.

``` r
apps_outcome <- mutate(apps_outcome,
                       Assessment_c_quad = Assessment_c * Assessment_c)
```

# Item/scale adjustments

## SIP-AD

For all items, if \<50% of items were missing, responses were imputed
with the mean of all non-missing item. For the SIP-AD, this measure is a
count (yes/no of whether or not a specific consequence was endorsed).
For participant 3352’s 12MFU, SIP item 4 was imputed as 0.36. Instead,
we modified this value to be the median value among the non-missing
responses so we could preserve the dichotomous nature of scale responses
for the SIP-AD. For this participant, this median value for the SIP-AD
at 12MFU is 0. Thus, this item was changed and 3352’s SIP_sum score for
the 12MFU is changed from 5.36 to 5 to reflect this change.

``` r
full <- full |> 
  mutate(SIP_4 = ifelse((ParticipantID == 3352 & Assessment == 3), 0, SIP_4),
         SIP_sum = ifelse((ParticipantID == 3352 & Assessment == 3), 5, SIP_sum))

apps_outcome <- apps_outcome |> 
  mutate(SIP_sum = ifelse((ParticipantID == 3352 & Assessment == 3), 5, SIP_sum))
```

## SIDAS

Given the potential for SIDAS scores to be highly skewed due to low
prevalence of suicidality in this sample, we’ll also create a variable
(0 = no, 1 = yes) as to whether or not suicidality was endorsed at that
wave.

``` r
apps_outcome <- apps_outcome |> 
  mutate(SIDAS_yn = ifelse(SIDAS_sum > 0, 1, 0))  
```

## Missing

We noticed that there were a couple of participants that scored a 999
(SPSS missing value) instead of `NA` (R missing value) on the `SIP_sum`
and `Sx_totalacts_sum` measure. This code will recode those values as
missing.

``` r
apps_outcome <- apps_outcome |> 
    mutate(SIP_sum = ifelse(SIP_sum == 999, NA, SIP_sum),
           Sx_totalacts_sum = ifelse(Sx_totalacts_sum == 999, NA, Sx_totalacts_sum))
```

# Covariates

Recode to 0/1 structure for models.

``` r
apps_outcome <- apps_outcome |> 
  mutate(CMHT = recode_factor(as.factor(CMHT), `0` = "No or no data", `1` = "Continued"),
         MHgen = recode_factor(as.factor(MHgen), `0` = "No or no data", `1` = "Utilized"),
         Medication = recode_factor(as.factor(Medication), `0` = "No or no data", `1` = "Taken"),
         Site = recode_factor(as.factor((Site - 1)), `0` = "NYC", `1` = "Miami")
         )
```

Continued CMHT treatment:

``` r
# overall
apps_outcome |> frq(CMHT)
```

    CMHT <categorical> 
    # total N=228 valid N=228 mean=1.05 sd=0.22

    Value         |   N | Raw % | Valid % | Cum. %
    ----------------------------------------------
    No or no data | 216 | 94.74 |   94.74 |  94.74
    Continued     |  12 |  5.26 |    5.26 | 100.00
    <NA>          |   0 |  0.00 |    <NA> |   <NA>

``` r
# by time point
apps_outcome |>  group_by(Assessment) |>  frq(CMHT)
```

    CMHT <categorical> 
    # grouped by: 0
    # total N=57 valid N=57 mean=1.00 sd=0.00

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 57 |   100 |     100 |    100
    Continued     |  0 |     0 |       0 |    100
    <NA>          |  0 |     0 |    <NA> |   <NA>

    CMHT <categorical> 
    # grouped by: 1
    # total N=57 valid N=57 mean=1.02 sd=0.13

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 56 | 98.25 |   98.25 |  98.25
    Continued     |  1 |  1.75 |    1.75 | 100.00
    <NA>          |  0 |  0.00 |    <NA> |   <NA>

    CMHT <categorical> 
    # grouped by: 2
    # total N=57 valid N=57 mean=1.09 sd=0.29

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 52 | 91.23 |   91.23 |  91.23
    Continued     |  5 |  8.77 |    8.77 | 100.00
    <NA>          |  0 |  0.00 |    <NA> |   <NA>

    CMHT <categorical> 
    # grouped by: 3
    # total N=57 valid N=57 mean=1.11 sd=0.31

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 51 | 89.47 |   89.47 |  89.47
    Continued     |  6 | 10.53 |   10.53 | 100.00
    <NA>          |  0 |  0.00 |    <NA> |   <NA>

MH treatment elsewhere:

``` r
# overall
apps_outcome |> frq(MHgen)
```

    MHgen <categorical> 
    # total N=228 valid N=228 mean=1.14 sd=0.34

    Value         |   N | Raw % | Valid % | Cum. %
    ----------------------------------------------
    No or no data | 197 | 86.40 |   86.40 |  86.40
    Utilized      |  31 | 13.60 |   13.60 | 100.00
    <NA>          |   0 |  0.00 |    <NA> |   <NA>

``` r
# by time point
apps_outcome |>  group_by(Assessment) |>  frq(MHgen)
```

    MHgen <categorical> 
    # grouped by: 0
    # total N=57 valid N=57 mean=1.00 sd=0.00

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 57 |   100 |     100 |    100
    Utilized      |  0 |     0 |       0 |    100
    <NA>          |  0 |     0 |    <NA> |   <NA>

    MHgen <categorical> 
    # grouped by: 1
    # total N=57 valid N=57 mean=1.11 sd=0.31

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 51 | 89.47 |   89.47 |  89.47
    Utilized      |  6 | 10.53 |   10.53 | 100.00
    <NA>          |  0 |  0.00 |    <NA> |   <NA>

    MHgen <categorical> 
    # grouped by: 2
    # total N=57 valid N=57 mean=1.21 sd=0.41

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 45 | 78.95 |   78.95 |  78.95
    Utilized      | 12 | 21.05 |   21.05 | 100.00
    <NA>          |  0 |  0.00 |    <NA> |   <NA>

    MHgen <categorical> 
    # grouped by: 3
    # total N=57 valid N=57 mean=1.23 sd=0.42

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 44 | 77.19 |   77.19 |  77.19
    Utilized      | 13 | 22.81 |   22.81 | 100.00
    <NA>          |  0 |  0.00 |    <NA> |   <NA>

Medication:

``` r
# overall
apps_outcome |> frq(Medication)
```

    Medication <categorical> 
    # total N=228 valid N=228 mean=1.12 sd=0.32

    Value         |   N | Raw % | Valid % | Cum. %
    ----------------------------------------------
    No or no data | 201 | 88.16 |   88.16 |  88.16
    Taken         |  27 | 11.84 |   11.84 | 100.00
    <NA>          |   0 |  0.00 |    <NA> |   <NA>

``` r
# by time point
apps_outcome |>  group_by(Assessment) |>  frq(Medication)
```

    Medication <categorical> 
    # grouped by: 0
    # total N=57 valid N=57 mean=1.00 sd=0.00

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 57 |   100 |     100 |    100
    Taken         |  0 |     0 |       0 |    100
    <NA>          |  0 |     0 |    <NA> |   <NA>

    Medication <categorical> 
    # grouped by: 1
    # total N=57 valid N=57 mean=1.14 sd=0.35

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 49 | 85.96 |   85.96 |  85.96
    Taken         |  8 | 14.04 |   14.04 | 100.00
    <NA>          |  0 |  0.00 |    <NA> |   <NA>

    Medication <categorical> 
    # grouped by: 2
    # total N=57 valid N=57 mean=1.18 sd=0.38

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 47 | 82.46 |   82.46 |  82.46
    Taken         | 10 | 17.54 |   17.54 | 100.00
    <NA>          |  0 |  0.00 |    <NA> |   <NA>

    Medication <categorical> 
    # grouped by: 3
    # total N=57 valid N=57 mean=1.16 sd=0.37

    Value         |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
    No or no data | 48 | 84.21 |   84.21 |  84.21
    Taken         |  9 | 15.79 |   15.79 | 100.00
    <NA>          |  0 |  0.00 |    <NA> |   <NA>

# Descriptives + missing data

Here is code that will pull descriptive information (including % missing
data) for all of the variables of interest - this averages across all
time points.

``` r
apps_outcome |>  
  select(Sx_CASriskacts, Sx_totalacts_sum, AUDIT_sum, SIP_sum, HAMD_sum, BAI_sum, SIDAS_sum, SIDAS_yn, IHS_mean, RS_mean, SOC_concealment_mean, LGBIS_identaffirm, Average_Affirm, WAI_avg) |>  
  summarize(across(everything(), list(mean = ~ mean(.x, na.rm = T),
                                      median = ~ median(.x, na.rm = T),
                                      sd = ~ sd(.x, na.rm = T),
                                      min = ~ min(.x, na.rm = T),
                                      max = ~ max(.x, na.rm = T),
                                      na_n = ~ sum(is.na(.x)),
                                      na_perc = ~ (sum(is.na(.x) / n_distinct(apps_outcome))) * 100),
                   .names = "{.col}-{.fn}")
            ) |>  
  pivot_longer(cols = everything(),
               names_to = c("ColNames", ".value"), 
               names_sep = "-",
               ) |>  
  mutate(max = round(max, digits = 2)) |>  
  unite(range, c("min", "max"), sep = "-")
```

    # A tibble: 14 × 7
       ColNames               mean median     sd range    na_n na_perc
       <chr>                 <dbl>  <dbl>  <dbl> <chr>   <int>   <dbl>
     1 Sx_CASriskacts        4.47    2     6.41  0-45       11    4.82
     2 Sx_totalacts_sum     28.6    23    27.6   0-206      11    4.82
     3 AUDIT_sum             8.16    7     6.28  0-29       12    5.26
     4 SIP_sum               3.19    2     3.68  0-15       13    5.70
     5 HAMD_sum             12.3    12     6.66  0-32       11    4.82
     6 BAI_sum              18.8    18    11.4   0-55       11    4.82
     7 SIDAS_sum             2.48    0     5.70  0-39       11    4.82
     8 SIDAS_yn              0.300   0     0.459 0-1        11    4.82
     9 IHS_mean              1.66    1.44  0.628 1-3.56     16    7.02
    10 RS_mean              12.2    11.2   7.28  1-36       12    5.26
    11 SOC_concealment_mean  1.79    1.8   0.687 1-3.4      12    5.26
    12 LGBIS_identaffirm     4.56    4.67  1.28  1-6        13    5.70
    13 Average_Affirm        1.38    1.25  1.20  0-4.75      0    0   
    14 WAI_avg              70.4    73    10.9   31.5-84     0    0   

Then do this by the specific time point, excluding `Average_Affirm` and
`WAI_avg` since these variables did not differ by time point:

``` r
apps_outcome |>  
  select(Assessment, Sx_CASriskacts, Sx_totalacts_sum, AUDIT_sum, SIP_sum, HAMD_sum, BAI_sum, SIDAS_sum, SIDAS_yn, IHS_mean, RS_mean, SOC_concealment_mean, LGBIS_identaffirm) |>  
  group_by(Assessment) |> 
  summarize(across(everything(), list(mean = ~ mean(.x, na.rm = T),
                                      median = ~ median(.x, na.rm = T),
                                      sd = ~ sd(.x, na.rm = T),
                                      min = ~ min(.x, na.rm = T),
                                      max = ~ max(.x, na.rm = T),
                                      na_n = ~ sum(is.na(.x)),
                                      na_perc = ~ (sum(is.na(.x) / n_distinct(apps_outcome))) * 100),
                   .names = "{.col}-{.fn}")
            ) %>% 
  pivot_longer(cols = -Assessment,
               names_to = c("ColNames", ".value"), 
               names_sep = "-",
               ) %>% 
  mutate(max = round(max, digits = 2),
         min = round(min, digits = 2)) %>% 
  unite(range, c("min", "max"), sep = "-") |> 
  arrange(ColNames)
```

    # A tibble: 48 × 8
       Assessment ColNames   mean median    sd range  na_n na_perc
            <dbl> <chr>     <dbl>  <dbl> <dbl> <chr> <int>   <dbl>
     1          0 AUDIT_sum  9.49      8  7.27 0-29      0    0   
     2          1 AUDIT_sum  8.13      7  6.13 0-23      4    1.75
     3          2 AUDIT_sum  8.23      8  5.77 0-24      4    1.75
     4          3 AUDIT_sum  6.70      6  5.57 0-22      4    1.75
     5          0 BAI_sum   23.4      24 11.1  2-47      0    0   
     6          1 BAI_sum   17.7      18 10.1  0-45      3    1.32
     7          2 BAI_sum   16.7      15 10.5  0-46      4    1.75
     8          3 BAI_sum   17.1      14 12.6  0-55      4    1.75
     9          0 HAMD_sum  15.6      16  6.28 3-32      0    0   
    10          1 HAMD_sum  10.8       9  6.60 1-28      3    1.32
    # ℹ 38 more rows

# Scale reliability

For reporting purposes, we’re going to pull Cronbach’s alpha overall and
for each time point. These functions will allow us to do that easily for
each time point.

``` r
# Define a function to just get the raw Cronbach's alpha from the psych package for each measure (aggregated across time point)
get_raw_alpha <- function(data, measure_cols) {
  alpha <- data %>%
    select({{measure_cols}}) %>%
    alpha(warnings = FALSE)
  
  raw_alpha <- alpha$total$raw_alpha
  return(raw_alpha)
}

# Define a function to calculate Cronbach's alpha for each assessment time
calculate_alpha_map <- function(data, assessment_col, measure_cols) {
  # Get unique assessment times
  assessment_times <- unique(data[[assessment_col]])
  
  # Define a function to calculate alpha for a specific assessment time
  alpha_for_time <- function(assessment_time) {
    data_filtered <- data %>% filter(!!sym(assessment_col) == assessment_time)
    alpha_result <- data_filtered %>% select({{measure_cols}}) %>% alpha(warnings = F)
    return(alpha_result$total$raw_alpha)
  }
  
  # Use map to apply the alpha function to each assessment time
  alpha_results <- map(assessment_times, alpha_for_time)
  names(alpha_results) <- paste("Assessment", assessment_times)
  
  return(alpha_results)
}
```

## AUDIT

``` r
get_raw_alpha(full, AUDIT_1:AUDIT_9)
```

    [1] 0.8568083

``` r
calculate_alpha_map(full, "Assessment", AUDIT_1:AUDIT_9) 
```

    $`Assessment 0`
    [1] 0.8660913

    $`Assessment 1`
    [1] 0.8475727

    $`Assessment 2`
    [1] 0.8367354

    $`Assessment 3`
    [1] 0.8510599

## SIP-AD

Note this will give you the KR-20 coefficient when the `alpha` function
is run with dichotomous items.

``` r
get_raw_alpha(full, SIP_1:SIP_15)
```

    [1] 0.8897054

``` r
calculate_alpha_map(full, "Assessment", SIP_1:SIP_15) 
```

    $`Assessment 0`
    [1] 0.9060917

    $`Assessment 1`
    [1] 0.8807223

    $`Assessment 2`
    [1] 0.8937622

    $`Assessment 3`
    [1] 0.8452679

## HAMD

``` r
get_raw_alpha(full, HAMD_1:HAMD_17)
```

    [1] 0.7993866

``` r
calculate_alpha_map(full, "Assessment", HAMD_1:HAMD_17) 
```

    $`Assessment 0`
    [1] 0.7376662

    $`Assessment 1`
    [1] 0.8107525

    $`Assessment 2`
    [1] 0.7873542

    $`Assessment 3`
    [1] 0.8127561

## BAI

``` r
get_raw_alpha(full, BAI_1:BAI_21)
```

    [1] 0.9148405

``` r
calculate_alpha_map(full, "Assessment", BAI_1:BAI_21) 
```

    $`Assessment 0`
    [1] 0.8930645

    $`Assessment 1`
    [1] 0.8945386

    $`Assessment 2`
    [1] 0.9039994

    $`Assessment 3`
    [1] 0.9401277

## SIDAS (continuous)

Note for the SIDAS the missing values were not coded as `NA` for item 2
so those need to be recoded.

``` r
full <- full |> 
    mutate(SIDAS_2r = ifelse(SIDAS_2r == 999, NA, SIDAS_2r))
```

Then calculate reliability:

``` r
get_raw_alpha(full, c(SIDAS_1, SIDAS_2r, SIDAS_3, SIDAS_4, SIDAS_5))
```

    [1] 0.8967987

``` r
calculate_alpha_map(full, "Assessment", c(SIDAS_2r, SIDAS_3, SIDAS_4, SIDAS_5))
```

    $`Assessment 0`
    [1] 0.7785711

    $`Assessment 1`
    [1] 0.8251342

    $`Assessment 2`
    [1] 0.9417549

    $`Assessment 3`
    [1] 0.8144857

## IHS

``` r
get_raw_alpha(full, IHS_1:IHS_9)
```

    [1] 0.8817837

``` r
calculate_alpha_map(full, "Assessment", IHS_1:IHS_9) 
```

    $`Assessment 0`
    [1] 0.8806606

    $`Assessment 1`
    [1] 0.8595452

    $`Assessment 2`
    [1] 0.8890164

    $`Assessment 3`
    [1] 0.8946053

## RS

``` r
get_raw_alpha(full, RS_1:RS_14)
```

    [1] 0.9173433

``` r
calculate_alpha_map(full, "Assessment", RS_1:RS_14) 
```

    $`Assessment 0`
    [1] 0.914077

    $`Assessment 1`
    [1] 0.8867508

    $`Assessment 2`
    [1] 0.9313299

    $`Assessment 3`
    [1] 0.9207489

## Concealment

``` r
get_raw_alpha(full, c(SOC_FamilyR, SOC_GLBFriendsR, SOC_Str8FriendsR, SOC_CoworkersR, SOC_HCPR))
```

    [1] 0.7853371

``` r
calculate_alpha_map(full, "Assessment", c(SOC_FamilyR, SOC_GLBFriendsR, SOC_Str8FriendsR, SOC_CoworkersR, SOC_HCPR)) 
```

    $`Assessment 0`
    [1] 0.7929316

    $`Assessment 1`
    [1] 0.8101648

    $`Assessment 2`
    [1] 0.8328684

    $`Assessment 3`
    [1] 0.8122468

## Identity affirmation

``` r
get_raw_alpha(full, c(LGBIS_6, LGBIS_13, LGBIS_26))
```

    [1] 0.9224737

``` r
calculate_alpha_map(full, "Assessment", c(LGBIS_6, LGBIS_13, LGBIS_26)) 
```

    $`Assessment 0`
    [1] 0.954971

    $`Assessment 1`
    [1] 0.9549971

    $`Assessment 2`
    [1] 0.8878999

    $`Assessment 3`
    [1] 0.8804268

## Working alliance

Because different sessions were coded, some wrangling is needed to get
the item-level alliance data for the session where the LADI was coded.
This first chunk will do that and store them in separate dataframes to
calculate reliability.

``` r
### code to pull in the item-by-item data for the alliance
alliance_items <- read_spss("data/WIDE_Therapy Post-Session Surveys for ESTEEM study 2021.10.21.sav")

# create a data frame that denotes which session the LADI was coded for each participant 
apps_alliance_rel <- 
  apps |> 
  select(ParticipantID, Early_Session, Mid_Session)

####### EARLY #######

# select only relevant columns from the alliance data - items 4 & 10 are reverse scored so this code replaces the original with the rescored
early_alliance_items <- alliance_items |> 
  select(ParticipantID, contains("P.2"), contains("P.3"), 
         -WAI4_P.2, -WAI4_P.3, -WAI10_P.2, -WAI10_P.3, 
         WAI4_P_R.2, WAI4_P_R.3, WAI10_P_R.2, WAI10_P_R.3) |> 
  select(ParticipantID, contains("WAI"), 
         WAI4_P_R.2, WAI4_P_R.3, WAI10_P_R.2, WAI10_P_R.3) |> 
  relocate(WAI4_P_R.2, .after = WAI3_P.2) |> 
  relocate(WAI4_P_R.3, .after = WAI3_P.3) |> 
  relocate(WAI10_P_R.2, .after = WAI9_P.2) |> 
  relocate(WAI10_P_R.3, .after = WAI9_P.3) 

# merge the two datasets together
apps_alliance_early <- apps_alliance_rel |> 
  left_join(early_alliance_items, by = "ParticipantID")

# pull the early session 2s
early_s2 <- apps_alliance_early |> 
  filter(Early_Session == 2) |> 
  select(ParticipantID, Early_Session, contains(".2")) |> 
  rename_with(~ str_replace_all(., c("_P.2" = "_early", "_P_R.2" = "_early")))

# pull the early session 3s
early_s3 <- apps_alliance_early |> 
  filter(Early_Session == 3) |> 
  select(ParticipantID, Early_Session, contains(".3")) |> 
  rename_with(~ str_replace_all(., c("_P.3" = "_early", "_P_R.3" = "_early")))

# merge early together
early_alliance <- rbind(early_s2, early_s3)

####### MIDDLE #######

# create a data frame that denotes which session the LADI was coded for each participant 
apps_alliance_rel <- 
  apps |> 
  select(ParticipantID, Early_Session, Mid_Session)

# select only relevant columns from the alliance data - items 4 & 10 are reverse scored so this code replaces the original with the rescored
mid_alliance_items <- alliance_items |> 
  select(ParticipantID, contains("P.5"), contains("P.6"), 
         -WAI4_P.5, -WAI4_P.6, -WAI10_P.5, -WAI10_P.6, 
         WAI4_P_R.5, WAI4_P_R.6, WAI10_P_R.5, WAI10_P_R.6) |> 
  select(ParticipantID, contains("WAI"), 
         WAI4_P_R.5, WAI4_P_R.6, WAI10_P_R.5, WAI10_P_R.6) |> 
  relocate(WAI4_P_R.5, .after = WAI3_P.5) |> 
  relocate(WAI4_P_R.6, .after = WAI3_P.6) |> 
  relocate(WAI10_P_R.5, .after = WAI9_P.5) |> 
  relocate(WAI10_P_R.6, .after = WAI9_P.6) 

# merge the two datasets together
apps_alliance_mid <- apps_alliance_rel |> 
  left_join(mid_alliance_items, by = "ParticipantID")

# pull the mid session 5s
mid_s2 <- apps_alliance_mid |> 
  filter(Mid_Session == 5) |> 
  select(ParticipantID, Mid_Session, contains(".5")) |> 
  rename_with(~ str_replace_all(., c("_P.5" = "_mid", "_P_R.5" = "_mid")))

# pull the mid session 6s
mid_s3 <- apps_alliance_mid |> 
  filter(Mid_Session == 6) |> 
  select(ParticipantID, Mid_Session, contains(".6")) |> 
  rename_with(~ str_replace_all(., c("_P.6" = "_mid", "_P_R.6" = "_mid")))

# merge early together
mid_alliance <- rbind(mid_s2, mid_s3)
```

Reliability for early alliance scores:

``` r
# calculate reliability for early 
early_alliance |> 
  select(-ParticipantID, -Early_Session) |> 
  alpha()
```


    Reliability analysis   
    Call: alpha(x = select(early_alliance, -ParticipantID, -Early_Session))

      raw_alpha std.alpha G6(smc) average_r S/N   ase mean   sd median_r
          0.91      0.93    0.96      0.53  14 0.018  5.8 0.97     0.67

        95% confidence boundaries 
             lower alpha upper
    Feldt     0.88  0.91  0.94
    Duhachek  0.88  0.91  0.95

     Reliability if an item is dropped:
                raw_alpha std.alpha G6(smc) average_r S/N alpha se var.r med.r
    WAI1_early       0.90      0.92    0.95      0.51  12    0.021 0.079  0.66
    WAI2_early       0.90      0.93    0.96      0.53  13    0.020 0.084  0.67
    WAI3_early       0.90      0.92    0.95      0.52  12    0.020 0.082  0.67
    WAI4_early       0.93      0.94    0.96      0.59  16    0.015 0.072  0.70
    WAI5_early       0.90      0.92    0.95      0.51  11    0.021 0.081  0.64
    WAI6_early       0.90      0.92    0.95      0.52  12    0.021 0.085  0.67
    WAI7_early       0.90      0.92    0.95      0.51  11    0.021 0.079  0.64
    WAI8_early       0.90      0.92    0.95      0.52  12    0.020 0.081  0.67
    WAI9_early       0.91      0.93    0.95      0.54  13    0.019 0.081  0.68
    WAI10_early      0.94      0.95    0.97      0.62  18    0.013 0.045  0.70
    WAI11_early      0.90      0.92    0.95      0.52  12    0.021 0.082  0.67
    WAI12_early      0.90      0.92    0.95      0.51  11    0.022 0.080  0.64

     Item statistics 
                 n raw.r std.r r.cor r.drop mean  sd
    WAI1_early  56  0.87  0.89  0.89   0.84  6.0 1.1
    WAI2_early  56  0.74  0.75  0.73   0.68  5.8 1.5
    WAI3_early  56  0.79  0.81  0.80   0.74  6.0 1.2
    WAI4_early  56  0.46  0.41  0.34   0.33  5.5 1.8
    WAI5_early  56  0.88  0.89  0.88   0.85  5.8 1.2
    WAI6_early  56  0.86  0.86  0.86   0.82  6.0 1.2
    WAI7_early  56  0.91  0.92  0.93   0.89  5.9 1.2
    WAI8_early  56  0.82  0.84  0.84   0.78  6.0 1.1
    WAI9_early  56  0.69  0.71  0.70   0.62  5.8 1.4
    WAI10_early 56  0.29  0.23  0.15   0.14  5.3 1.8
    WAI11_early 56  0.84  0.86  0.86   0.80  5.8 1.3
    WAI12_early 56  0.90  0.91  0.91   0.88  5.8 1.3

    Non missing response frequency for each item
                   1    2    3    4    5    6    7 miss
    WAI1_early  0.00 0.00 0.00 0.14 0.16 0.29 0.41 0.02
    WAI2_early  0.04 0.02 0.00 0.11 0.16 0.29 0.39 0.02
    WAI3_early  0.00 0.00 0.04 0.12 0.16 0.16 0.52 0.02
    WAI4_early  0.05 0.05 0.07 0.09 0.02 0.32 0.39 0.02
    WAI5_early  0.00 0.02 0.02 0.12 0.12 0.38 0.34 0.02
    WAI6_early  0.02 0.00 0.00 0.09 0.20 0.21 0.48 0.02
    WAI7_early  0.00 0.02 0.02 0.11 0.18 0.29 0.39 0.02
    WAI8_early  0.00 0.00 0.04 0.07 0.18 0.25 0.46 0.02
    WAI9_early  0.02 0.00 0.05 0.07 0.20 0.23 0.43 0.02
    WAI10_early 0.07 0.04 0.04 0.11 0.14 0.32 0.29 0.02
    WAI11_early 0.00 0.02 0.04 0.12 0.18 0.23 0.41 0.02
    WAI12_early 0.02 0.00 0.02 0.12 0.18 0.32 0.34 0.02

Reliability for middle alliance scores:

``` r
# calculate reliability for middle 
mid_alliance |> 
  select(-ParticipantID, -Mid_Session) |> 
  alpha()
```


    Reliability analysis   
    Call: alpha(x = select(mid_alliance, -ParticipantID, -Mid_Session))

      raw_alpha std.alpha G6(smc) average_r S/N   ase mean sd median_r
          0.93      0.95    0.97      0.62  19 0.014  5.9  1      0.7

        95% confidence boundaries 
             lower alpha upper
    Feldt      0.9  0.93  0.96
    Duhachek   0.9  0.93  0.96

     Reliability if an item is dropped:
              raw_alpha std.alpha G6(smc) average_r S/N alpha se var.r med.r
    WAI1_mid       0.92      0.94    0.96      0.60  17   0.0165 0.066  0.69
    WAI2_mid       0.92      0.94    0.97      0.60  17   0.0169 0.064  0.70
    WAI3_mid       0.92      0.95    0.97      0.62  18   0.0160 0.067  0.71
    WAI4_mid       0.93      0.95    0.97      0.66  21   0.0141 0.066  0.73
    WAI5_mid       0.92      0.94    0.96      0.60  17   0.0170 0.065  0.70
    WAI6_mid       0.92      0.94    0.97      0.60  16   0.0168 0.066  0.69
    WAI7_mid       0.92      0.95    0.97      0.61  17   0.0162 0.067  0.70
    WAI8_mid       0.92      0.94    0.96      0.61  17   0.0164 0.065  0.70
    WAI9_mid       0.92      0.95    0.97      0.61  17   0.0164 0.068  0.71
    WAI10_mid      0.96      0.96    0.98      0.71  27   0.0084 0.023  0.73
    WAI11_mid      0.92      0.94    0.96      0.60  17   0.0165 0.063  0.69
    WAI12_mid      0.92      0.94    0.96      0.60  16   0.0174 0.063  0.69

     Item statistics 
               n raw.r std.r r.cor r.drop mean  sd
    WAI1_mid  57  0.88  0.89  0.89   0.86  6.1 1.1
    WAI2_mid  57  0.89  0.90  0.90   0.86  6.1 1.3
    WAI3_mid  57  0.80  0.82  0.82   0.76  6.0 1.2
    WAI4_mid  57  0.63  0.59  0.54   0.54  5.9 1.7
    WAI5_mid  57  0.88  0.90  0.90   0.86  5.9 1.4
    WAI6_mid  57  0.90  0.91  0.91   0.88  6.1 1.2
    WAI7_mid  57  0.82  0.84  0.84   0.79  6.0 1.2
    WAI8_mid  57  0.86  0.88  0.88   0.83  6.0 1.1
    WAI9_mid  57  0.82  0.83  0.82   0.78  5.9 1.4
    WAI10_mid 57  0.36  0.29  0.22   0.20  5.1 2.1
    WAI11_mid 57  0.87  0.89  0.89   0.84  6.0 1.1
    WAI12_mid 57  0.91  0.93  0.93   0.89  5.9 1.4

    Non missing response frequency for each item
                 1    2    3    4    5    6    7 miss
    WAI1_mid  0.00 0.02 0.02 0.04 0.14 0.37 0.42    0
    WAI2_mid  0.02 0.02 0.00 0.09 0.07 0.32 0.49    0
    WAI3_mid  0.00 0.00 0.04 0.09 0.18 0.23 0.47    0
    WAI4_mid  0.05 0.04 0.00 0.09 0.07 0.21 0.54    0
    WAI5_mid  0.02 0.04 0.02 0.04 0.18 0.30 0.42    0
    WAI6_mid  0.02 0.00 0.00 0.09 0.09 0.33 0.47    0
    WAI7_mid  0.00 0.02 0.02 0.07 0.18 0.26 0.46    0
    WAI8_mid  0.00 0.02 0.00 0.11 0.12 0.32 0.44    0
    WAI9_mid  0.04 0.00 0.02 0.04 0.21 0.26 0.44    0
    WAI10_mid 0.14 0.05 0.05 0.07 0.04 0.35 0.30    0
    WAI11_mid 0.00 0.00 0.04 0.09 0.19 0.25 0.44    0
    WAI12_mid 0.04 0.00 0.02 0.11 0.14 0.25 0.46    0

Remove dataframes no longer needed:

``` r
rm(apps_alliance_early)
rm(apps_alliance_mid)
rm(apps_alliance_rel)
rm(early_alliance)
rm(early_alliance_items)
rm(early_s2)
rm(early_s3)
rm(mid_alliance)
rm(mid_alliance_items)
rm(mid_s2)
rm(mid_s3)
rm(alliance_items)
```

# Overall histograms

``` r
plot_outcomes_sum <- function(data, outcome_vars) {
  plots <- map(outcome_vars, function(var) {
    var_sym <- sym(var)
    data %>%
      ggplot(aes(x = !!var_sym)) +
      geom_histogram(binwidth = 1)
  })
  
  return(plots)
}

outcome_vars_sum <- c("Sx_CASriskacts", "AUDIT_sum", "SIP_sum", "HAMD_sum", "BAI_sum", "SIDAS_sum", "SIDAS_yn")
plot_outcomes_sum(apps_outcome, outcome_vars_sum)
```

    [[1]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-1.png)


    [[2]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-2.png)


    [[3]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-3.png)


    [[4]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-4.png)


    [[5]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-5.png)


    [[6]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-6.png)


    [[7]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-7.png)

``` r
plot_outcomes_mean <- function(data, outcome_vars) {
  plots <- map(outcome_vars, function(var) {
    var_sym <- sym(var)
    data %>%
      ggplot(aes(x = !!var_sym)) +
      geom_histogram(binwidth = 0.3)
  })
  
  return(plots)
}

outcome_vars_mean <- c("IHS_mean", "RS_mean", "SOC_concealment_mean", "LGBIS_identaffirm")
plot_outcomes_mean(apps_outcome, outcome_vars_mean)
```

    [[1]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-8.png)


    [[2]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-9.png)


    [[3]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-10.png)


    [[4]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-28-11.png)

# Change trajectories

## Individual plots

``` r
individ_trajs <- function(data, outcome_vars) {
  plots <- map(outcome_vars, function(var) {
    var_sym <- sym(var)
    data |> 
      ggplot(aes(x = Assessment, y = !!var_sym)) +
      geom_point() +
      stat_smooth(method = "loess", se = FALSE) +
      facet_wrap(~ ParticipantID) 
  })
  
  return(plots)
}

outcome_vars <- c("Sx_CASriskacts", "AUDIT_sum", "SIP_sum", "HAMD_sum", "BAI_sum", "SIDAS_sum", "IHS_mean", "RS_mean", "SOC_concealment_mean")
individ_trajs(apps_outcome, outcome_vars)
```

    [[1]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-29-1.png)


    [[2]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-29-2.png)


    [[3]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-29-3.png)


    [[4]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-29-4.png)


    [[5]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-29-5.png)


    [[6]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-29-6.png)


    [[7]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-29-7.png)


    [[8]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-29-8.png)


    [[9]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-29-9.png)

### Spaghetti plots

``` r
spaghetti <- function(data, outcome_vars) {
  plots <- map(outcome_vars, function(var) {
    var_sym <- sym(var)
    data |> 
      ggplot(aes(x = Assessment, y = !!var_sym)) +
      stat_smooth(aes(group = ParticipantID),
              method = "loess", se = F, size = 1/16) +
      stat_smooth(method = "loess", se = F, size = 2, color = "purple")
  })
  
  return(plots)
}

spaghetti(apps_outcome, outcome_vars)
```

    [[1]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-30-1.png)


    [[2]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-30-2.png)


    [[3]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-30-3.png)


    [[4]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-30-4.png)


    [[5]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-30-5.png)


    [[6]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-30-6.png)


    [[7]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-30-7.png)


    [[8]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-30-8.png)


    [[9]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-30-9.png)

# Scatterplots

Overall scores:

``` r
#apps_outcome |> 
#  ggplot(aes(x = Average_Affirm, y = HAMD_sum)) + 
#  geom_point() + 
#  geom_jitter() +
#  stat_smooth(method = "loess")

scatter <- function(data, outcome_vars) {
  plots <- map(outcome_vars, function(var) {
    var_sym <- sym(var)
    data |> 
      ggplot(aes(x = Average_Affirm, y = !!var_sym)) +
      geom_point() +
      geom_jitter() +
      stat_smooth(method = "loess")
  })
  
  return(plots)
}

scatter(apps_outcome, outcome_vars)
```

    [[1]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-31-1.png)


    [[2]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-31-2.png)


    [[3]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-31-3.png)


    [[4]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-31-4.png)


    [[5]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-31-5.png)


    [[6]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-31-6.png)


    [[7]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-31-7.png)


    [[8]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-31-8.png)


    [[9]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-31-9.png)

Just post-treatment scores to make this clearer. Note that we wrote over
the function above.

``` r
scatter <- function(data, outcome_vars) {
  plots <- map(outcome_vars, function(var) {
    var_sym <- sym(var)
    data |> 
      filter(Assessment == 1) |> 
      ggplot(aes(x = Average_Affirm, y = !!var_sym)) +
      geom_point() +
      geom_jitter() +
      stat_smooth(method = "loess")
  })
  
  return(plots)
}

scatter(apps_outcome, outcome_vars)
```

    [[1]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-32-1.png)


    [[2]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-32-2.png)


    [[3]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-32-3.png)


    [[4]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-32-4.png)


    [[5]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-32-5.png)


    [[6]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-32-6.png)


    [[7]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-32-7.png)


    [[8]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-32-8.png)


    [[9]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-32-9.png)

Inspection of the bivariate scatterplot shows that there are a couple of
points at the far range of the LADI that might be influencing the
bivariate relationship. Let’s identify what those scores are:

``` r
apps_outcome |> 
  frq(Average_Affirm)
```

    Average of early and mid-sessions LGBTQ-APPS Score (Average_Affirm) <numeric> 
    # total N=228 valid N=228 mean=1.38 sd=1.20

    Value |  N | Raw % | Valid % | Cum. %
    -------------------------------------
     0.00 | 44 | 19.30 |   19.30 |  19.30
     0.25 | 12 |  5.26 |    5.26 |  24.56
     0.50 | 16 |  7.02 |    7.02 |  31.58
     0.75 | 24 | 10.53 |   10.53 |  42.11
     1.00 | 16 |  7.02 |    7.02 |  49.12
     1.25 | 20 |  8.77 |    8.77 |  57.89
     1.50 | 12 |  5.26 |    5.26 |  63.16
     1.75 | 24 | 10.53 |   10.53 |  73.68
     2.25 | 16 |  7.02 |    7.02 |  80.70
     2.50 |  4 |  1.75 |    1.75 |  82.46
     2.75 |  4 |  1.75 |    1.75 |  84.21
     3.00 |  8 |  3.51 |    3.51 |  87.72
     3.25 |  8 |  3.51 |    3.51 |  91.23
     3.50 | 12 |  5.26 |    5.26 |  96.49
     4.00 |  4 |  1.75 |    1.75 |  98.25
     4.75 |  4 |  1.75 |    1.75 | 100.00
     <NA> |  0 |  0.00 |    <NA> |   <NA>

Ok we see here that those are likely the one participant whose therapist
scored scored a 4.75 on the LADI. Let’s exclude those to see how
sensitive the plots are to those two points. Note that we wrote over the
function above to make this easier.

``` r
scatter <- function(data, outcome_vars) {
  plots <- map(outcome_vars, function(var) {
    var_sym <- sym(var)
    data |> 
      filter(Assessment == 1) |> 
      filter(Average_Affirm != 4.75) |> 
      ggplot(aes(x = Average_Affirm, y = !!var_sym)) +
      geom_point() +
      geom_jitter() +
      stat_smooth(method = "loess")
  })
  
  return(plots)
}

scatter(apps_outcome, outcome_vars)
```

    [[1]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-34-1.png)


    [[2]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-34-2.png)


    [[3]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-34-3.png)


    [[4]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-34-4.png)


    [[5]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-34-5.png)


    [[6]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-34-6.png)


    [[7]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-34-7.png)


    [[8]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-34-8.png)


    [[9]]

![](analyses-basic_files/figure-commonmark/unnamed-chunk-34-9.png)

By and large the scatterplots are similar with and without the excluded
folks. However, the HAMD, IHS, and SOC are all measures that seem to
change. Here is the before and after for each to showcase this.

For HAMD:

``` r
hamd.1 <- apps_outcome |> 
      filter(Assessment == 1) |> 
      ggplot(aes(x = Average_Affirm, y = HAMD_sum)) +
      geom_point() +
      geom_jitter() +
      stat_smooth(method = "loess")

hamd.2 <- apps_outcome |> 
      filter(Assessment == 1) |> 
      filter(Average_Affirm != 4.75) |> 
      ggplot(aes(x = Average_Affirm, y = HAMD_sum)) +
      geom_point() +
      geom_jitter() +
      stat_smooth(method = "loess")

grid.arrange(hamd.1, hamd.2, ncol = 2)
```

![](analyses-basic_files/figure-commonmark/unnamed-chunk-35-1.png)

For IHS:

``` r
ihs.1 <- apps_outcome |> 
      filter(Assessment == 1) |> 
      ggplot(aes(x = Average_Affirm, y = IHS_mean)) +
      geom_point() +
      geom_jitter() +
      stat_smooth(method = "loess")

ihs.2 <- apps_outcome |> 
      filter(Assessment == 1) |> 
      filter(Average_Affirm != 4.75) |> 
      ggplot(aes(x = Average_Affirm, y = IHS_mean)) +
      geom_point() +
      geom_jitter() +
      stat_smooth(method = "loess")

grid.arrange(ihs.1, ihs.2, ncol = 2)
```

![](analyses-basic_files/figure-commonmark/unnamed-chunk-36-1.png)

For SOC:

``` r
soc.1 <- apps_outcome |> 
      filter(Assessment == 1) |> 
      ggplot(aes(x = Average_Affirm, y = SOC_concealment_mean)) +
      geom_point() +
      geom_jitter() +
      stat_smooth(method = "loess")

soc.2 <- apps_outcome |> 
      filter(Assessment == 1) |> 
      filter(Average_Affirm != 4.75) |> 
      ggplot(aes(x = Average_Affirm, y = SOC_concealment_mean)) +
      geom_point() +
      geom_jitter() +
      stat_smooth(method = "loess")

grid.arrange(soc.1, soc.2, ncol = 2)
```

![](analyses-basic_files/figure-commonmark/unnamed-chunk-37-1.png)

# Correlations

We’ll pull a correlation matrix for variables across all time points.

``` r
corr_results <- correlation(data = apps_outcome, 
            select = c("Sx_CASriskacts", "Sx_totalacts_sum", "AUDIT_sum", "SIP_sum", "HAMD_sum", "BAI_sum", "SIDAS_sum", "IHS_mean", "RS_mean", "SOC_concealment_mean", "Average_Affirm", "WAI_avg")) |> 
  summary(redundant = T) 
corr_results
```

    # Correlation Matrix (pearson-method)

    Parameter            | Sx_CASriskacts | Sx_totalacts_sum | AUDIT_sum | SIP_sum
    ------------------------------------------------------------------------------
    Sx_CASriskacts       |                |          0.60*** |     0.22* |    0.20
    Sx_totalacts_sum     |        0.60*** |                  |      0.19 |    0.18
    AUDIT_sum            |          0.22* |             0.19 |           | 0.64***
    SIP_sum              |           0.20 |             0.18 |   0.64*** |        
    HAMD_sum             |           0.04 |             0.08 |      0.22 |  0.26**
    BAI_sum              |          -0.02 |             0.02 |     0.24* |  0.29**
    SIDAS_sum            |          -0.06 |             0.01 |      0.16 |    0.16
    IHS_mean             |           0.06 |             0.10 |    0.27** |    0.21
    RS_mean              |           0.04 |             0.16 |      0.18 |  0.27**
    SOC_concealment_mean |          -0.10 |             0.04 |     -0.07 |   -0.05
    Average_Affirm       |           0.13 |             0.08 |     -0.02 |   -0.12
    WAI_avg              |          -0.01 |             0.01 |      0.08 |    0.12

    Parameter            | HAMD_sum | BAI_sum | SIDAS_sum | IHS_mean | RS_mean
    --------------------------------------------------------------------------
    Sx_CASriskacts       |     0.04 |   -0.02 |     -0.06 |     0.06 |    0.04
    Sx_totalacts_sum     |     0.08 |    0.02 |      0.01 |     0.10 |    0.16
    AUDIT_sum            |     0.22 |   0.24* |      0.16 |   0.27** |    0.18
    SIP_sum              |   0.26** |  0.29** |      0.16 |     0.21 |  0.27**
    HAMD_sum             |          | 0.37*** |   0.29*** |    0.24* |  0.26**
    BAI_sum              |  0.37*** |         |      0.15 |     0.10 |   0.23*
    SIDAS_sum            |  0.29*** |    0.15 |           |     0.13 |    0.14
    IHS_mean             |    0.24* |    0.10 |      0.13 |          | 0.34***
    RS_mean              |   0.26** |   0.23* |      0.14 |  0.34*** |        
    SOC_concealment_mean |     0.02 |   -0.06 |      0.09 |  0.35*** |    0.18
    Average_Affirm       |    -0.12 |   -0.02 |     -0.15 |    -0.02 |   -0.01
    WAI_avg              |    -0.13 |   -0.09 |     -0.05 |    -0.14 |    0.09

    Parameter            | SOC_concealment_mean | Average_Affirm | WAI_avg
    ----------------------------------------------------------------------
    Sx_CASriskacts       |                -0.10 |           0.13 |   -0.01
    Sx_totalacts_sum     |                 0.04 |           0.08 |    0.01
    AUDIT_sum            |                -0.07 |          -0.02 |    0.08
    SIP_sum              |                -0.05 |          -0.12 |    0.12
    HAMD_sum             |                 0.02 |          -0.12 |   -0.13
    BAI_sum              |                -0.06 |          -0.02 |   -0.09
    SIDAS_sum            |                 0.09 |          -0.15 |   -0.05
    IHS_mean             |              0.35*** |          -0.02 |   -0.14
    RS_mean              |                 0.18 |          -0.01 |    0.09
    SOC_concealment_mean |                      |        -0.27** |  -0.23*
    Average_Affirm       |              -0.27** |                |    0.15
    WAI_avg              |               -0.23* |           0.15 |        

    p-value adjustment method: Holm (1979)

# ICCs

We’ll calculate ICCs from a random intercepts MLM. First, we’ll evaluate
nesting of timepoints within participants.

``` r
# function to calculate ICCs for all outcomes
icc_pt <- apps_outcome %>% 
  select(ParticipantID, Assessment_c, Sx_CASriskacts, Sx_totalacts_sum, AUDIT_sum, SIP_sum, HAMD_sum, BAI_sum, SIDAS_sum, IHS_mean, RS_mean, SOC_concealment_mean, LGBIS_identaffirm) %>% 
  pivot_longer(c(Sx_CASriskacts, Sx_totalacts_sum, AUDIT_sum, SIP_sum, HAMD_sum, BAI_sum, SIDAS_sum, IHS_mean, RS_mean, SOC_concealment_mean, LGBIS_identaffirm),
               names_to = "variable",
               values_to = "outcome") %>% 
  mutate(variable = as.factor(variable)) %>% 
  group_by(variable) %>% 
  group_split() %>% 
  map_dfr(.f = function(df) {
    lmer(outcome ~ 1 + (1 | ParticipantID), REML = T, data = df) %>% 
      broom.mixed::tidy(effects = "ran_pars", ddf.method = "Kenward-Roger", scales = "vcov") %>% 
      add_column(variable = unique(df$variable), .before = 1) %>% 
      group_by(variable) %>%
      summarize(pt_var = first(estimate),
                total_var = sum(estimate),
                icc = pt_var/total_var)
  })
icc_pt %>% select(variable, icc) %>% mutate_if(is.numeric, round, digits = 3)
```

    # A tibble: 11 × 2
       variable               icc
       <fct>                <dbl>
     1 AUDIT_sum            0.721
     2 BAI_sum              0.506
     3 HAMD_sum             0.236
     4 IHS_mean             0.677
     5 LGBIS_identaffirm    0.776
     6 RS_mean              0.624
     7 SIDAS_sum            0.328
     8 SIP_sum              0.63 
     9 SOC_concealment_mean 0.872
    10 Sx_CASriskacts       0.305
    11 Sx_totalacts_sum     0.232

Second, we’ll consider nesting of participants within therapists.

``` r
icc_ther <- apps_outcome |>  
  select(ParticipantID, TherapistID, Sx_CASriskacts, Sx_totalacts_sum, AUDIT_sum, SIP_sum, HAMD_sum, BAI_sum, SIDAS_sum, IHS_mean, RS_mean, SOC_concealment_mean, LGBIS_identaffirm) |>  
  pivot_longer(c(Sx_CASriskacts, Sx_totalacts_sum, AUDIT_sum, SIP_sum, HAMD_sum, BAI_sum, SIDAS_sum, IHS_mean, RS_mean, SOC_concealment_mean, LGBIS_identaffirm),
               names_to = "variable",
               values_to = "outcome") |> 
  mutate(variable = as.factor(variable)) |>  
  group_by(variable) |> 
  group_split() |>  
  map_dfr(.f = function(df) {
    lmer(outcome ~ 1 + (1 | ParticipantID) + (1 | TherapistID), REML = T, data = df) |>  
      broom.mixed::tidy(effects = "ran_pars", ddf.method = "Kenward-Roger", scales = "vcov") |>  
      add_column(variable = unique(df$variable), .before = 1) |>  
      group_by(variable) |>  
      summarize(pt_var = first(estimate),
                ther_var = nth(estimate, 2),
                total_var = sum(estimate),
                icc_ther = ther_var / total_var,
                icc_pt = (pt_var + ther_var) / total_var)
  })
```

    boundary (singular) fit: see help('isSingular')
    boundary (singular) fit: see help('isSingular')
    boundary (singular) fit: see help('isSingular')
    boundary (singular) fit: see help('isSingular')
    boundary (singular) fit: see help('isSingular')
    boundary (singular) fit: see help('isSingular')
    boundary (singular) fit: see help('isSingular')

``` r
icc_ther |> select(variable, icc_pt, icc_ther) |> mutate_if(is.numeric, round, digits = 3)
```

    # A tibble: 11 × 3
       variable             icc_pt icc_ther
       <fct>                 <dbl>    <dbl>
     1 AUDIT_sum             0.721    0    
     2 BAI_sum               0.507    0.014
     3 HAMD_sum              0.236    0    
     4 IHS_mean              0.677    0    
     5 LGBIS_identaffirm     0.776    0    
     6 RS_mean               0.624    0    
     7 SIDAS_sum             0.328    0    
     8 SIP_sum               0.63     0.02 
     9 SOC_concealment_mean  0.872    0    
    10 Sx_CASriskacts        0.323    0.141
    11 Sx_totalacts_sum      0.232    0.021

Of note, there were 24 unique therapists and many only had 1 participant
in the sample.

``` r
apps_outcome |> 
  group_by(TherapistID) |> 
  summarize(n_distinct(ParticipantID)) |> 
  rename(number_pts = "n_distinct(ParticipantID)")
```

    # A tibble: 24 × 2
       TherapistID number_pts
       <dbl+lbl>        <int>
     1  1 [AG]              1
     2  2 [AS]              1
     3  3 [BC]              1
     4  4 [CT]              2
     5  5 [EMB]             1
     6  6 [FK]              2
     7  7 [FN]              4
     8  8 [FP]              2
     9  9 [ID]              1
    10 10 [JC]              1
    # ℹ 14 more rows

# Demographics

The demographics that were reported in the primary outcome paper
include: age, race, ethnicity, sex assigned at birth, gender identity,
sexual orientation, education (coded as less than college vs college
degree), employment status, past-year personal income (coded as up to
\$39k, 30-75k or more), and relationship status.

Let’s pull those variables:

``` r
demos <- full |> 
  filter(Assessment == 0) |> 
  select(DQ_age, # age
       DQ_raceEthnicity, # race
       DQ_raceEthnicity.text, # race free response 
       DQ_latinx, # ethnicity 
       DQ_birth.sex, # assigned sex at birth
       DQ_GenID.man:DQ_GenID.other, # gender identity
       DQ_SexualOrientation, # sexual orientation
       DQ_education, # education
       DQ_employment, # employment status
       DQ_income, # personal income past year
       DQ_relationshipStatus# relationship status 
       )
```

Then recode relevant variables to be consistent with main protocol paper
reporting:

``` r
# sexual orientation
demos <- demos |> mutate(DQ_SOr = case_when(DQ_SexualOrientation == 1 ~ 1,
                                   DQ_SexualOrientation == 2 ~ 2,
                                   DQ_SexualOrientation == 3 ~ 2,
                                   DQ_SexualOrientation == 4 ~ 2,
                                   DQ_SexualOrientation == 6 ~ 3,
                                   DQ_SexualOrientation == 7 ~ 4),
                DQ_SOr = recode_factor(DQ_SOr, `1` = "Gay", `2` = "Bi", `3` = "Queer", `4` = "Uncertain")) |> 
  # education
  mutate(DQ_educr = case_when(DQ_education == 1 ~ 1,
                              DQ_education == 2 ~ 1, 
                             .default = 2),
       DQ_educr = recode_factor(DQ_educr, `1` = "Less than college", `2` = "College degree")) |> 
  # income
  mutate(DQ_incr = case_when(DQ_income == 1 ~ 1,
                             DQ_income == 2 ~ 1,
                             DQ_income == 3 ~ 1,
                             .default = 2),
         DQ_incr = recode_factor(DQ_incr, `1` = "up to $29k", `2` = "$30-75k or more")) |> 
  # relationship status
  mutate(DQ_relstatusr = case_when(DQ_relationshipStatus == 4 ~ 1,
                                   DQ_relationshipStatus == 5 ~ 2,
                                   .default = 3),
         DQ_relstatusr = recode_factor(DQ_relstatusr, `1` = "casually dating", `2` = "single", `3` = "in a relationship")) 
```

Pull demographic information for reporting:

``` r
demos |> 
  select(DQ_age, DQ_raceEthnicity, DQ_latinx, DQ_birth.sex, DQ_GenID.man:DQ_GenID.other, DQ_SOr, DQ_educr, DQ_employment, DQ_incr, DQ_relstatusr) |> 
  frq()
```

    What is your age? (DQ_age) <numeric> 
    # total N=57 valid N=57 mean=26.67 sd=4.09

    Value | Label | N | Raw % | Valid % | Cum. %
    --------------------------------------------
       18 |    18 | 0 |  0.00 |    0.00 |   0.00
       19 |    19 | 0 |  0.00 |    0.00 |   0.00
       20 |    20 | 4 |  7.02 |    7.02 |   7.02
       21 |    21 | 1 |  1.75 |    1.75 |   8.77
       22 |    22 | 5 |  8.77 |    8.77 |  17.54
       23 |    23 | 7 | 12.28 |   12.28 |  29.82
       24 |    24 | 3 |  5.26 |    5.26 |  35.09
       25 |    25 | 2 |  3.51 |    3.51 |  38.60
       26 |    26 | 9 | 15.79 |   15.79 |  54.39
       27 |    27 | 1 |  1.75 |    1.75 |  56.14
       28 |    28 | 4 |  7.02 |    7.02 |  63.16
       29 |    29 | 5 |  8.77 |    8.77 |  71.93
       30 |    30 | 6 | 10.53 |   10.53 |  82.46
       31 |    31 | 1 |  1.75 |    1.75 |  84.21
       32 |    32 | 4 |  7.02 |    7.02 |  91.23
       33 |    33 | 2 |  3.51 |    3.51 |  94.74
       34 |    34 | 2 |  3.51 |    3.51 |  98.25
       35 |    35 | 1 |  1.75 |    1.75 | 100.00
     <NA> |  <NA> | 0 |  0.00 |    <NA> |   <NA>

    What racial or ethnic group do you belong to? - Selected Choice (DQ_raceEthnicity) <numeric> 
    # total N=57 valid N=57 mean=1.46 sd=2.06

    Value |                                     Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------------------------------------------
        0 |                                     White | 29 | 50.88 |   50.88 |  50.88
        1 |                    Black/African American | 13 | 22.81 |   22.81 |  73.68
        2 |                                     Asian |  2 |  3.51 |    3.51 |  77.19
        3 |          American Indian or Alaska Native |  0 |  0.00 |    0.00 |  77.19
        4 | Native Hawaiian or Other Pacific Islander |  1 |  1.75 |    1.75 |  78.95
        5 |                               Multiracial | 10 | 17.54 |   17.54 |  96.49
        6 |                    Other (please specify) |  2 |  3.51 |    3.51 | 100.00
     <NA> |                                      <NA> |  0 |  0.00 |    <NA> |   <NA>

    Do you consider yourself Hispanic or Latino? (DQ_latinx) <numeric> 
    # total N=57 valid N=57 mean=0.33 sd=0.48

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 38 | 66.67 |   66.67 |  66.67
        1 |   Yes | 19 | 33.33 |   33.33 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    What was your sex assigned at birth? (DQ_birth.sex) <numeric> 
    # total N=57 valid N=57 mean=1.02 sd=0.13

    Value |    Label |  N | Raw % | Valid % | Cum. %
    ------------------------------------------------
        1 |     Male | 56 | 98.25 |   98.25 |  98.25
        2 |   Female |  1 |  1.75 |    1.75 | 100.00
        3 | Intersex |  0 |  0.00 |    0.00 | 100.00
     <NA> |     <NA> |  0 |  0.00 |    <NA> |   <NA>

    Gender Identity-Man (DQ_GenID.man) <numeric> 
    # total N=57 valid N=57 mean=1.00 sd=0.00

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        1 |   Man | 57 |   100 |     100 |    100
     <NA> |  <NA> |  0 |     0 |    <NA> |   <NA>

    Gender Identity-Woman (DQ_GenID.woman) <numeric> 
    # total N=57 valid N=57 mean=0.00 sd=0.00

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |     0 | 57 |   100 |     100 |    100
        1 | Woman |  0 |     0 |       0 |    100
     <NA> |  <NA> |  0 |     0 |    <NA> |   <NA>

    Gender Identity-Transgender Man (DQ_GenID.transman) <numeric> 
    # total N=57 valid N=57 mean=0.00 sd=0.00

    Value |                 Label |  N | Raw % | Valid % | Cum. %
    -------------------------------------------------------------
        0 |                     0 | 57 |   100 |     100 |    100
        1 | Transgender Man (FtM) |  0 |     0 |       0 |    100
     <NA> |                  <NA> |  0 |     0 |    <NA> |   <NA>

    Gender Identity-Transgender Woman (DQ_GenID.transwoman) <numeric> 
    # total N=57 valid N=57 mean=0.00 sd=0.00

    Value |                   Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------------------------
        0 |                       0 | 57 |   100 |     100 |    100
        1 | Transgender Woman (MtF) |  0 |     0 |       0 |    100
     <NA> |                    <NA> |  0 |     0 |    <NA> |   <NA>

    Gender Identity-Gender Queer (DQ_GenID.genderqueer) <numeric> 
    # total N=57 valid N=57 mean=0.00 sd=0.00

    Value |        Label |  N | Raw % | Valid % | Cum. %
    ----------------------------------------------------
        0 |            0 | 57 |   100 |     100 |    100
        1 | Gender Queer |  0 |     0 |       0 |    100
     <NA> |         <NA> |  0 |     0 |    <NA> |   <NA>

    Gender Identity-Gender Non-Conforming (DQ_GenID.GNC) <numeric> 
    # total N=57 valid N=57 mean=0.00 sd=0.00

    Value |                 Label |  N | Raw % | Valid % | Cum. %
    -------------------------------------------------------------
        0 |                     0 | 57 |   100 |     100 |    100
        1 | Gender Non-Conforming |  0 |     0 |       0 |    100
     <NA> |                  <NA> |  0 |     0 |    <NA> |   <NA>

    Gender Identity-Two Spirit (DQ_GenID.twospirit) <numeric> 
    # total N=57 valid N=57 mean=0.00 sd=0.00

    Value |      Label |  N | Raw % | Valid % | Cum. %
    --------------------------------------------------
        0 |          0 | 57 |   100 |     100 |    100
        1 | Two Spirit |  0 |     0 |       0 |    100
     <NA> |       <NA> |  0 |     0 |    <NA> |   <NA>

    Gender Identity-Hijra (DQ_GenID.hijra) <numeric> 
    # total N=57 valid N=57 mean=0.00 sd=0.00

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |     0 | 57 |   100 |     100 |    100
        1 | Hijra |  0 |     0 |       0 |    100
     <NA> |  <NA> |  0 |     0 |    <NA> |   <NA>

    Gender Identity-Other (DQ_GenID.other) <numeric> 
    # total N=57 valid N=57 mean=0.00 sd=0.00

    Value |                  Label |  N | Raw % | Valid % | Cum. %
    --------------------------------------------------------------
        0 |                      0 | 57 |   100 |     100 |    100
        1 | Other (Please specify) |  0 |     0 |       0 |    100
     <NA> |                   <NA> |  0 |     0 |    <NA> |   <NA>

    DQ_SOr <categorical> 
    # total N=57 valid N=57 mean=1.37 sd=0.67

    Value     |  N | Raw % | Valid % | Cum. %
    -----------------------------------------
    Gay       | 41 | 71.93 |   71.93 |  71.93
    Bi        | 12 | 21.05 |   21.05 |  92.98
    Queer     |  3 |  5.26 |    5.26 |  98.25
    Uncertain |  1 |  1.75 |    1.75 | 100.00
    <NA>      |  0 |  0.00 |    <NA> |   <NA>

    DQ_educr <categorical> 
    # total N=57 valid N=57 mean=1.93 sd=0.26

    Value             |  N | Raw % | Valid % | Cum. %
    -------------------------------------------------
    Less than college |  4 |  7.02 |    7.02 |   7.02
    College degree    | 53 | 92.98 |   92.98 | 100.00
    <NA>              |  0 |  0.00 |    <NA> |   <NA>

    Which of the following best describes your current employment status? (DQ_employment) <numeric> 
    # total N=57 valid N=57 mean=2.42 sd=1.91

    Value |                                                         Label |  N
    --------------------------------------------------------------------------
        1 |                                 Full-time (40 hours per week) | 26
        2 |                       Part-time (less than 40 hours per week) | 13
        3 |                            Part-time work - full time student |  8
        4 |           Permanently or temporarily disabled and NOT working |  0
        5 | Permanently or temporarily disabled BUT working off the books |  1
        6 |                                          Unemployed - Student |  6
        7 |                                            Unemployed - Other |  3
     <NA> |                                                          <NA> |  0

    Value | Raw % | Valid % | Cum. %
    --------------------------------
        1 | 45.61 |   45.61 |  45.61
        2 | 22.81 |   22.81 |  68.42
        3 | 14.04 |   14.04 |  82.46
        4 |  0.00 |    0.00 |  82.46
        5 |  1.75 |    1.75 |  84.21
        6 | 10.53 |   10.53 |  94.74
        7 |  5.26 |    5.26 | 100.00
     <NA> |  0.00 |    <NA> |   <NA>

    DQ_incr <categorical> 
    # total N=57 valid N=57 mean=1.46 sd=0.50

    Value           |  N | Raw % | Valid % | Cum. %
    -----------------------------------------------
    up to $29k      | 31 | 54.39 |   54.39 |  54.39
    $30-75k or more | 26 | 45.61 |   45.61 | 100.00
    <NA>            |  0 |  0.00 |    <NA> |   <NA>

    DQ_relstatusr <categorical> 
    # total N=57 valid N=57 mean=2.07 sd=0.42

    Value             |  N | Raw % | Valid % | Cum. %
    -------------------------------------------------
    casually dating   |  3 |  5.26 |    5.26 |   5.26
    single            | 47 | 82.46 |   82.46 |  87.72
    in a relationship |  7 | 12.28 |   12.28 | 100.00
    <NA>              |  0 |  0.00 |    <NA> |   <NA>

# MINI diagnoses

``` r
full |>
  select(ParticipantID, Assessment, contains("MINI")) |> 
  filter(Assessment == 0) |> 
  select(MINI_Depression, MINI_Dysthymia, MINI_PanicDisorder, MINI_Agoraphobia, MINI_SocialAnxiety, MINI_OCD, MINI_PTSD, MINI_GAD, MINI_AUD, MINI_SUD) |> 
    frq()
```

    Did the pt meet criteria for Depression on the MINI? (MINI_Depression) <numeric> 
    # total N=57 valid N=57 mean=0.81 sd=0.40

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 11 | 19.30 |   19.30 |  19.30
        1 |   Yes | 46 | 80.70 |   80.70 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    Did the pt meet criteria for Dysthymia on the MINI? (MINI_Dysthymia) <numeric> 
    # total N=57 valid N=57 mean=0.39 sd=0.49

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 35 | 61.40 |   61.40 |  61.40
        1 |   Yes | 22 | 38.60 |   38.60 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    Did the pt meet criteria for PanicDisorder on the MINI? (MINI_PanicDisorder) <numeric> 
    # total N=57 valid N=57 mean=0.30 sd=0.46

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 40 | 70.18 |   70.18 |  70.18
        1 |   Yes | 17 | 29.82 |   29.82 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    Did the pt meet criteria for Agoraphobia on the MINI? (MINI_Agoraphobia) <numeric> 
    # total N=57 valid N=57 mean=0.28 sd=0.45

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 41 | 71.93 |   71.93 |  71.93
        1 |   Yes | 16 | 28.07 |   28.07 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    Did the pt meet criteria for SocialAnxiety on the MINI? (MINI_SocialAnxiety) <numeric> 
    # total N=57 valid N=57 mean=0.49 sd=0.50

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 29 | 50.88 |   50.88 |  50.88
        1 |   Yes | 28 | 49.12 |   49.12 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    Did the pt meet criteria for OCD on the MINI? (MINI_OCD) <numeric> 
    # total N=57 valid N=57 mean=0.23 sd=0.42

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 44 | 77.19 |   77.19 |  77.19
        1 |   Yes | 13 | 22.81 |   22.81 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    Did the pt meet criteria for PTSD on the MINI? (MINI_PTSD) <numeric> 
    # total N=57 valid N=57 mean=0.18 sd=0.38

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 47 | 82.46 |   82.46 |  82.46
        1 |   Yes | 10 | 17.54 |   17.54 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    Did the pt meet criteria for GAD on the MINI? (MINI_GAD) <numeric> 
    # total N=57 valid N=57 mean=0.44 sd=0.50

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 32 | 56.14 |   56.14 |  56.14
        1 |   Yes | 25 | 43.86 |   43.86 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    Did the pt meet criteria for Alcohol Use Disorder on the MINI? (MINI_AUD) <numeric> 
    # total N=57 valid N=57 mean=0.23 sd=0.42

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 44 | 77.19 |   77.19 |  77.19
        1 |   Yes | 13 | 22.81 |   22.81 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

    Did the pt meet criteria for any Substance Use Disorder on the MINI? (MINI_SUD) <numeric> 
    # total N=57 valid N=57 mean=0.56 sd=0.50

    Value | Label |  N | Raw % | Valid % | Cum. %
    ---------------------------------------------
        0 |    No | 25 | 43.86 |   43.86 |  43.86
        1 |   Yes | 32 | 56.14 |   56.14 | 100.00
     <NA> |  <NA> |  0 |  0.00 |    <NA> |   <NA>

# Save data file

``` r
saveRDS(apps_outcome, file = "data/apps_outcome.rds")
```

# Output

## Correlation matrix

``` r
# computes correlation matrix in a dataframe that can be saved and easily reformated for publication
corr_results <- correlation(data = apps_outcome, select = c("Sx_CASriskacts", "AUDIT_sum", "SIP_sum", "HAMD_sum", "BAI_sum", "SIDAS_sum", "IHS_mean", "RS_mean", "SOC_concealment_mean", "LGBIS_identaffirm", "Average_Affirm")) |> 
  as.data.frame() |> 
  select(Parameter1, Parameter2, r, p) |> 
  mutate(r = sprintf("%.2f", r), # keeps the trailing zeros for correlations
         ps = ifelse(p < .10, '+', ''),
         ps = ifelse(p < .05, '*', ps),
         ps = ifelse(p < .01, '**', ps),
         ps = ifelse(p < .001, '***', ps),
         r = as.character(r)) |> 
  tidyr::unite(r_p, c("r", "ps"), sep = "") |> 
  select(-p)

write.csv(corr_results, "output/corr_table.csv")

corr_results
```

                 Parameter1           Parameter2      r_p
    1        Sx_CASriskacts            AUDIT_sum    0.22*
    2        Sx_CASriskacts              SIP_sum     0.20
    3        Sx_CASriskacts             HAMD_sum     0.04
    4        Sx_CASriskacts              BAI_sum    -0.02
    5        Sx_CASriskacts            SIDAS_sum    -0.06
    6        Sx_CASriskacts             IHS_mean     0.06
    7        Sx_CASriskacts              RS_mean     0.04
    8        Sx_CASriskacts SOC_concealment_mean    -0.10
    9        Sx_CASriskacts    LGBIS_identaffirm    -0.06
    10       Sx_CASriskacts       Average_Affirm     0.13
    11            AUDIT_sum              SIP_sum  0.64***
    12            AUDIT_sum             HAMD_sum    0.22*
    13            AUDIT_sum              BAI_sum    0.24*
    14            AUDIT_sum            SIDAS_sum     0.16
    15            AUDIT_sum             IHS_mean   0.27**
    16            AUDIT_sum              RS_mean     0.18
    17            AUDIT_sum SOC_concealment_mean    -0.07
    18            AUDIT_sum    LGBIS_identaffirm    -0.10
    19            AUDIT_sum       Average_Affirm    -0.02
    20              SIP_sum             HAMD_sum   0.26**
    21              SIP_sum              BAI_sum  0.29***
    22              SIP_sum            SIDAS_sum     0.16
    23              SIP_sum             IHS_mean    0.21+
    24              SIP_sum              RS_mean   0.27**
    25              SIP_sum SOC_concealment_mean    -0.05
    26              SIP_sum    LGBIS_identaffirm    -0.04
    27              SIP_sum       Average_Affirm    -0.12
    28             HAMD_sum              BAI_sum  0.37***
    29             HAMD_sum            SIDAS_sum  0.29***
    30             HAMD_sum             IHS_mean    0.24*
    31             HAMD_sum              RS_mean   0.26**
    32             HAMD_sum SOC_concealment_mean     0.02
    33             HAMD_sum    LGBIS_identaffirm    -0.11
    34             HAMD_sum       Average_Affirm    -0.12
    35              BAI_sum            SIDAS_sum     0.15
    36              BAI_sum             IHS_mean     0.10
    37              BAI_sum              RS_mean    0.23*
    38              BAI_sum SOC_concealment_mean    -0.06
    39              BAI_sum    LGBIS_identaffirm    -0.06
    40              BAI_sum       Average_Affirm    -0.02
    41            SIDAS_sum             IHS_mean     0.13
    42            SIDAS_sum              RS_mean     0.14
    43            SIDAS_sum SOC_concealment_mean     0.09
    44            SIDAS_sum    LGBIS_identaffirm     0.02
    45            SIDAS_sum       Average_Affirm    -0.15
    46             IHS_mean              RS_mean  0.34***
    47             IHS_mean SOC_concealment_mean  0.35***
    48             IHS_mean    LGBIS_identaffirm -0.56***
    49             IHS_mean       Average_Affirm    -0.02
    50              RS_mean SOC_concealment_mean     0.18
    51              RS_mean    LGBIS_identaffirm    -0.02
    52              RS_mean       Average_Affirm    -0.01
    53 SOC_concealment_mean    LGBIS_identaffirm -0.31***
    54 SOC_concealment_mean       Average_Affirm  -0.27**
    55    LGBIS_identaffirm       Average_Affirm     0.05

## Descriptives

Pull the descriptives for main study variables for reporting:

``` r
descriptives_results <- apps_outcome |>  
  select(Sx_CASriskacts, AUDIT_sum, SIP_sum, HAMD_sum, BAI_sum, SIDAS_sum, IHS_mean, RS_mean, SOC_concealment_mean, LGBIS_identaffirm, Average_Affirm) |>  
  summarize(across(everything(), list(mean = ~ mean(.x, na.rm = T),
                                      median = ~ median(.x, na.rm = T),
                                      sd = ~ sd(.x, na.rm = T),
                                      min = ~ min(.x, na.rm = T),
                                      max = ~ max(.x, na.rm = T)),
                   .names = "{.col}-{.fn}")
            ) |>  
  pivot_longer(cols = everything(),
               names_to = c("ColNames", ".value"), 
               names_sep = "-",
               ) |> 
  mutate_if(is.numeric, round, digits = 2) |>
  mutate(sd = paste0('(', sd) %>% paste0(')')) |> 
  tidyr::unite(mean_sd, c("mean", "sd"), sep = " ") |> 
  tidyr::unite(range, c("min", "max"), sep = "-") |> 
  select(-median)

write.csv(descriptives_results, "output/descriptives.csv")

descriptives_results
```

    # A tibble: 11 × 3
       ColNames             mean_sd      range 
       <chr>                <chr>        <chr> 
     1 Sx_CASriskacts       4.47 (6.41)  0-45  
     2 AUDIT_sum            8.16 (6.28)  0-29  
     3 SIP_sum              3.19 (3.68)  0-15  
     4 HAMD_sum             12.27 (6.66) 0-32  
     5 BAI_sum              18.8 (11.38) 0-55  
     6 SIDAS_sum            2.48 (5.7)   0-39  
     7 IHS_mean             1.66 (0.63)  1-3.56
     8 RS_mean              12.2 (7.28)  1-36  
     9 SOC_concealment_mean 1.79 (0.69)  1-3.4 
    10 LGBIS_identaffirm    4.56 (1.28)  1-6   
    11 Average_Affirm       1.38 (1.2)   0-4.75

# Session info

``` r
sessionInfo()
```

    R version 4.4.3 (2025-02-28 ucrt)
    Platform: x86_64-w64-mingw32/x64
    Running under: Windows 11 x64 (build 22631)

    Matrix products: default


    locale:
    [1] LC_COLLATE=English_United States.utf8 
    [2] LC_CTYPE=English_United States.utf8   
    [3] LC_MONETARY=English_United States.utf8
    [4] LC_NUMERIC=C                          
    [5] LC_TIME=English_United States.utf8    

    time zone: America/New_York
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] lmerTest_3.1-3     lme4_1.1-37        Matrix_1.7-2       gridExtra_2.3     
     [5] psych_2.5.3        see_0.11.0         report_0.6.1       parameters_0.24.2 
     [9] performance_0.13.0 modelbased_0.10.0  insight_1.1.0      effectsize_1.0.0  
    [13] datawizard_1.0.2   correlation_0.8.7  bayestestR_0.15.2  easystats_0.7.4   
    [17] haven_2.5.4        sjmisc_2.8.10      lubridate_1.9.4    forcats_1.0.0     
    [21] stringr_1.5.1      dplyr_1.1.4        purrr_1.0.4        readr_2.1.5       
    [25] tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0   

    loaded via a namespace (and not attached):
     [1] sjlabelled_1.2.0    tidyselect_1.2.1    farver_2.1.2       
     [4] fastmap_1.2.0       digest_0.6.37       estimability_1.5.1 
     [7] timechange_0.3.0    lifecycle_1.0.4     magrittr_2.0.3     
    [10] compiler_4.4.3      rlang_1.1.5         tools_4.4.3        
    [13] utf8_1.2.4          yaml_2.3.10         knitr_1.50         
    [16] labeling_0.4.3      mnormt_2.1.1        withr_3.0.2        
    [19] numDeriv_2016.8-1.1 grid_4.4.3          xtable_1.8-4       
    [22] colorspace_2.1-1    future_1.34.0       emmeans_1.11.0     
    [25] globals_0.16.3      scales_1.3.0        MASS_7.3-64        
    [28] cli_3.6.4           mvtnorm_1.3-3       crayon_1.5.3       
    [31] rmarkdown_2.29      reformulas_0.4.0    generics_0.1.3     
    [34] rstudioapi_0.17.1   tzdb_0.5.0          minqa_1.2.8        
    [37] splines_4.4.3       parallel_4.4.3      vctrs_0.6.5        
    [40] boot_1.3-31         jsonlite_2.0.0      hms_1.1.3          
    [43] listenv_0.9.1       parallelly_1.43.0   glue_1.8.0         
    [46] nloptr_2.2.1        codetools_0.2-20    stringi_1.8.7      
    [49] gtable_0.3.6        broom.mixed_0.2.9.6 munsell_0.5.1      
    [52] pillar_1.10.2       furrr_0.3.1         htmltools_0.5.8.1  
    [55] R6_2.6.1            Rdpack_2.6.3        evaluate_1.0.3     
    [58] lattice_0.22-6      rbibutils_2.3       backports_1.5.0    
    [61] broom_1.0.8         snakecase_0.11.1    Rcpp_1.0.14        
    [64] coda_0.19-4.1       nlme_3.1-167        mgcv_1.9-1         
    [67] xfun_0.52           pkgconfig_2.0.3    
