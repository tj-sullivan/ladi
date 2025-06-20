---
title: "Reliability Analyses"
format: html
keep-md: true 
embed-resources: true
---




# Load packages


::: {.cell}

```{.r .cell-code}
library(tidyverse)
library(sjmisc)
library(haven)
library(irr)
```
:::



# Load data [not shown for peer review]






# Calculate ICCs

First, make an overall variable for the average rating (early + mid) for each coder:


::: {.cell}

```{.r .cell-code}
apps <- apps |> 
  row_means(Early_Affirm_Coder1, Mid_Affirm_Coder1, n = Inf, var = "Average_Affirm_Coder1") |> 
  row_means(Early_Affirm_Coder2, Mid_Affirm_Coder2, n = Inf, var = "Average_Affirm_Coder2")
```
:::



The `irr` package was used to calculate ICCs. The first part of the code specifies the type of model, which here is set two "twoway" b/c the design is fully crossed, such that we have 20% of videos rated by all coders to compare for reliability -- it would be oneway if not all coders rated all participants used for reliability. The second part specifies the type of analyses, which here is set to "agreement" b/c we want people to be closer in absolute value, whereas consistency will calculate this based on rank order. The third part specifies the unit, which here is set to "average" b/c we took the average ratings across coders. The last two options (r0) and conf.level are for null hypothesis significance testing. 

ICC for early sessions:


::: {.cell}

```{.r .cell-code}
select(apps, Early_Affirm_Coder1, Early_Affirm_Coder2) |> 
  data.matrix() |> 
  icc(c("twoway"), type = c("agreement"), unit = c("average"), r0 = 0, conf.level = 0.95)
```

::: {.cell-output .cell-output-stdout}

```
 Average Score Intraclass Correlation

   Model: twoway 
   Type : agreement 

   Subjects = 57 
     Raters = 2 
   ICC(A,2) = 0.886

 F-Test, H0: r0 = 0 ; H1: r0 > 0 
 F(56,50.1) = 9.19 , p = 1.96e-13 

 95%-Confidence Interval for ICC Population Values:
  0.803 < ICC < 0.933
```


:::
:::



ICC for mid sessions:



::: {.cell}

```{.r .cell-code}
select(apps, Mid_Affirm_Coder1, Mid_Affirm_Coder2) |> 
  data.matrix() |> 
  icc(c("twoway"), type = c("agreement"), unit = c("average"), r0 = 0, conf.level = 0.95)
```

::: {.cell-output .cell-output-stdout}

```
 Average Score Intraclass Correlation

   Model: twoway 
   Type : agreement 

   Subjects = 57 
     Raters = 2 
   ICC(A,2) = 0.764

 F-Test, H0: r0 = 0 ; H1: r0 > 0 
 F(56,56.1) = 4.18 , p = 1.48e-07 

 95%-Confidence Interval for ICC Population Values:
  0.598 < ICC < 0.861
```


:::
:::



ICC for average sessions:


::: {.cell}

```{.r .cell-code}
select(apps, Average_Affirm_Coder1, Average_Affirm_Coder2) |> 
  data.matrix() |> 
  icc(c("twoway"), type = c("agreement"), unit = c("average"), r0 = 0, conf.level = 0.95)
```

::: {.cell-output .cell-output-stdout}

```
 Average Score Intraclass Correlation

   Model: twoway 
   Type : agreement 

   Subjects = 57 
     Raters = 2 
   ICC(A,2) = 0.877

 F-Test, H0: r0 = 0 ; H1: r0 > 0 
 F(56,56.3) = 8.28 , p = 1.46e-13 

 95%-Confidence Interval for ICC Population Values:
  0.792 < ICC < 0.928
```


:::
:::




# Wrangle data for item-by-item code 


::: {.cell}

```{.r .cell-code}
data <- apps |> 
  select(ParticipantID, 
         Early_Session_Coder1_Affirm_Item1:Early_Session_Coder1_Affirm_Item10,
         Early_Session_Coder2_Affirm_Item1:Early_Session_Coder2_Affirm_Item10,
         Mid_Session_Coder1_Affirm_Item1:Mid_Session_Coder1_Affirm_Item10,
         Mid_Session_Coder2_Affirm_Item1:Mid_Session_Coder2_Affirm_Item10) |> 
  rename(Early1_r1 = Early_Session_Coder1_Affirm_Item1,
         Early2_r1 = Early_Session_Coder1_Affirm_Item2,
         Early3_r1 = Early_Session_Coder1_Affirm_Item3,
         Early4_r1 = Early_Session_Coder1_Affirm_Item4,
         Early5_r1 = Early_Session_Coder1_Affirm_Item5,
         Early6_r1 = Early_Session_Coder1_Affirm_Item6,
         Early7_r1 = Early_Session_Coder1_Affirm_Item7,
         Early8_r1 = Early_Session_Coder1_Affirm_Item8,
         Early9_r1 = Early_Session_Coder1_Affirm_Item9,
         Early10_r1 = Early_Session_Coder1_Affirm_Item10) |> 
  rename(Early1_r2 = Early_Session_Coder2_Affirm_Item1,
         Early2_r2 = Early_Session_Coder2_Affirm_Item2,
         Early3_r2 = Early_Session_Coder2_Affirm_Item3,
         Early4_r2 = Early_Session_Coder2_Affirm_Item4,
         Early5_r2 = Early_Session_Coder2_Affirm_Item5,
         Early6_r2 = Early_Session_Coder2_Affirm_Item6,
         Early7_r2 = Early_Session_Coder2_Affirm_Item7,
         Early8_r2 = Early_Session_Coder2_Affirm_Item8,
         Early9_r2 = Early_Session_Coder2_Affirm_Item9,
         Early10_r2 = Early_Session_Coder2_Affirm_Item10) |> 
  rename(Mid1_r1 = Mid_Session_Coder1_Affirm_Item1,
         Mid2_r1 = Mid_Session_Coder1_Affirm_Item2,
         Mid3_r1 = Mid_Session_Coder1_Affirm_Item3,
         Mid4_r1 = Mid_Session_Coder1_Affirm_Item4,
         Mid5_r1 = Mid_Session_Coder1_Affirm_Item5,
         Mid6_r1 = Mid_Session_Coder1_Affirm_Item6,
         Mid7_r1 = Mid_Session_Coder1_Affirm_Item7,
         Mid8_r1 = Mid_Session_Coder1_Affirm_Item8,
         Mid9_r1 = Mid_Session_Coder1_Affirm_Item9,
         Mid10_r1 = Mid_Session_Coder1_Affirm_Item10) |> 
  rename(Mid1_r2 = Mid_Session_Coder2_Affirm_Item1,
         Mid2_r2 = Mid_Session_Coder2_Affirm_Item2,
         Mid3_r2 = Mid_Session_Coder2_Affirm_Item3,
         Mid4_r2 = Mid_Session_Coder2_Affirm_Item4,
         Mid5_r2 = Mid_Session_Coder2_Affirm_Item5,
         Mid6_r2 = Mid_Session_Coder2_Affirm_Item6,
         Mid7_r2 = Mid_Session_Coder2_Affirm_Item7,
         Mid8_r2 = Mid_Session_Coder2_Affirm_Item8,
         Mid9_r2 = Mid_Session_Coder2_Affirm_Item9,
         Mid10_r2 = Mid_Session_Coder2_Affirm_Item10) |> 
  zap_labels() |> 
  zap_label() |> 
  zap_formats() |> 
  zap_widths()
```
:::



# Holley & Guilford's G calculation



::: {.cell}

```{.r .cell-code}
# Tell R which columns correspond to rater 1 and rater 2 respectively
#grep is cool for pattern identification
rater1_cols <- grep("_r1$", names(data), value = TRUE)
#but make sure that column names are the same
#because R is finding and replacing here
rater2_cols <- gsub("_r1$", "_r2", rater1_cols)

# Compute Guilford's G for each item with a loop
guilford_g <- numeric(length(rater1_cols))
item_names <- gsub("_r1$", "", rater1_cols)
# 
for (i in seq_along(rater1_cols)) {
  r1 <- data[[rater1_cols[i]]]
  r2 <- data[[rater2_cols[i]]]

  agreements <- sum(r1 == r2, na.rm = TRUE)
  disagreements <- sum(r1 != r2, na.rm = TRUE)

  if ((agreements + disagreements) > 0) {
    guilford_g[i] <- ((agreements / (agreements + disagreements))-0.5)/(1-0.5)
    #using 0.5 here because 1-Pc, chance is (1/2) 0.5 with dichotomous vars
  } else {
    guilford_g[i] <- NA #to avoid dividing by 0
  }
}
```
:::



# Percent agreement on presence and absence

Code created from calculations presented by Xu & Lorber (2014)


::: {.cell}

```{.r .cell-code}
ppos <- function(df) {
  cts <- table(df)
  
  if (sum(dim(cts)) != 4){
  df[,1] <- factor(df[,1], levels = c("0", "1"))
  df[,2] <- factor(df[,2], levels = c("0", "1"))
  
  cts <- table(df)
  }
  
  a <- cts[2,2] # represents count of agreement of presence
  b <- cts[1,2] # represents column-wise disagreement
  c <- cts[2,1] # represents row-wise disagreement
  d <- cts[1,1] # represents count of agreement on absence

  ppos <- (2*a)/((a+b)+(a+c))
  
  print(ppos)
}


pneg <- function(df) {
  cts <- table(df)
  
  if (sum(dim(cts)) != 4){
  df[,1] <- factor(df[,1], levels = c("0", "1"))
  df[,2] <- factor(df[,2], levels = c("0", "1"))
  
  cts <- table(df)
  }
  
  a <- cts[2,2] # represents count of agreement of presence
  b <- cts[1,2] # represents column-wise disagreement
  c <- cts[2,1] # represents row-wise disagreement
  d <- cts[1,1] # represents count of agreement on absence

  pneg <- (2*d)/((b+d)+(c+d))

  print(pneg)
}


ppos_df <- numeric(length(rater1_cols))
for (i in seq_along(rater1_cols)) {
  r1 <- data[[rater1_cols[i]]]
  r2 <- data[[rater2_cols[i]]]
  
  tb <- as.data.frame(cbind(r1, r2))

  ppos_df[i] <- ppos(tb)
}
```

::: {.cell-output .cell-output-stdout}

```
[1] 0.5945946
[1] 0.7692308
[1] 0.7692308
[1] 0.7777778
[1] 0.2
[1] 0.4
[1] 1
[1] 0.5
[1] 0.5
[1] 0.7692308
[1] 0.5641026
[1] 0.7118644
[1] 0.5
[1] 0.6666667
[1] 0.3333333
[1] 0
[1] 0
[1] 0
[1] 0.5
[1] 0
```


:::

```{.r .cell-code}
pneg_df <- numeric(length(rater1_cols))
for (i in seq_along(rater1_cols)) {
  r1 <- data[[rater1_cols[i]]]
  r2 <- data[[rater2_cols[i]]]
  
  tb <- as.data.frame(cbind(r1, r2))

  pneg_df[i] <- pneg(tb)
}
```

::: {.cell-output .cell-output-stdout}

```
[1] 0.8051948
[1] 0.8064516
[1] 0.970297
[1] 0.9583333
[1] 0.9230769
[1] 0.9724771
[1] 1
[1] 0.9818182
[1] 0.9818182
[1] 0.970297
[1] 0.7733333
[1] 0.6909091
[1] 0.9622642
[1] 0.9247312
[1] 0.9215686
[1] 0.9636364
[1] 0.9821429
[1] 0.9911504
[1] 0.9818182
[1] 0.9636364
```


:::
:::



# Output



::: {.cell}

```{.r .cell-code}
#Create and print summary table with Guilford's G for each item
g_table <- tibble(
  Item = item_names,
  Guilford_G = round(guilford_g, 3),
  ppos = round(ppos_df, 4),
  pneg = round(pneg_df, 4)
)
g_table
```

::: {.cell-output .cell-output-stdout}

```
# A tibble: 20 × 4
   Item    Guilford_G  ppos  pneg
   <chr>        <dbl> <dbl> <dbl>
 1 Early1       0.474 0.595 0.805
 2 Early2       0.579 0.769 0.806
 3 Early3       0.895 0.769 0.970
 4 Early4       0.86  0.778 0.958
 5 Early5       0.719 0.2   0.923
 6 Early6       0.895 0.4   0.972
 7 Early7       1     1     1    
 8 Early8       0.93  0.5   0.982
 9 Early9       0.93  0.5   0.982
10 Early10      0.895 0.769 0.970
11 Mid1         0.404 0.564 0.773
12 Mid2         0.404 0.712 0.691
13 Mid3         0.86  0.5   0.962
14 Mid4         0.754 0.667 0.925
15 Mid5         0.719 0.333 0.922
16 Mid6         0.86  0     0.964
17 Mid7         0.93  0     0.982
18 Mid8         0.965 0     0.991
19 Mid9         0.93  0.5   0.982
20 Mid10        0.86  0     0.964
```


:::
:::



Average G:


::: {.cell}

```{.r .cell-code}
#Average G: just averages across all 10 Gs above
average_g <- mean(guilford_g, na.rm = TRUE)
cat("\nAverage Guilford's G:", round(average_g, 3), "\n")
```

::: {.cell-output .cell-output-stdout}

```

Average Guilford's G: 0.793 
```


:::
:::



Early G average:


::: {.cell}

```{.r .cell-code}
early_g <- mean(guilford_g[1:10], na.rm = T)
cat("\nEarly Sessions Guilford's G:", round(early_g, 3), "\n")
```

::: {.cell-output .cell-output-stdout}

```

Early Sessions Guilford's G: 0.818 
```


:::
:::



Mid G average:


::: {.cell}

```{.r .cell-code}
mid_g <- mean(guilford_g[11:20], na.rm = T)
cat("\nMid Sessions Guilford's G:", round(mid_g, 3), "\n")
```

::: {.cell-output .cell-output-stdout}

```

Mid Sessions Guilford's G: 0.768 
```


:::
:::



Save to a .csv file:


::: {.cell}

```{.r .cell-code}
write.csv(g_table, "output/reliability.csv")
```
:::



# Session info



::: {.cell}

```{.r .cell-code}
sessionInfo()
```

::: {.cell-output .cell-output-stdout}

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
 [1] irr_0.84.1      lpSolve_5.6.23  haven_2.5.4     sjmisc_2.8.10  
 [5] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
 [9] purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
[13] ggplot2_3.5.1   tidyverse_2.0.0

loaded via a namespace (and not attached):
 [1] gtable_0.3.6      jsonlite_2.0.0    compiler_4.4.3    tidyselect_1.2.1 
 [5] scales_1.3.0      yaml_2.3.10       fastmap_1.2.0     R6_2.6.1         
 [9] generics_0.1.3    sjlabelled_1.2.0  knitr_1.50        insight_1.1.0    
[13] munsell_0.5.1     pillar_1.10.2     tzdb_0.5.0        rlang_1.1.5      
[17] utf8_1.2.4        stringi_1.8.7     xfun_0.52         timechange_0.3.0 
[21] cli_3.6.4         withr_3.0.2       magrittr_2.0.3    digest_0.6.37    
[25] grid_4.4.3        rstudioapi_0.17.1 hms_1.1.3         lifecycle_1.0.4  
[29] vctrs_0.6.5       evaluate_1.0.3    glue_1.8.0        colorspace_2.1-1 
[33] rmarkdown_2.29    tools_4.4.3       pkgconfig_2.0.3   htmltools_0.5.8.1
```


:::
:::
