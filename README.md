Varroa control limits virus spillover: data files and scripts
================

This repository contains supporting files and code associated with the
submitted manuscript entitled “Controlling Varroa destructor
infestations in managed honey bee colonies limits virus spillover into
wild pollinator communities.” Three csv files provide qPCR data for DWV
screens of (1) honey bees, (2) goldenrod flowers, and (3) all wild
pollinators. There are additionally two metadata csv files that
characterize (1) data on apiary size (number of colonies) and varroa
loads from all colonies sampled (even if not screened for DWV) and (2)
data for all floral quadrats on honey bee visitation and distance from
the apiary.

## Description of the data and file structure

### Data files:

#### (1) **honey.bee.metadata.csv** : varroa loads in each colony sampled (4 or 10 colonies per apiary, based on the size of the apiary)

    - Sample_ID : unique identification code for each colony sampled 
    - Sampling_Date : the date each colony was sampled to measure varroa loads
    - Day_of_Year : the 2020 ordinal day number for the date the sample was collected
    - Apiary_ID : identification code for the apiary at which the sample was collected
    - Num_colonies : the number of colonies in the apiary where the colony was located
    - Percent_varroa : the number of mites per 100 bees (percent infestation level) in the colony
    - Virus_Screening : denotes whether the colony was included in subsequent DWV screening (Yes = bees from the colony were screened for DWV, No = bees from the colony were not screened for DWV)

*Note that only a subset of colonies included in this dataset were
selected for DWV screening, but data from all colonies were included to
calculate apiary-average varroa levels.*

#### (2) **honey.bee.dwv.data.csv** : qPCR loads of DWV and reference gene in honey bee samples

    - Sample_ID : unique identification code for each colony sampled 
    - 28S : average Cq value of the 28S reference gene
    - DWV : average Cq value of DWV 
    - Apiary_ID : identification code for the apiary at which the sample was collected

#### (3) **flower.metadata.csv** : information on all flowers collected for the duration of the project

    - Sample_ID : unique identification code for each flower collected
    - Apiary_ID : identification code for the apiary at which the sample was collected
    - Sampling_Date : the date the flower was collected and observed
    - Day_of_Year : the 2020 ordinal day number for the date the sample was collected
    - `Quadrat #` : below the site level, this denotes the quadrat in which the flower was collected (4-10 quadrats per sampling date)
    - `# of HB Visits` :  the number of honey bee visits to goldenrod flowers within the quadrat from which the flower was collected, during a five-minute period (NA = no honey bee observation data was collected from the quadrat)
    - `Distance to colonies (m)` : the distance of the quadrat in which the flower was collected from the nearest honey bee colony in the apiary (in meters)
    - Virus_Screening : denotes whether the flower was included in subsequent DWV screening (Yes = flowers were screened for DWV, No = flowers were not screened for DWV)

*Note that only a subset of flowers included in this dataset were
selected for DWV screening, but data from all quadrats were included to
calculate apiary-average honey bee visitation estimates.*

#### (4) **flower.dwv.data.csv** : presence / absence of DWV on goldenrod flowers

    - Sample_ID : unique identification code for each flower collected
    - Apiary_ID : identification code for the apiary at which the sample was collected
    - DWV_Cq : the Cq value from qPCR screening for the DWV target in the sample
    - DWV_presence : the presence or absence of DWV in each sample (0 = absent, 1 = present; Cq values above 35 were treated as non-detections and coded as "0")

#### (5) **all.wild.pollinators.dwv.data.csv** : presence / absence and load (when applicable) of DWV in wild pollinators

    - Sample_ID : unique identification code for each individual pollinator (in the format concatenating Apiary_ID, sample date (MMDDYY), and pollinator number)
    - Taxa : specifies whether the sample is a bumble bee (target genus), other wild bee, or a syrphid fly 
    - Family : specifies the family to which the sample belongs (specimens categorized as "unknown" if family-level identification could not be determined)
    - Subfamily : specifies the subfamily to which the sample belongs, used for syrphid flies only (specimens categorized as "unknown" if subfamily-level identification could not be determined)
    - Genus : specifies the genus to which the sample belongs (specimens categorized as "unknown" if genus-level identification could not be determined)
    - Species : specifies the species to which the sample belongs (specimens described as "Genus sp." if species-level identification could not be determined and categorized as "unknown" if genus-level identification could not be determined)
    - Sex : sex of the specimen (M = male, F = female, U = undetermined, Q = queen (bumble bee only), W = worker (bumble bee only))
    - Apiary_ID : identification code for the apiary at which the sample was collected
    - Day_of_Year : the 2020 ordinal day number for the date the sample was collected
    - DWV_Cq : the Cq value from qPCR screening for the DWV target in the sample (average value from 3 technical replicates for syrphid flies and wild bees, single value for bumble bees; NA denotes non-detections)
    - Reference_Cq : the Cq value from qPCR screening for the reference target in the sample (average value from 3 technical replicates for syrphid flies and wild bees; bumble bees and two Xylocopa virginica were not screened for a reference gene and this is denoted "none")
    - DWV_presence : the presence or absence of DWV in each sample (0 = absent, 1 = present; Cq values above 35 were treated as non-detections and coded as "0")

## Code/Software

All statistical tests were performed using R Statistical Software
(v4.3.1, R Core Team, 2023) using the packages lme4 (Bates et al.,
2015), emmeans (Lenth, 2023), and glmmTMB (Brooks et al., 2017). All
model fits were assessed using the DHARMa package (Hartig, 2022). The
code to reproduce all the statistical tests and figures included in this
manuscript can be found in **varroa.virus.spillover.final.script.R**.

The full list of packages used for statistical analysis, data carpentry
and visualization as well as their version info are listed below:

``` r
library(dplyr)
library(tidyverse)
library(forcats)
library(glmmTMB)
library(lme4)
library(lmerTest)
library(emmeans)
library(DHARMa)
library(ggplot2)
library(PNWColors)
library(NatParksPalettes)
library(ggbeeswarm)
library(gridExtra)
library(ggpubr)
library(cowplot)
library(sessioninfo)

sessioninfo::session_info()
```

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.1 (2025-06-13 ucrt)
    ##  os       Windows 11 x64 (build 22631)
    ##  system   x86_64, mingw32
    ##  ui       RTerm
    ##  language (EN)
    ##  collate  English_United States.utf8
    ##  ctype    English_United States.utf8
    ##  tz       America/New_York
    ##  date     2026-07-09
    ##  pandoc   3.8.3 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)
    ##  quarto   1.9.38 @ C:\\PROGRA~1\\RStudio\\RESOUR~1\\app\\bin\\quarto\\bin\\quarto.exe
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  ! package          * version    date (UTC) lib source
    ##    abind              1.4-8      2024-09-12 [1] CRAN (R 4.5.0)
    ##    backports          1.5.0      2024-05-23 [1] CRAN (R 4.5.0)
    ##    beeswarm           0.4.0      2021-06-01 [1] CRAN (R 4.5.0)
    ##    boot               1.3-31     2024-08-28 [1] CRAN (R 4.5.1)
    ##    broom              1.0.8      2025-03-28 [1] CRAN (R 4.5.0)
    ##    car                3.1-3      2024-09-27 [1] CRAN (R 4.5.1)
    ##    carData            3.0-5      2022-01-06 [1] CRAN (R 4.5.1)
    ##    cli                3.6.5      2025-04-23 [1] CRAN (R 4.5.0)
    ##    cowplot          * 1.2.0      2025-07-07 [1] CRAN (R 4.5.1)
    ##    DHARMa           * 0.4.7      2024-10-18 [1] CRAN (R 4.5.0)
    ##    digest             0.6.37     2024-08-19 [1] CRAN (R 4.5.0)
    ##    dplyr            * 1.1.4      2023-11-17 [1] CRAN (R 4.5.0)
    ##    emmeans          * 1.11.1     2025-05-04 [1] CRAN (R 4.5.0)
    ##    estimability       1.5.1      2024-05-12 [1] CRAN (R 4.5.0)
    ##    evaluate           1.0.4      2025-06-18 [1] CRAN (R 4.5.0)
    ##    farver             2.1.2      2024-05-13 [1] CRAN (R 4.5.0)
    ##    fastmap            1.2.0      2024-05-15 [1] CRAN (R 4.5.0)
    ##    forcats          * 1.0.0      2023-01-29 [1] CRAN (R 4.5.0)
    ##    Formula            1.2-5      2023-02-24 [1] CRAN (R 4.5.0)
    ##    generics           0.1.4      2025-05-09 [1] CRAN (R 4.5.0)
    ##    ggbeeswarm       * 0.7.2      2023-04-29 [1] CRAN (R 4.5.0)
    ##    ggplot2          * 3.5.2      2025-04-09 [1] CRAN (R 4.5.0)
    ##    ggpubr           * 0.6.1      2025-06-27 [1] CRAN (R 4.5.1)
    ##    ggsignif           0.6.4      2022-10-13 [1] CRAN (R 4.5.1)
    ##    glmmTMB          * 1.1.11     2025-04-02 [1] CRAN (R 4.5.0)
    ##    glue               1.8.0      2024-09-30 [1] CRAN (R 4.5.0)
    ##    gridExtra        * 2.3        2017-09-09 [1] CRAN (R 4.5.1)
    ##    gtable             0.3.6      2024-10-25 [1] CRAN (R 4.5.0)
    ##    hms                1.1.3      2023-03-21 [1] CRAN (R 4.5.0)
    ##    htmltools          0.5.8.1    2024-04-04 [1] CRAN (R 4.5.0)
    ##    knitr              1.50       2025-03-16 [1] CRAN (R 4.5.0)
    ##    lattice            0.22-7     2025-04-02 [1] CRAN (R 4.5.1)
    ##    lifecycle          1.0.4      2023-11-07 [1] CRAN (R 4.5.0)
    ##    lme4             * 1.1-37     2025-03-26 [1] CRAN (R 4.5.0)
    ##    lmerTest         * 3.1-3      2020-10-23 [1] CRAN (R 4.5.0)
    ##    lubridate        * 1.9.4      2024-12-08 [1] CRAN (R 4.5.0)
    ##    magrittr           2.0.3      2022-03-30 [1] CRAN (R 4.5.0)
    ##    MASS               7.3-65     2025-02-28 [1] CRAN (R 4.5.1)
    ##    Matrix           * 1.7-3      2025-03-11 [1] CRAN (R 4.5.1)
    ##    mgcv               1.9-3      2025-04-04 [1] CRAN (R 4.5.1)
    ##    minqa              1.2.8      2024-08-17 [1] CRAN (R 4.5.0)
    ##    mvtnorm            1.3-3      2025-01-10 [1] CRAN (R 4.5.0)
    ##    NatParksPalettes * 0.2.0      2022-10-09 [1] CRAN (R 4.5.0)
    ##    nlme               3.1-168    2025-03-31 [1] CRAN (R 4.5.1)
    ##    nloptr             2.2.1      2025-03-17 [1] CRAN (R 4.5.0)
    ##    numDeriv           2016.8-1.1 2019-06-06 [1] CRAN (R 4.5.0)
    ##    pillar             1.10.2     2025-04-05 [1] CRAN (R 4.5.0)
    ##    pkgconfig          2.0.3      2019-09-22 [1] CRAN (R 4.5.0)
    ##    PNWColors        * 0.1.0      2020-06-12 [1] CRAN (R 4.5.3)
    ##    purrr            * 1.0.4      2025-02-05 [1] CRAN (R 4.5.0)
    ##    R6                 2.6.1      2025-02-15 [1] CRAN (R 4.5.0)
    ##    rbibutils          2.3        2024-10-04 [1] CRAN (R 4.5.0)
    ##    RColorBrewer       1.1-3      2022-04-03 [1] CRAN (R 4.5.0)
    ##    Rcpp               1.0.14     2025-01-12 [1] CRAN (R 4.5.0)
    ##    Rdpack             2.6.4      2025-04-09 [1] CRAN (R 4.5.0)
    ##    readr            * 2.1.5      2024-01-10 [1] CRAN (R 4.5.0)
    ##    reformulas         0.4.1      2025-04-30 [1] CRAN (R 4.5.0)
    ##    rlang              1.1.6      2025-04-11 [1] CRAN (R 4.5.0)
    ##    rmarkdown          2.29       2024-11-04 [1] CRAN (R 4.5.0)
    ##    rstatix            0.7.2      2023-02-01 [1] CRAN (R 4.5.1)
    ##    rstudioapi         0.19.0     2026-06-11 [1] CRAN (R 4.5.3)
    ##    scales             1.4.0      2025-04-24 [1] CRAN (R 4.5.0)
    ##    sessioninfo      * 1.2.4      2026-06-04 [1] CRAN (R 4.5.3)
    ##    stringi            1.8.7      2025-03-27 [1] CRAN (R 4.5.0)
    ##    stringr          * 1.5.1      2023-11-14 [1] CRAN (R 4.5.0)
    ##    tibble           * 3.3.0      2025-06-08 [1] CRAN (R 4.5.0)
    ##    tidyr            * 1.3.1      2024-01-24 [1] CRAN (R 4.5.0)
    ##    tidyselect         1.2.1      2024-03-11 [1] CRAN (R 4.5.0)
    ##    tidyverse        * 2.0.0      2023-02-22 [1] CRAN (R 4.5.0)
    ##    timechange         0.3.0      2024-01-18 [1] CRAN (R 4.5.0)
    ##  D TMB                1.9.17     2025-03-10 [1] CRAN (R 4.5.0)
    ##    tzdb               0.5.0      2025-03-15 [1] CRAN (R 4.5.0)
    ##    vctrs              0.6.5      2023-12-01 [1] CRAN (R 4.5.0)
    ##    vipor              0.4.7      2023-12-18 [1] CRAN (R 4.5.0)
    ##    withr              3.0.2      2024-10-28 [1] CRAN (R 4.5.0)
    ##    xfun               0.52       2025-04-02 [1] CRAN (R 4.5.0)
    ##    xtable             1.8-4      2019-04-21 [1] CRAN (R 4.5.0)
    ##    yaml               2.3.10     2024-07-26 [1] CRAN (R 4.5.0)
    ## 
    ##  [1] C:/Users/krdeutsch/AppData/Local/Programs/R/R-4.5.1/library
    ## 
    ##  * ── Packages attached to the search path.
    ##  D ── DLL MD5 mismatch, broken installation.
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────
