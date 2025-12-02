
# Varroa control limits virus spillover: data files and scripts


[Upon publication, data will be deposited in  public repository on Dryad]

This repository contains supporting files and code associated with the submitted manuscript entitled "Controlling Varroa destructor infestations in managed honey bee colonies limits virus spillover into wild pollinator communities." Three csv files provide qPCR data for DWV screens of (1) honey bees, (2) goldenrod flowers, and (3) all wild pollinators. There are additionally two metadata csv files that characterize (1) data on apiary size (number of colonies) and varroa loads from all colonies sampled (even if not screened for DWV) and (2) data  for all floral quadrats on honey bee visitation and distance from the apiary.


## Description of the data and file structure

### Data files: 

#### (1) **honey.bee.metadata.csv** : varroa loads in each colony sampled (4 or 10 colonies per apiary, based on the size of the apiary)
    - Sample_ID : unique identification code for each colony sampled 
    - Sampling_Date : the date each colony was sampled to measure varroa loads
    - Julian_Day : the 2020 Julian day for the date the sample was collected
    - Apiary_ID : identification code for the apiary at which the sample was collected
    - Num_colonies : the number of colonies in the apiary where the colony was located
    - Percent_varroa : the number of mites per 100 bees (percent infestation level) in the colony
    - Virus_Screening : denotes whether the colony was included in subsequent DWV screening (Yes = bees from the colony were screened for DWV, No = bees from the colony were not screened for DWV)
*Note that only a subset of colonies included in this dataset were selected for DWV screening, but data from all colonies were included to calculate apiary-average varroa levels.*

   
#### (2) **honey.bee.dwv.data.csv** : qPCR loads of DWV and reference gene in honey bee samples 
    - Sample_ID : unique identification code for each colony sampled 
    - 28S : average Cq value of the 28S reference gene
    - DWV : average Cq value of DWV 
    - Apiary_ID : identification code for the apiary at which the sample was collected


#### (3) **flower.metadata.csv** : information on all flowers collected for the duration of the project
    - Sample_ID : unique identification code for each flower collected
    - Apiary_ID : identification code for the apiary at which the sample was collected
    - Sampling_Date : the date the flower was collected and observed
    - Julian_Day : the 2020 Julian day for the date the flower  was collected
    - `Quadrat #` : below the site level, this denotes the quadrat in which the flower was collected (4-10 quadrats per sampling date)
    - `# of HB Visits` :  the number of honey bee visits to goldenrod flowers within the quadrat from which the flower was collected, during a five-minute period (NA = no honey bee observation data was collected from the quadrat)
    - `Distance to colonies (m)` : the distance of the quadrat in which the flower was collected from the nearest honey bee colony in the apiary (in meters)
    - Virus_Screening : denotes whether the flower was included in subsequent DWV screening (Yes = flowers were screened for DWV, No = flowers were not screened for DWV)
*Note that only a subset of flowers included in this dataset were selected for DWV screening, but data from all quadrats were included to calculate apiary-average honey bee visitation estimates.*

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
    - Julian_Day : the 2020 Julian day for the date the sample was collected
    - DWV_Cq : the Cq value from qPCR screening for the DWV target in the sample (average value from 3 technical replicates for syrphid flies and wild bees, single value for bumble bees; NA denotes non-detections)
    - Reference_Cq : the Cq value from qPCR screening for the reference target in the sample (average value from 3 technical replicates for syrphid flies and wild bees; bumble bees and two Xylocopa virginica
 were not screened for a reference gene and this is denoted "none")
    - DWV_presence : the presence or absence of DWV in each sample (0 = absent, 1 = present; Cq values above 35 were treated as non-detections and coded as "0")

## Code/Software
All statistical tests were performed using R Statistical Software (v4.3.1, R Core Team, 2023) using the packages lme4 (Bates et al., 2015), emmeans (Lenth, 2023), and glmmTMB (Brooks et al., 2017). All model fits were assessed using the DHARMa package (Hartig, 2022). The code to reproduce all the statistical tests and figures included in this manuscript can be found in **varroa.virus.spillover.full.script.R**, which also includes the full list of packages used for data carpentry and visualization.
