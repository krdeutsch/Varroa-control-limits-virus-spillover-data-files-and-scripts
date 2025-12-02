# Load required packages ----
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

# Read in data ----
# Honey bee colony metadata
full.hb.meta <-read_csv("FILE-PATH/honey.bee.metadata.csv")

# Honey bee colony virus loads 
HB.dwv.data <-read_csv("FILE-PATH/honey.bee.dwv.data.csv") 

# Flower metadata
full.flower.meta <- read_csv("FILE-PATH/flower.metadata.csv")

# DWV presence on flowers 
flower.dwv.data <- read_csv("FILE-PATH/flower.dwv.data.csv")

# DWV presence in wild pollinators 
poll_DWV <- read_csv("FILE-PATH/all.wild.pollinators.dwv.data.csv") 

## Create a dataframe of apiary-level averages for varroa load and honey bee visitation rates ----
# First we create a dataframe of apiary-average varroa loads from all colonies sampled for varroa at an apiary, even if honey bees from those colonies were not screened for DWV
apiary.varroa.data <-full.hb.meta %>% 
  group_by(Apiary_ID) %>% 
  summarize(varroa_avg = mean(Percent_varroa),Num_colonies=mean(Num_colonies))
# the colony density (number of colonies) is the same for every colony sampled from the apiary, so taking the mean just condenses this to a single apiary-level value

apiary.visitation.data <- full.flower.meta %>% 
# first we get a single value of honey bee visitation for each unique quadrat observation period, since multiple flowers were collected from the same quadrat - they all have the same visitation so averaging just condenses this to a single quadrat visitation value
  group_by(Apiary_ID,Sampling_Date,`Quadrat #`) %>%
  summarize(quadrat.hb.visitiation=mean(`# of HB Visits`,na.rm=TRUE)) %>% 
# then, for each apiary, we calculate an average visitation value from all unique quadrat observation periods at an apiary
    group_by(Apiary_ID) %>% 
  summarize(avg.hb.visitation = mean(quadrat.hb.visitiation,na.rm=TRUE)) %>% 
  mutate_if(is.numeric, round, 2) 

# Join these two datasets together for a dataframe that has apiary: (1) colony density; (2) average varroa load; and (3) average honey bee visitation rate
apiary.averages <- apiary.varroa.data %>% 
  full_join(apiary.visitation.data,by="Apiary_ID")

## Figure S1: variation in honey bee colony varroa levels ---- 
# This code creates a boxplot of varroa levels across the 27 apiaries sampled for this study. 
varroa_col <- c("#dd4124","#00496F")

Fig_S1 <- full.hb.meta %>% 
  left_join(apiary.averages, by="Apiary_ID") %>% 
  select(-Num_colonies.y) %>% 
  rename(Num_colonies=Num_colonies.x) %>% 
  mutate(thresh =ifelse(varroa_avg > 3, "Above","Below")) %>% 
  mutate(apiary.number=sapply(strsplit(Apiary_ID, "_"), `[`,2)) %>%
  mutate(Apiary.Name=paste("Apiary",apiary.number,sep=" ")) %>% 
  ggplot(aes(y=Percent_varroa, x=fct_reorder(Apiary.Name, Percent_varroa), fill=thresh))+
# The mean infestation level of fourteen apiaries was above the 3% threshold (red), and the mean infestation level of thirteen apiaries fell below the 3% threshold (blue). 
  geom_hline(yintercept=3, linetype="solid", color = "darkgray",lwd=2, alpha=0.75)+
# The 3% infestation threshold at which beekeepers are recommended to treat is plotted as a gray vertical line. 
  geom_boxplot(col= "black",alpha=0.6)+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, fill="black",alpha=.8,show.legend = F)+
# Each point is the varroa level in an individual colony, plotted on top of the summary boxplot 
  theme_minimal()+
  scale_fill_manual(name = "Varroa threshold",values=varroa_col)+
  coord_flip()+
  ylab("Varroa load (mites per 100 bees)")+
  theme(axis.text.x = element_text(color = "black", angle = 35, size = 22, hjust=1)) +
  theme(axis.text.y = element_text(color = "black", size = 15, hjust=1)) +
  theme(axis.title.x = element_text(color = "black", size = 25)) +
  # theme(axis.title.y = element_text(color = "black", size = 25)) +
  theme(axis.title.y = element_blank()) +
  theme(legend.background = element_rect(color = 'black', fill = 'white', linetype='solid'))+
  theme(legend.position = c(0.75, 0.25))+
  theme(legend.title = element_text(size=22))+
  theme(legend.text = element_text(size=22))
Fig_S1

# Honey Bee Analyses and Figures ----
## Calculate relative load values using ddCT methods ----
HB_loads <- HB.dwv.data %>% 
  mutate(deltaCT = DWV-`28S`) %>% 
# we ignore the DWV Cq values >35 for now (will set them to 0 in line 96), as setting them to 0 here will lead to non-finite values & nonsensical calculations
# the largest dCT value is 25.02667, so we use Apiary_16-212 as the reference sample
  mutate(deltadeltaCT = deltaCT - 25.02667) %>%
  mutate(twoddCT =  2^(-deltadeltaCT)) %>% 
  mutate(logtwoddCT =  log10(2^(-deltadeltaCT))) %>% 
# now we set all values above 35 to non-detections 
  mutate(logtwoddCT = ifelse(Sample_ID %in% c("Apiary_22-229", "Apiary_12-120","Apiary_25-218","Apiary_12-117","Apiary_19-412","Apiary_20-401"), 0, logtwoddCT)) %>% 
# we want to join each colony to its varroa load value from the metadata file
  left_join(full.hb.meta, by="Sample_ID") %>% 
  select(-c(Apiary_ID.y,Virus_Screening)) %>% 
  rename(Apiary_ID=Apiary_ID.x) %>% 
# some varroa loads are 0, which leads to non-finite logged values (necessary for fulfilling assumptions of linearity in analysis of individual colonies), so we added the smallest varroa infestation level (0.31) to ALL varroa levels
  mutate(logv = log10(Percent_varroa+0.31))

## Site-level models ----
# We chose to use an apiary-average varroa load for the main analysis (rather than each colony's individual varroa level - see lines 184-191 and Figure S2) as all other analyses for flowers and wild pollinators must be done at the apiary level
HB_site2 <-HB_loads %>% 
  select(-c(Num_colonies)) %>% 
  left_join(apiary.averages, by="Apiary_ID") %>% 
# to fulfill assumptions of linearity for subsequent models, we log-transform the average varroa load for each apiary
  mutate(logv_avg = log10(varroa_avg)) %>%
  mutate(thresh =ifelse(varroa_avg > 3, "Above Threshold","Below Threshold"))

# Here we use a LMM to test the effect of varroa level on DWV loads in honey bees, by correlating virus loads from individual colonies with apiary-level varroa averages
site_model.whole <- lmer(logtwoddCT~logv_avg+Julian_Day+Num_colonies+(1|Apiary_ID), data=HB_site2)
drop1(site_model.whole)
# Evaluate residuals for the model
resids_site<- simulateResiduals(site_model.whole)
plot(resids_site)

# Here we use a LMM to test the effect of binary apiary-level varroa control (on average, the apiary is above vs below treatment threshold) on DWV loads in honey bees
site_model.thresh <-lmer(logtwoddCT~thresh +Num_colonies+Julian_Day+(1|Apiary_ID), data=HB_site2)
drop1(site_model.thresh)
# Evaluate residuals for the model
resids_site.thresh<- simulateResiduals(site_model.thresh)
plot(resids_site.thresh)

### Figure 1A ----
# Here we construct a model-estimated best-fit line for the relationship between honey bee viral loads and (log-transformed) apiary-level varroa infestation levels 
model.emp.site <- emmip(site_model.whole,  ~logv_avg,at=list(logv_avg= c(-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3)), type="response", CIs = TRUE, plotit = FALSE)

Fig_1A <- ggplot(data=model.emp.site, aes(x=logv_avg, y=yvar)) +
  xlim(-0.45,1.4)+
  geom_line(color="#00496F",lwd=1.75) +
  geom_vline(xintercept=log10(3), linetype="dashed", color = "gray",lwd=1.5)+
  ylim(0,8)+
  geom_ribbon(aes(ymin=LCL, ymax=UCL), fill="#00496F", show.legend = FALSE,  alpha=0.35) + 
# We add in our raw data points of colony-level DWV loads for all 4 colonies from an apiary, plotted against the apiary average varroa load 
  geom_point(data=HB_site2, aes(x=logv_avg, y=logtwoddCT), alpha=0.85,size=3)+ 
  scale_fill_manual(values=varroa_col)+
  theme_minimal() +
  xlab(expression(~Average~varroa~load~per~apiary~(log[10]~mites~per~100~bees))) +  
  ylab(expression(atop(DWV~load~"in"~honey~bees, paste((log[10]~2^{"-ΔΔCT"})))))  +
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0.5)) +
  annotate(geom="text", x=1.25,y=0.15, label="p<0.001",size=5)+
  annotate(geom="text", x=-0.35,y=7.75, label="(A)",size=14)
Fig_1A

### Figure 1B ----
# Here we construct model-estimated values for honey bee viral loads from apiaries, on average, above or below the binary varroa treatment threshold 
model.emp.small.hb.model <- emmip(site_model.thresh, ~thresh, type="response", CIs = TRUE, plotit = FALSE)
#  thresh          yvar    SE df lower.CL upper.CL
# Above Threshold   4.25 0.280 23     3.67     4.83
# Below Threshold   2.01 0.291 23     1.40     2.61

varroa_col <- c("#dd4124","#00496F")

Fig_1B <- ggplot()+
# First we plot the raw data points of colony-level DWV loads for colonies from apiaries (on average) either above or below the varroa treatment threshold 
  geom_point(data=HB_site2, aes(x=thresh, y=logtwoddCT, fill=thresh),
             shape = 21,size=5, position=position_jitter(width=0.15, height=0), alpha=.4,show.legend = F)+
# We then plot model-estimated confidence intervals for honey bee DWV loads for colonies from apiaries above or below the varroa treatment threshold 
  geom_errorbar(data=model.emp.small.hb.model,aes(y=yvar, x=fct_reorder(thresh, yvar), ymin=yvar-SE, ymax=yvar+SE),
                show.legend = F,colour = "black", stat = "identity", width = 0.15,linewidth=1)+
# Finally, we plot the model-estimated honey bee DWV loads for colonies from apiaries above or below the varroa treatment threshold 
  geom_point(data=model.emp.small.hb.model,aes(y=yvar, x=fct_reorder(thresh, yvar), fill=thresh),show.legend = F,col="black",size=9) +
  geom_point(data=model.emp.small.hb.model,aes(y=yvar, x=fct_reorder(thresh, yvar), fill=thresh),show.legend = F,col=varroa_col,size=7) +
  scale_fill_manual(values=varroa_col)+
  ylab(expression(atop(DWV~load~"in"~honey~bees, paste((log[10]~2^{"-ΔΔCT"})))))  +
  xlab(expression(~Average~apiary~mite~level)) +  
  theme_minimal()+
  ylim(0,8)+
  theme(legend.position="none")+
  geom_bracket(xmin = 1, xmax = 2, y.position = 6.75, 
              label = "****",size=1,label.size = 6) +
  annotate(geom="text", x=2.4,y=0.15, label="p<0.0001",size=5)+
  annotate(geom="text", x=0.75,y=7.75, label="(B)",size=14)+
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0.5))+
  theme(axis.title.x = element_text(margin = margin(t = 9)))
Fig_1B 

### Figure 1 Whole ----
Figure.1 <- grid.arrange(Fig_1A, Fig_1B, ncol=3, widths = c(0.85,0.04, 0.6), layout_matrix = rbind(c(1,NA,2)))

## Individual colony models ----
# Now we correlate virus and varroa loads from the SAME colony 
# Here we use an LMM to test the effect of varroa level on DWV loads in honey bees at the individual colony level
individ_load_model.whole <- lmer(logtwoddCT~logv+Num_colonies +Julian_Day+(1|Apiary_ID), data=HB_loads)
drop1(individ_load_model.whole)
# Evaluate residuals for the model
resids.hb.indiv <- simulateResiduals(individ_load_model.whole)
plot(resids.hb.indiv)

### Figure S2 ----
# Here we construct a model-estimated best-fit line for the relationship between individual honey bee colony viral loads and (log-transformed) varroa infestation levels in the same colony
model.emp.indiv <- emmip(individ_load_model.whole,  ~logv,at=list(logv= c(-0.75,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6)), type="response", CIs = TRUE, plotit = FALSE)

Fig_S2 <- ggplot(data=model.emp.indiv, aes(x=logv, y=yvar)) +
  geom_line(color="#00496F",lwd=1.75) +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), fill="#00496F", show.legend = FALSE,  alpha=0.35) + 
  ylim(0,7)+  geom_vline(xintercept=log10(3), linetype="dashed", color = "gray",lwd=1.5)+
# We also plot the raw data points of colony-level DWV loads for all colonies across varroa infestation levels
  geom_point(data=HB_loads, aes(x=logv, y=logtwoddCT), size=3.5)+ 
  theme_minimal() +
  xlab(expression(Colony~varroa~load~(log[10]~mites~per~100~honey~"bees"))) +  
  ylab(expression(atop(DWV~load~"in"~honey~bees, paste((log[10]~2^{"-ΔΔCT"})))))  +
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0.5))+
  annotate(geom="text", x=1.5,y=0.25, label="p<0.001",size=6)
Fig_S2

# Flower Analyses and Figures ----
# Calculate overall DWV prevalence on goldenrod flowers across all apiaries
flower.dwv.data %>% 
  summarize(total_flowers=n(),overall_prev=(sum(DWV_presence)/n()),DWV_positives = sum(DWV_presence))
# total_flowers  overall_prev  DWV_positives
#           342         0.263             90

# Calculate apiary-level DWV prevalence on goldenrod flowers (for future visualization)
flower_prev = flower.dwv.data %>% 
  group_by(Apiary_ID) %>%
  summarize(DWV_positives = sum(DWV_presence), total_flowers=n()) %>% 
  mutate(DWV_prev=DWV_positives/total_flowers) %>% 
  mutate(DWV_negatives=(total_flowers-DWV_positives)) %>% 
  left_join(apiary.averages, by="Apiary_ID") %>% 
  mutate(thresh=ifelse(varroa_avg > 3, "Above Threshold","Below Threshold"))

## Varroa & Visitation Analyses ----
# We want to combine our data on DWV presence on flowers with our apiary-level metadata (average varroa load, honey bee visitation, colony density) and Julian day 
flower_dwv <- flower.dwv.data %>% 
  left_join(full.flower.meta, by="Sample_ID") %>%
  select(-c(Apiary_ID.y,Virus_Screening)) %>% 
  rename(Apiary_ID=Apiary_ID.x) %>% 
  left_join(apiary.averages, by="Apiary_ID") %>% 
  mutate(thresh=ifelse(varroa_avg > 3, "Above Threshold","Below Threshold")) 
View(flower_dwv)

# Here we use a GLMM to test the effect of varroa level and honey bee visitation rate on flower DWV prevalence 
flower.model.whole <-glmmTMB(DWV_presence~varroa_avg +avg.hb.visitation+Julian_Day+(1|Apiary_ID), family=binomial, data=flower_dwv)
drop1(flower.model.whole,test="Chisq")
# Evaluate residuals for the model
resids.flower <- simulateResiduals(flower.model.whole)
plot(resids.flower)

### Here we use a GLMM  to test the effect of binary apiary-level varroa control on flower DWV prevalence (above or below treatment threshold)
flower.model.thresh<- glmmTMB(DWV_presence~thresh +avg.hb.visitation+Julian_Day+ (1|Apiary_ID), family="binomial", data=flower_dwv)
drop1(flower.model.thresh, test="Chisq")
# Evaluate residuals for the model
resids.flower.thresh <- simulateResiduals(flower.model.thresh)
plot(resids.flower.thresh)

# Calculate model-estimated DWV prevalence on flowers from apiaries above vs. below varroa treatment threshold
emmeans(flower.model.thresh, ~thresh, type="response")
# thresh           prob     SE  df asymp.LCL asymp.UCL
# Above Threshold 0.284 0.0528 Inf     0.193     0.398
# Below Threshold 0.179 0.0465 Inf     0.105     0.289


### Figure 2A ----
# Here we construct a model-estimated best-fit line for the relationship between DWV prevalence on flowers and apiary-average varroa loads
model.emp.flower.v <- emmip(flower.model.whole, ~varroa_avg, at=list(varroa_avg= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)), type="response", CIs = TRUE, plotit = FALSE)

Fig.2A <-ggplot(data=model.emp.flower.v, aes(x=varroa_avg, y=yvar)) +
  geom_line(color="#00496F",lwd=1.75) +
# The 3% infestation threshold at which beekeepers are recommended to treat is plotted as a gray vertical line. 
  geom_vline(xintercept=3, linetype="dashed", color = "gray",lwd=1.5)+
  ylim(0,1)+
  geom_ribbon(aes(ymin=LCL, ymax=UCL), fill="#00496F", show.legend = FALSE,  alpha=0.35) + 
# We also plot the raw data points of apiary-level DWV prevalence on flowers against apiary-average varroa load
  geom_point(data=flower_prev, aes(x=varroa_avg, y=DWV_prev), size=5)+ 
  theme_minimal() +
  xlab(expression(atop(Average~varroa~load, paste("(mites"~per~100~"bees)"))))  +
  ylab("DWV prevalence on flowers") +  
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0.5))+
  annotate(geom="text", x=15,y=0.02, label="p=0.076",size=5)+
  annotate(geom="text", x=15,1, label="(A)",size=14)
Fig.2A

### Figure 2B ----
# Here we construct a model-estimated best-fit line for the relationship between DWV prevalence on flowers and apiary-average honey bee visitation rate
model.emp.col <- emmip(flower.model.whole, ~avg.hb.visitation, at=list(avg.hb.visitation= c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32)), type="response", CIs = TRUE, plotit = FALSE)

Fig.2B <- ggplot(data=model.emp.col, aes(x=avg.hb.visitation, y=yvar)) +
  geom_line(color="#00496F",lwd=1.75) +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), fill="#00496F", alpha=0.35, show.legend = FALSE) + 
# We also plot the raw data points of apiary-level DWV prevalence on flowers against apiary-average honey bee visitation rate
  geom_point(data=flower_prev, aes(x=avg.hb.visitation, y=DWV_prev), size=5)+ 
  theme_minimal() +
  ylim(0,1)+
  xlab(expression(atop(Average~honey~bee~visitation~rate, paste("(visits"~per~5~"min)"))))  +
  ylab("DWV prevalence on flowers") +
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0.5))+
  annotate(geom="text", x=30,y=0.02, label="p=0.004",size=5)+
  annotate(geom="text", x=30,y=1, label="(B)",size=14)
Fig.2B

## Distance Analysis ----
# We first exclude 10 flowers for which distance data was not collected, as well as create a concatenated identifier for unique quadrat x date combinations
flower_dwv_dist <- flower_dwv %>% 
  mutate(`Distance to colonies (m)` = as.numeric(`Distance to colonies (m)`)) %>% 
  mutate(`# of HB Visits`=as.numeric(`# of HB Visits`)) %>% 
  mutate(date_numeric=sapply(strsplit(Sample_ID, "_"), `[`,3)) %>% 
  mutate(quadrat_date = paste(`Quadrat #`,date_numeric, sep="_")) %>% 
  filter(!c(Apiary_ID=="Apiary_27"&date_numeric=="090720")) %>% 
# excluding flowers for which no distance data was collected (only at Apiary_27 on 09/07/2020)
  mutate(thresh=ifelse(varroa_avg > 3, "above","below"))

# Here we use a GLMM to test the effect of distance and apiary-average varroa level on the likelihood of DWV presence on flowers
flower.dist.model.whole <- glmmTMB(DWV_presence~ varroa_avg+`Distance to colonies (m)`+Julian_Day+(1|Apiary_ID/quadrat_date), family=binomial, data=flower_dwv_dist)
drop1(flower.dist.model.whole, test="Chisq")
# Evaluate residuals for the model
resids.flower.dist <- simulateResiduals(flower.dist.model.whole)
plot(resids.flower.dist)

### ### Figure 2C ----
# Here we construct a model-estimated best-fit line for the relationship between DWV presence on flowers and distance from the apiary
model.emp.flower.dist <- emmip(flower.dist.model.whole, ~`Distance to colonies (m)`, at=list(`Distance to colonies (m)`= c(0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350,400,450,500,550,600,650)), type="response", CIs = TRUE, plotit = FALSE)

Fig.2C <- ggplot(data=model.emp.flower.dist, aes(x=`Distance to colonies (m)`, y=yvar)) +
  geom_line(color="#00496F",lwd=1.75) +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), fill="#00496F", alpha=0.35, show.legend = FALSE) + 
# We also plot the raw data points of DWV presence on flowers across distances from the apiary
  geom_point(data=flower_dwv_dist, aes(x=`Distance to colonies (m)`, y=DWV_presence), alpha=.75,size=5)+ 
  theme_minimal() +
  xlim(0,650)+
  ylab("Probability of DWV on flowers") +
  xlab(expression(atop(Distance~from~apiary~(m), paste(""))))  +
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0.5))+
  annotate(geom="text", x=585,y=0.125, label="p=0.039",size=5)+
  annotate(geom="text", x=600,y=1, label="(C)",size=14)
Fig.2C

### Figure 2 Whole -----
Figure.2 <- grid.arrange(Fig.2A, Fig.2B, Fig.2C, ncol=3, widths = c(1,1,1))

# Bumble bee analyses -----
# First we select only bumble bees from our overall pollinator dataframe
bumble_dwv <- poll_DWV %>% 
  filter(Taxa=="bumble") %>% 
  left_join(apiary.averages,by="Apiary_ID") %>% 
  mutate(thresh=ifelse(varroa_avg > 3, "Above Threshold","Below Threshold")) 

# Calculate overall DWV prevalence in bumble bees across all apiaries
bumble_dwv %>% 
  summarize(total=n(),prevalence=sum(DWV_presence)/n())
# total (n)   prevalence
#       592        0.635

# Apiary-level DWV prevalence in bumble bees (for future visualization)
bumble_prev <- poll_DWV %>% 
  filter(Taxa=="bumble") %>% 
  group_by(Apiary_ID) %>% 
  summarize(DWV_positives = sum(DWV_presence), total_bumbles=n()) %>% 
  mutate(DWV_negatives = total_bumbles-DWV_positives, DWV_prev=DWV_positives/total_bumbles) %>% 
  left_join(apiary.averages,by="Apiary_ID") %>% 
  mutate(thresh=ifelse(varroa_avg > 3, "Above Threshold","Below Threshold"))

## Here we use a GLMM to test the effect of varroa levels on DWV presence in bumble bees by correlating bumble bee DWV prevalence across with apiary-average varroa loads and colony density 
bumble.model.whole <-glmmTMB(DWV_presence~varroa_avg+ Julian_Day+Num_colonies +(1|Apiary_ID), family=binomial, data=bumble_dwv)
drop1(bumble.model.whole, test="Chisq")
# Evaluate residuals for the model
resids.bumble <- simulateResiduals(bumble.model.whole)
plot(resids.bumble)

# Here we use a GLMM to test the effect of binary apiary-level varroa control (on average, the apiary is above vs below treatment threshold) on DWV presence in bumble bees
bumble.thresh.model <-  glmmTMB(DWV_presence~ thresh+Julian_Day+ Num_colonies+ (1|Apiary_ID), family=binomial, data=bumble_dwv)
drop1(bumble.thresh.model, test="Chisq")
# Evaluate residuals for the model
resids.bumble.thresh <- simulateResiduals(bumble.thresh.model)
plot(resids.bumble.thresh)

# Calculate model-estimated DWV prevalence in bumble bees from apiaries above vs. below varroa treatment threshold
emmeans(bumble.thresh.model, ~thresh, type="response")
#  thresh           prob     SE  df asymp.LCL asymp.UCL
# Above Threshold 0.798 0.0433 Inf     0.700     0.870
# Below Threshold 0.484 0.0647 Inf     0.361     0.609

### Figure 3A----
# Here we construct a model-estimated best-fit line for the relationship between bumble bee DWV prevalence and apiary-level varroa infestation levels 
model.emp.bumble <- emmip(bumble.model.whole, ~varroa_avg, at=list(varroa_avg= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)), type="response", CIs = TRUE, plotit = FALSE)

Fig.3A <- ggplot(data=model.emp.bumble, aes(x=varroa_avg, y=yvar)) +
  geom_line(color="#00496F",lwd=1.75) +
  geom_ribbon(aes(ymin=LCL, ymax=UCL), fill="#00496F", show.legend = FALSE,  alpha=0.35) +
# The 3% infestation threshold at which beekeepers are recommended to treat is plotted as a gray vertical line. 
  geom_vline(xintercept=3, linetype="dashed", color = "gray",lwd=1.5)+
  ylim(0,1)+
# We also plot the raw data points of bumble bee DWV prevalence against apiary-average varroa loads
  geom_point(data=bumble_prev, aes(x=varroa_avg, y=DWV_prev), size=5)+ 
  theme_minimal() +
  xlab("Average varroa load per apiary (mites per 100 honey bees)") +
  ylab("DWV prevalence in bumble bees") +
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0.5)) +
  annotate(geom="text", x=17,y=0.025, label="p<0.0001",size=5)+
  annotate(geom="text", x=0.5,y=1, label="(A)",size=14)+
  theme(axis.title.x = element_text(margin = margin(t = 10)))+
  theme(axis.title.y = element_text(margin = margin(r = 10)))
Fig.3A

### Figure 3B ----
# Calculate model-estimated DWV prevalence in bumble bees from apiaries above vs. below varroa treatment threshold
model.thresh.bumble <- as.data.frame(emmip(bumble.thresh.model, ~thresh, type="response", CIs = TRUE, plotit = FALSE))

varroa_col <- c("#dd4124","#00496F")

Fig.3B <- ggplot()+
# First we plot the raw data points of bumble bee DWV prevalence from apiaries either above or below the varroa treatment threshold 
  geom_point(data=bumble_prev, aes(x=thresh, y=DWV_prev, fill=thresh),
             shape = 21,size=5, position=position_jitter(width=0.15, height=0), alpha=1,show.legend = F)+
# We then plot the model-estimated prevalence of DWV in bumble bees from apiaries above or below the varroa treatment threshold 
    geom_col(data=model.thresh.bumble,aes(y=yvar, x=fct_reorder(thresh, yvar), fill=thresh),show.legend = F,col=varroa_col, alpha=0.65) +  
# Finally, we plot model-estimated confidence intervals for bumble bee DWV prevalence for apiaries above or below the varroa treatment threshold 
  geom_errorbar(data=model.thresh.bumble,aes(y=yvar, x=fct_reorder(thresh, yvar), ymin=yvar-SE, ymax=yvar+SE),show.legend = F,colour = "black", stat = "identity", lwd=1, width = 0.15)+
  scale_fill_manual(values=varroa_col)+
  ylab("DWV prevalence in bumble bees")  +
  xlab("Average apiary mite level")  +
  theme_minimal() +
  theme(legend.position="none")+
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0.5)) +
  geom_bracket(xmin = 1, xmax = 2, y.position = 1, 
               label = "***",size=1,label.size = 8) +
  annotate(geom="text", x=2.4,y=0.85, label="p<0.001",size=5)+
  annotate(geom="text", x=0.65,y=1, label="(B)",size=14)+
  theme(axis.title.x = element_text(margin = margin(t = 10)))+
  theme(axis.title.y = element_text(margin = margin(r = 10)))
Fig.3B

### Figure 3 Whole ----
Figure.3 <- grid.arrange(Fig.3A, Fig.3B, ncol=3, widths = c(0.8,0.01, 0.5), layout_matrix = rbind(c(1,NA,2)))

# Wild pollinators ----
## Prevalence model ----
# First we select all non-Bombus pollinators from our overall pollinator dataframe
poll_nb <- poll_DWV %>% 
  filter(!Taxa=="bumble") %>% 
# For statistical analyses, we grouped wild non-Bombus bees by families and hover flies by subfamilies
# Bee families with a sample size of less than 10, or bees lacking family identification, were grouped into a broader category of "Other bees"
  mutate(stat_grouping=ifelse(Family %in% c("Colletidae", "Andrenidae","Unknown"),"Other bees",
                         ifelse(Subfamily %in% c("Syrphinae","Eristalinae","unknown"),Subfamily,Family)))%>% 
  mutate(stat_grouping=ifelse(stat_grouping=="unknown","Unknown syrphid",stat_grouping)) %>% 
  relocate(stat_grouping, .after=Species) %>% 
  left_join(apiary.averages,by="Apiary_ID")  %>% 
  mutate(thresh =ifelse(varroa_avg > 3, "Above","Below"))

poll_nb %>% 
  group_by(Taxa) %>% 
  summarize(total=n(), DWV_positives=sum(DWV_presence), DWV_prevalence=(sum(DWV_presence)/n()))
#  Taxa      n() DWV_positives DWV_prevalence
# syrphid    128            28          0.219
# wild bee   116            46          0.397

# # Here we use a GLMM to test the effect of pollinator identity and binary apiary-level varroa control (on average, the apiary is above vs below treatment threshold) on DWV presence in wild non-Bombus pollinators
model_thresh <-  glmmTMB(DWV_presence~ thresh*stat_grouping+ Julian_Day+ (1|Apiary_ID), family=binomial, data=poll_nb)
drop1(model_thresh, test="Chisq")
# Evaluate residuals for the model
resids_thresh <- simulateResiduals(model_thresh)
plot(resids_thresh)

# We use the emmeans package to investigate pairwise differences in DWV prevalence between pollinators from apiaries above/below varroa treatment threshold and between pollinator taxa 
emmeans(model_thresh, pairwise ~  thresh  | stat_grouping, type="response", adjust="tukey")
emmeans(model_thresh, pairwise ~  stat_grouping, type="response")

### Figure 4A -----
# Calculate samples sizes for each taxonomic grouping collected from apiaries above/below threshold
thresh.sample_size <- poll_nb %>% 
  group_by(stat_grouping,thresh) %>% 
  summarize(total=n())

# Construct model-estimated DWV prevalence in wild non-Bombus pollinators from apiaries above vs. below varroa treatment threshold
model.thresh.nbs <- as.data.frame(emmip(model_thresh, ~thresh | stat_grouping, type="response", CIs = TRUE, plotit = FALSE)) %>% 
# this is to keep hover flies together while visualizing differences in DWV prevalence  
  mutate(order =ifelse(stat_grouping %in% c("Eristalinae","Syrphinae"),1,2)) %>% 
# this is to have apiaries "below" threshold be shown first on plots
  mutate(viz =ifelse(thresh=="Below",1,2)) %>% 
# for simplicity, we will not shown unknown syrphids (n=3) on our plot
  filter(!stat_grouping %in% c("Unknown syrphid"))%>% 
  left_join(thresh.sample_size)
  
varroa_col_reverse=c("#00496F","#dd4124")

Fig.4A <-ggplot2::ggplot(data=model.thresh.nbs, aes(x=fct_reorder(stat_grouping, order), y=yvar, ymin=yvar-SE,ymax=yvar+SE, fill=fct_reorder(thresh, viz)))+
  geom_col(position=position_dodge(width=0.95),show.legend = F, alpha=0.88) +    
  geom_errorbar(col="black",stat = "identity", lwd=0.5,width = 0.25,position=position_dodge(width=0.95))+
  ylab("DWV prevalence in wild pollinators") +
  theme_minimal()+
  scale_fill_manual(name = "Varroa threshold",values=varroa_col_reverse)+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 330, color="grey29", size = 24, hjust=0)) +
  theme(axis.title.x = element_blank())+ 
  theme(text = element_text(size=26)) +
  geom_bracket(data=model.thresh.nbs, xmin = 1.75, xmax = 2.25, y.position = 0.35, 
                                                      label = "*",size=0.75,label.size = 7) + 
  geom_bracket(data=model.thresh.nbs, xmin = 3.75, xmax = 4.25, y.position = 0.73, 
               label = "*",size=0.75,label.size = 7)+
  annotate(geom="text", x=1,y=1, label="(A)",size=14)+
  geom_text(aes(label = total),
            position = position_dodge(width = 1),
            vjust = 1, size = 5,color="white",  fontface = "bold")
Fig.4A

## DWV Load models ----
### Hover fly load model ----
# First we calculate DWV relative load values for hover flies using ddCT methods 
hf_load <- poll_DWV %>% 
  filter(Taxa=="syrphid") %>% 
  filter(!DWV_presence == 0) %>% 
  mutate(Reference_Cq=as.numeric(Reference_Cq)) %>% 
  filter(Reference_Cq <= 30 & Reference_Cq >= 15) %>% 
  mutate(deltaCT = DWV_Cq-Reference_Cq) %>% 
# the largest dCT value is 14.910000, so we use Apiary_14_090920_P2 as the reference sample
  mutate(deltadeltaCT =  deltaCT - 14.910000) %>% 
  mutate(twoddCT =  2^(-deltadeltaCT)) %>% 
  mutate(logtwoddCT= log10(twoddCT)) %>% 
  left_join(apiary.averages) %>% 
  mutate(thresh =ifelse(varroa_avg > 3, "Above","Below")) 

# Here we use a LMM to test the effect of binary apiary-level varroa control (on average, the apiary is above vs below treatment threshold) on DWV loads in hover flies
load.model.hf <- glmmTMB(logtwoddCT~ Subfamily+thresh + Julian_Day +(1|Apiary_ID), data=hf_load)
drop1(load.model.hf,  test="Chisq")
# Evaluate residuals for the model
resids.load.hf <-simulateResiduals(load.model.hf)
plot(resids.load.hf)

### Bee load model ----
bee_load <- poll_DWV %>% 
  filter(Taxa=="wild bee") %>% 
  filter(!DWV_Cq == "NaN") %>% 
# Bee families with a sample size of less than 10, or bees lacking family identification, were grouped into a broader category of "Other bees" for statistical analyses
  mutate(stat_grouping=ifelse(Family %in% c("Colletidae", "Andrenidae","Unknown"),"Other bees",Family)) %>% 
  mutate(Reference_Cq=as.numeric(Reference_Cq)) %>% 
  filter(DWV_Cq < 35) %>% 
  filter(Reference_Cq <= 30 & Reference_Cq >= 15) %>% 
  mutate(deltaCT = DWV_Cq-Reference_Cq) %>% 
# the largest dCT value is 19.928333, so we use Apiary_9_091220_P14 as the reference sample
  mutate(deltadeltaCT =  deltaCT - 19.928333) %>% 
  mutate(twoddCT =  2^(-deltadeltaCT)) %>% 
  mutate(logtwoddCT= log10(twoddCT)) %>% 
  left_join(apiary.averages) %>% 
  mutate(thresh =ifelse(varroa_avg > 3, "Above","Below"))

# Here we use a LMM to test the effect of binary apiary-level varroa control (on average, the apiary is above vs below treatment threshold) on DWV loads in wild non-Bombus bees
load.model.bee <- lmer(logtwoddCT~ stat_grouping+thresh+Julian_Day +(1|Apiary_ID), data=bee_load)
drop1(load.model.bee, test="Chisq")
# Evaluate residuals for the model
resids.load <-simulateResiduals(load.model.bee)
plot(resids.load)
# We use emmeans package to investigate pairwise differences in DWV loads between bee families 
emmeans(load.model.bee, pairwise ~  stat_grouping, type="response")

### Fig 4B ----
# Construct model-estimated DWV loads in wild (1) hover flies and (2) non-Bombus bees from apiaries above vs. below varroa treatment threshold
model.thresh.load.hf <- as.data.frame(emmip(load.model.hf, ~thresh, type="response", CIs = TRUE, plotit = FALSE)) %>% 
  mutate(taxa="syrphid")
model.thresh.load.bee <- as.data.frame(emmip(load.model.bee, ~thresh, type="response", CIs = TRUE, plotit = FALSE)) %>% 
  mutate(taxa="bee")
# Combine model-estimated DWV loads for hover flies and non-Bombus bees into a single dataframe 
load.estimates.full <- rbind(model.thresh.load.bee, model.thresh.load.hf) %>% 
  mutate(taxa=ifelse(taxa=="bee","Bee", "Hover fly")) %>% 
# this line is silly but is included so apiaries "below" threshold are shown first on plots
  mutate(thresh =ifelse(thresh=="Below","A-Below","Above"))

# Combine raw DWV loads in non-Bombus bees and hover flies into a single dataframe
bee_load.fig <- bee_load %>% 
  select(-stat_grouping)
raw.load.full<- rbind(hf_load, bee_load.fig) %>% 
# again, this line is silly but is included so apiaries "below" threshold are shown first on plots
  mutate(thresh =ifelse(varroa_avg > 3, "Above","A-Below")) %>% 
  mutate(Taxa=ifelse(Taxa=="wild bee","Bee", "Hover fly"))

Fig.4B <- ggplot()+
# First we plot the raw DWV loads for wild non-Bombus pollinators collected from apiaries (on average) either above or below the varroa treatment threshold 
  ggbeeswarm::geom_quasirandom(data=raw.load.full, aes(x=Taxa, y=logtwoddCT,
       col=thresh, fill=thresh), size=5, col="black",alpha=0.5, show.legend = T,shape = 21, dodge.width = 0.8)+
# We then plot model-estimated confidence intervals for wild non-Bombus pollinator DWV loads from apiaries above or below the varroa treatment threshold 
  geom_errorbar(data=load.estimates.full, aes(x=taxa, y=yvar, ymin=yvar-SE,ymax=yvar+SE,fill=thresh),col="black",stat = "identity", lwd=1,width = 0.35,position=position_dodge(width=0.8),show.legend = F)+
# Finally, we plot the model-estimated DWV loads for wild non-Bombus pollinators collected from apiaries above or below the varroa treatment threshold   
  geom_point(data=load.estimates.full, aes(x=fct_reorder(taxa, yvar), y=yvar, col=thresh, fill=thresh), size=9, col="black",show.legend = F,position=position_dodge(width=0.8))+
  geom_point(data=load.estimates.full, aes(x=fct_reorder(taxa, yvar), y=yvar, col=thresh, fill=thresh), size=7, position=position_dodge(width=0.8),show.legend = T)+
  scale_color_manual(name="Varroa threshold", values=varroa_col_reverse,labels = c("Below", "Above"),)+
  scale_fill_manual(name="Varroa threshold",values=varroa_col_reverse,labels = c("Below", "Above"))+
  ylab(expression(atop(DWV~load~"in"~wild~pollinators, paste((log[10]~2^{"-ΔΔCT"})))))  +
  theme_minimal()+
  theme(legend.title = element_text(size=26))+
  theme(legend.text = element_text(size=22))+
  theme(axis.text.x = element_text(angle = 330, color="grey29", size = 24, hjust=0)) +
  theme(axis.title.x = element_blank())+ 
  theme(text = element_text(size=26))+
  annotate(geom="text", x=0.75,y=8, label="(B)",size=14)
Fig.4B

### Figure 4 Whole ----
Figure.4 <-plot_grid(Fig.4A,NULL, Fig.4B, align="h",axis = "b", ncol=3,rel_widths = c(1, 0.01, 0.75))

## Figure S3: Julian Date & Bee DWV loads ----
### Figure S3A ----
# Construct model-estimated DWV loads in wild non-Bombus bees across bee taxa
load.model.bee.fam <- as.data.frame(emmip(load.model.bee, ~stat_grouping, type="response", CIs = TRUE, plotit = FALSE)) %>% 
# Indicate significance based on Table S3
    mutate(sig.fam=ifelse(stat_grouping=="Megachilidae","b","a"))

bee_fam_col <- pnw_palette("Starfish",4)

Fig.S3A<- bee_load %>% 
# For the purpose of data visualization, we treat the reference sample used to calculate load as "0" - otherwise it is a negative value  
  mutate(logtwoddCT=ifelse(logtwoddCT<0,0,logtwoddCT)) %>% 
  ggplot()+
# First plot raw DWV loads across different non-Bombus bee taxa  
  geom_point(aes(x=stat_grouping, y=logtwoddCT,fill=stat_grouping), shape = 21,size=5, position=position_jitter(width=0.15, height=0), alpha=.4,show.legend = F)+
# Then plot model-estimated confidence intervals for DWV loads across non-Bombus bee taxa
  geom_errorbar(data=load.model.bee.fam, aes(x=stat_grouping, y=yvar, ymin=yvar-SE, ymax=yvar+SE),show.legend = F,colour = "black", stat = "identity", width = 0.15,linewidth=1)+
# Finally, plot model-estimated DWV loads for different non-Bombus bee taxa   
  geom_point(data=load.model.bee.fam, aes(x=stat_grouping, y=yvar),show.legend = F,col="black",size=9)+
  geom_point(data=load.model.bee.fam, aes(x=stat_grouping, y=yvar),show.legend = F,col=bee_fam_col,size=7) +
  scale_fill_manual(values=bee_fam_col) + 
  ylab(expression(atop(DWV~load~"in"~wild~bees, paste((log[10]~2^{"-ΔΔCT"})))))  +
  theme_minimal()+
  ylim(0,8)+
  geom_text(data=load.model.bee.fam,aes(x=stat_grouping, y=yvar,label = sig.fam), position = position_dodge(width=0), vjust=-6, size=8, fontface="bold")+ 
  annotate(geom="text", x=0.85,y=8, label="(A)",size=14)+
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0,angle=330))+
  theme(axis.title.x = element_blank())
Fig.S3A

### Figure S3B -----
# Construct model-estimated best-fit line for the relationship between DWV loads in wild non-Bombus bees across time (Julian day)
load.model.bee.julian <- emmip(load.model.bee, ~Julian_Day, at=list(Julian_Day= c(245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,
                                                                          262,263,264,265,266,267,268,269,270,271,272,273,274,275,276,278)), type="response", CIs = TRUE, plotit = FALSE)
Fig.S3B <-  bee_load %>% 
# For the purpose of data visualization, we treat the reference sample used to calculate load as "0" - otherwise it is a negative value  
  mutate(logtwoddCT=ifelse(logtwoddCT<0,0,logtwoddCT)) %>% 
  ggplot()+
# Plot raw DWV loads for bees across time, sorted by taxonomic group  
  geom_point(aes(x=Julian_Day, y=logtwoddCT, col=stat_grouping),alpha=0.85,size=4.5)+
  geom_line(data=load.model.bee.julian, aes(x=Julian_Day, y=yvar),color="#00496F",lwd=1.75) +
# Plot the model-estimated best-fit line for the relationship between DWV loads in wild non-Bombus bees across Julian day
  geom_ribbon(data=load.model.bee.julian, aes(x=Julian_Day, y=yvar,ymin=LCL, ymax=UCL), fill="#00496F", show.legend = FALSE,  alpha=0.35) + 
  theme_minimal()+
  xlim(245,280)+
  scale_color_manual(values=bee_fam_col, name="Bee family")+
  xlab("Julian Day")+
  theme(axis.title.y = element_blank())+
#  ylab(expression(DWV~load~"in"~wild~bees~(log[10]~2^{"-ΔΔCT"})))  +
  theme(text = element_text(size=25), axis.text.x = element_text(hjust=0.5))+
  annotate(geom="text", x=278,y=0.15, label="p=0.046",size=6)+
  annotate(geom="text", x=247.5,y=8, label="(B)",size=14)+
  theme(axis.title.x = element_text(margin = margin(t = 9)))+
  theme(legend.position = c(0.78, 0.85))+
  theme(legend.background = element_rect(color = "black", fill = 'white'))
Fig.S3B

### Figure S3 Whole ----
Figure.S3 <-plot_grid(Fig.S3A,NULL, Fig.S3B, align="h",axis = "b", ncol=3,rel_widths = c(0.75, 0.1, 1))
Figure.S3

