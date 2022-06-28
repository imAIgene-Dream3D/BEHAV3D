library(dplyr)
library(scales)
library(MESS)
##Import dataframe live_deadROI7 for each experiment. Combine several if necessary
Combi<-readRDS("Full_well_death_dynamics")


## Rescale the mean dead cell dye per experiment
Combi<-Combi%>% 
  group_by(exp) %>% 
  mutate(red_sc=rescale(red, to=c(0,100)))%>%ungroup()
## for each well normalize the dead cell dye mean intensity to timepoint 1
Combi<-Combi%>%group_by(Org, type, exp)%>%arrange(Time2)%>%mutate(red_norm=red_sc - first(red_sc))%>%filter(Time2<23)

## calculate area under the curve for each donor/T cell co-culture
Area_under_Curve<-Combi%>%
  group_by(Org, type, exp)%>%summarise(auc=auc(Time2,red_norm))%>%ungroup()
## calculate delta increase for each donor/T cell co-culture
delta <-Combi%>%group_by(Org, type, exp)%>%arrange(Time2)%>%summarize(delta= last(red_norm) - first(red_norm))%>%ungroup()
## combine both delta and AUC in one dataframe
AUC_delta<-left_join(Area_under_Curve,delta,by=c("Org","type", "exp"))
library(xlsx)
## export calculation for plot 
write.xlsx(AUC_delta, file = "AUC_delta.xlsx")
