library(ggplot2)
library(gridExtra)
library(dplyr)
live_deadROI6_combi<-readRDS(file = "Individual_orgs_death_dynamics")  ### set here the direction of your df with processed individual organoids death dynamics

live_deadROI6_combi$dead<-ifelse(live_deadROI6_combi$red<7,0,1)
library(tidyr)
### fill missing values with NA
live_deadROI6_combi<-live_deadROI6_combi %>%group_by(Org, exp,type)%>%arrange(Time, .by_group = TRUE)%>% complete(Time, nesting(Track2))
## refill NA with previous value
live_deadROI6_combi<-live_deadROI6_combi %>%group_by(Org, exp,type, Track2)%>%arrange(Time, .by_group = TRUE)%>%fill(names(live_deadROI6_combi))  
## select the dead cells
live_deadROI4_dead <-live_deadROI6_combi%>%group_by(Track2,exp)%>%filter(max_red>7 )
live_deadROI4_dead$track_exp<-with(live_deadROI4_dead , interaction(Track2,  exp, Time))
live_deadROI4_live <-live_deadROI6_combi%>%group_by(Track2,exp)%>%filter(max_red<=7 ) ## select the live cells
temp2 <-live_deadROI4_dead%>%group_by(Track2, exp)## group by track ID
temp2 <- temp2%>%arrange(Time, .by_group = TRUE) 
live_deadROI5_dead <- temp2 %>% group_by(Track2,exp) %>% filter(row_number() <= which.max(dead)) ## stop at max dead
test1 <-live_deadROI4_dead%>%filter(!track_exp%in% live_deadROI5_dead$track_exp)  ###find all the rows that are beyond max red
test1$dead<-1 ## set as dead beyond max red
live_deadROI6_combi_2<-rbind(live_deadROI4_live,live_deadROI5_dead[,c(-16)],test1[,c(-16)])  ## bind dataframe again
## calculate the number of organoids at T0 to calculate the percentage
live_deadROI6_n_t0<-live_deadROI6_combi_2%>%group_by(Time, Org, exp,type)%>%summarise(n_first= n())%>%filter(Time==1)
live_deadROI6_combi_2<-left_join(live_deadROI6_combi_2,live_deadROI6_n_t0, by=c("Org", "exp","type"))
live_deadROI6_combi_2$n_first<-as.numeric(live_deadROI6_combi_2$n_first)
# for each dead calculate the percentage of dead cells, if no, then 0 
live_deadROI6_per<-live_deadROI6_combi_2%>%group_by(Time.x, Org, exp,type, .drop=FALSE)%>%summarize(n=sum(dead==1), n_first=mean(n_first))%>%mutate(perc = n*100 / n_first)



library(ggplot2)
live_deadROI6_per$Time<-live_deadROI6_per$Time.x/2
### normalize to the initial number of dead organoids
perc_dead4<-live_deadROI6_per%>%group_by(Org, type, exp)%>%arrange(Time)%>%mutate(perc_dead_norm=perc - first(perc))%>%filter(Time<=24)
p1 <-ggplot(perc_dead4, aes(Time, perc_dead_norm, group=Org, color = Org)) + 
  geom_smooth(method = "loess",size = 0.5, se = T, alpha=0.3, span=1)+
  theme_bw() + 
  ylab("percentage of dead organoids") + coord_cartesian(ylim=c(0,100))+scale_y_continuous(limit=c(0,100),oob=squish)+
  xlab("Time (hours)") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=5), axis.title.y = element_text(size = 10), axis.title.x = element_text(size = 15), legend.text=element_text(size= 10))+
  ggtitle("percenatge of dead organoids overtime n0-n4")+theme(aspect.ratio=1)

p1<-p1+facet_grid(.~type, scales = "free")
p1

