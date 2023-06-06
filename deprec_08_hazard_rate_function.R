
# env settings ------------------------------------------------------------

library(survival)
library(survminer)
library(openxlsx)
library(muhaz)
library(ggplot2)
library(cowplot)
library(dplyr)
load("01/tidy_data_diagnosis_after_attending.RData")

tidy_data_dia_after_att <- transform(
  tidy_data_dia_after_att,
  futime = as.numeric(Date_end - Date_of_diagnosis)
)
# calculating -------------------------------------------------------------

smoothhazp <- tidy_data_dia_after_att[, c("Platelet400", "futime", "OS")] %>%
  na.omit() %>% 
  group_by(Platelet400) %>% #根据自己的需要修改分组列名
  do(haz = muhaz(times = .$futime,delta = .$OS,
                 min.time = min(.$futime[.$OS==1]),
                 max.time = max(.$futime[.$OS==1]),
                 bw.grid = 7,bw.method = "g",b.cor = "b")) %>%
  do(data.frame(Hazard = .$haz$haz.est, 
                Days = .$haz$est.grid, 
                Subgroup = .$Platelet400))

orange <- "#EB292A"
blue <- "#2271B4"

p1 <- ggplot(smoothhazp, aes(x=Days, y=Hazard, colour= Subgroup)) +
  geom_line() + ggtitle("Kernel-smoothing hazard function plot")+
  #scale_color_brewer(palette = "Set1") + 
  scale_color_manual(values = c(blue, orange)) + 
  theme_classic() +                                                      
  theme(plot.title = element_text(hjust=0.5,size=8,vjust=0.2,face = "bold"), 
        axis.text = element_text(size=8,
                                 ,colour="black"),  
        axis.title = element_text(size=8,
                                  ,colour="black"), 
        legend.title = element_blank(),                                 
        legend.position = c(0.8,0.9)) #+
  #scale_x_continuous(name = "",breaks = seq(0,72,12)) 

p1
