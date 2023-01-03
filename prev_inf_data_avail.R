#Data availability#Figure 1

library(tidyverse)
library(reshape2)
library(data.table)
library(stringi)
library(ggplot2)

#PREVIOUS INFECTION by variant and outcome

#Fig1 Pan A: All studies
ve1 <- read.csv("FILE_PATH/pinf_final.csv")

ve2 <- subset(ve1, select = c(study_id, author, location_id, location_name, piv, prim_inf_var, variant, symptom_severity, severity, 
                              sample_size, efficacy_mean, efficacy_lower, efficacy_upper, pinf, age_start, age_end, start_interval, end_interval,
                              quality_rating, adjustment))

#creating a new varaible to have 3 categories of severity: infection, symptomatic and severe
ve2 <- mutate(ve2, sev_severity = if_else(symptom_severity == "Severe", "severe",
                                  if_else(symptom_severity == "severe", "severe",
                                  if_else(symptom_severity == "Infection", "infection",
                                  if_else(symptom_severity == "infection", "infection",
                                  if_else(symptom_severity == "asymptomatic, Mild, Moderate, Severe", "symptomatic", 
                                  if_else(symptom_severity == "Asymptompatic, Mild, Moderate + Severe", "symptomatic",
                                  if_else(symptom_severity == "Mild, Moderate", "symptomatic",
                                  if_else(symptom_severity == "Mild, Moderate + Severe", "symptomatic",
                                  if_else(symptom_severity == "Moderate, Severe", "symptomatic", "NA"))))))))))

#rename prim_variant_name primary variant name
ve2 <- mutate(ve2, prim_variant_name = if_else(prim_inf_var == "Ancestral", "Ancestral",
                                       if_else(prim_inf_var == "all variant", "Ancestral",          
                                       if_else(prim_inf_var == "Ancestral, B.1.1.7", "Ancestral",   
                                       if_else(prim_inf_var == "mixed variant", "Ancestral",        
                                       if_else(prim_inf_var == "mixed", "Ancestral",                
                                       if_else(prim_inf_var == "B.1.1.7", "Alpha",              
                                       if_else(prim_inf_var == "B.1.351", "Beta", 
                                       if_else(prim_inf_var == "B.1.617.2", "Delta",            
                                       if_else(prim_inf_var == "B.1.1.529", "Omicron BA.1",
                                       if_else(prim_inf_var == "B.1.1.529.1", "Omicron BA.1",            
                                       if_else(prim_inf_var == "B.1.1.529.2", "Omicron BA.2", "NA"))))))))))))

ve2 <- subset(ve2, prim_variant_name == "Ancestral" | prim_variant_name == "Alpha" | prim_variant_name == "Beta" |
                prim_variant_name == "Delta")

#rename variant (subsequent variant) variable
ve2 <- mutate(ve2, variant_name = if_else(variant == "Ancestral", "Ancestral",
                                  if_else(variant == "all variant", "Ancestral",           #original: "all variant", "all variant",
                                  if_else(variant == "Ancestral, B.1.1.7", "Ancestral",    #original: "Ancestral, B.1.1.7", "mixed",       	
                                  if_else(variant == "B.1.1.7", "Alpha",
                                  if_else(variant == "mixed", "Ancestral",                 #original: "mixed", "mixed",
                                  if_else(variant == "mixed variant", "Ancestral",         #original: "mixed variant", "mixed",
                                  if_else(variant == "B.1.351", "Beta",
                                  if_else(variant == "B.1.617.2", "Delta",
                                  if_else(variant == "B.1.1.529", "Omicron BA.1",
                                  if_else(variant == "B.1.1.529.1", "Omicron BA.1",
                                  if_else(variant == "B.1.1.529.2", "Omicron BA.2", 
                                  if_else(variant == "B.1.1.529.4 and B.1.1.529.5", "Omicron BA.4/BA.5", "NA")))))))))))))


ve2 <- subset(ve2, variant_name == "Ancestral" | variant_name == "Alpha" | variant_name == "Beta" |
                variant_name == "Delta" | variant_name == "Omicron BA.1")

ve2a <- distinct(ve2, study_id, sev_severity , variant_name,   )
ve2_t <- aggregate(ve2a$study_id, by=list(ve2a$sev_severity,ve2a$variant_name), FUN=length)

#Fig1 Pan A 
gg1 <- ggplot(ve2_t, aes(x = factor(Group.1, level = c("infection", "symptomatic", "severe")), y = Group.2,
                         fill = factor(x), color = factor(x), height=0.5, width = 0.9)) + #, width = 0.5
  geom_tile(colour="white",size=0.1) +
  scale_fill_manual(values = c("1" = "lightgoldenrodyellow", "2" = "lightgoldenrod2", "3" = "lightgoldenrod", "4" = "orange", 
                               "5" = "orange2", "6" = "orange3", "7" = "red", "10" = "red1", "11" = "red2", "12" = "red3", 
                                "31" = "red4")) + 
  labs(x= "Outcome", y= " ", fill = "Number of studies") +
  theme_grey(base_size=14)+
  theme(legend.text=element_text(face="bold"), 
        plot.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  scale_x_discrete(labels=c("Reinfection", "Symptomatic", "Severe")) +
#legend.position = "bottom")
  ggtitle("Panel A: All studies") + 
  theme(plot.title.position = "plot")
print(gg1)


#Fig1 Pan B: Studies with information on time since reinfection
ve1 <- read.csv("FILE_PATH/pinf_final.csv")

ve2 <- subset(ve1, select = c(exclude, study_id, author, location_id, location_name, piv, prim_inf_var, variant, symptom_severity, severity, 
                              sample_size, efficacy_mean, efficacy_lower, efficacy_upper, pinf, average_time_since_infection,
                              start_interval, end_interval))

ve2 <- mutate(ve2, sev_severity = if_else(symptom_severity == "Severe", "severe",
                                          if_else(symptom_severity == "severe", "severe",
                                                  if_else(symptom_severity == "Infection", "infection",
                                                          if_else(symptom_severity == "infection", "infection",
                                                                  if_else(symptom_severity == "asymptomatic, Mild, Moderate, Severe", "symptomatic", 
                                                                          if_else(symptom_severity == "Asymptompatic, Mild, Moderate + Severe", "symptomatic",
                                                                                  if_else(symptom_severity == "Mild, Moderate", "symptomatic",
                                                                                          if_else(symptom_severity == "Mild, Moderate + Severe", "symptomatic",
                                                                                                  if_else(symptom_severity == "Moderate, Severe", "symptomatic", "NA"))))))))))


#rename prim_variant_name primary variant name
ve2 <- mutate(ve2, prim_variant_name = if_else(prim_inf_var == "Ancestral", "Ancestral",
                                               if_else(prim_inf_var == "all variant", "Ancestral",          
                                                       if_else(prim_inf_var == "Ancestral, B.1.1.7", "Ancestral",   
                                                               if_else(prim_inf_var == "mixed variant", "Ancestral",        
                                                                       if_else(prim_inf_var == "mixed", "Ancestral",                
                                                                               if_else(prim_inf_var == "B.1.1.7", "Alpha",              
                                                                                       if_else(prim_inf_var == "B.1.351", "Beta", 
                                                                                               if_else(prim_inf_var == "B.1.617.2", "Delta",            
                                                                                                       if_else(prim_inf_var == "B.1.1.529", "Omicron BA.1",
                                                                                                               if_else(prim_inf_var == "B.1.1.529.1", "Omicron BA.1",            
                                                                                                                       if_else(prim_inf_var == "B.1.1.529.2", "Omicron BA.2", "NA"))))))))))))




ve2 <- subset(ve2, prim_variant_name == "Ancestral" | prim_variant_name == "Alpha" | prim_variant_name == "Beta" |
                prim_variant_name == "Delta")

#rename variant (subsequent variant) variable
ve2 <- mutate(ve2, variant_name = if_else(variant == "Ancestral", "Ancestral",
                                          if_else(variant == "all variant", "Ancestral",           #original: "all variant", "all variant",
                                                  if_else(variant == "Ancestral, B.1.1.7", "Ancestral",    #original: "Ancestral, B.1.1.7", "mixed",       	
                                                          if_else(variant == "B.1.1.7", "Alpha",
                                                                  if_else(variant == "mixed", "Ancestral",                 #original: "mixed", "mixed",
                                                                          if_else(variant == "mixed variant", "Ancestral",         #original: "mixed variant", "mixed",
                                                                                  if_else(variant == "B.1.351", "Beta",
                                                                                          if_else(variant == "B.1.617.2", "Delta",
                                                                                                  if_else(variant == "B.1.1.529", "Omicron BA.1",
                                                                                                          if_else(variant == "B.1.1.529.1", "Omicron BA.1",
                                                                                                                  if_else(variant == "B.1.1.529.2", "Omicron BA.2", 
                                                                                                                          if_else(variant == "B.1.1.529.4 and B.1.1.529.5", "Omicron BA.4/BA.5", "NA")))))))))))))



ve2 <- subset(ve2, variant_name == "Ancestral" | variant_name == "Alpha" | variant_name == "Beta" |
                variant_name == "Delta" | variant_name == "Omicron BA.1")

#keep week after infection studies (1: only time since infection and 2:average time since infection)
ve2 <- subset(ve2, pinf == 1 | pinf == 2) 

ve2a <- distinct(ve2, study_id, sev_severity , variant_name,   )
ve2_t <- aggregate(ve2a$study_id, by=list(ve2a$sev_severity,ve2a$variant_name), FUN=length)

#Fig1 Panel B
gg1 <- ggplot(ve2_t, aes(x = factor(Group.1, level = c("infection", "symptomatic", "severe")), y = Group.2,
                         fill = factor(x), color = factor(x), height=0.5, width = 0.9)) + #, width = 0.5
  geom_tile(colour="white",size=0.1) +
  scale_fill_manual(values = c("1" = "lightgoldenrodyellow", "2" = "lightgoldenrod2", "3" = "lightgoldenrod", "4" = "orange", 
                              "5" = "orange2","7" = "red", "15" = "darkred")) + 
  labs(x= "Outcome", y= " ", fill = "Number of studies") +
  theme_grey(base_size=14)+
  theme(legend.text=element_text(face="bold"), 
        plot.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  scale_x_discrete(labels=c("Reinfection", "Symptomatic", "Severe")) +
  #legend.position = "bottom")
  ggtitle("Panel B: Studies with information on time since reinfection") + 
  theme(plot.title.position = "plot")
print(gg1)


