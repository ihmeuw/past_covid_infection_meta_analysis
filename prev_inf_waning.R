#Previous infection waning immunity (time since infection)

library(tidyverse)
library(reshape2)
library(data.table)

ve1 <- read.csv("FILE_PATH/vaceff/pinf_final.csv")

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


#new study_id_2 variable (to differentiate studies with the same study_id but different variants)
ve2$study_id2 = paste(ve2$study_id,ve2$prim_variant_name, ve2$variant_name,  sep = "&")

#keep week after infection studies (1: only time since infection and 2:average time since infection)
ve2 <- subset(ve2, pinf == 1 | pinf == 2) 

#number of studies and countries
length(unique(ve2$study_id))
length(unique(ve2$location_name))
table(ve2$location_name)
table(ve2$study_id)

#weeks after previous infection
ve2$end_interval <- as.numeric(as.character(ve2$end_interval))
ve2$start_interval <- as.numeric(as.character(ve2$start_interval))
ve2 <- mutate(ve2, mid_point1 = ((end_interval - start_interval)/2) + start_interval)

#combine in "mid_point" two columns (mid_point1 and average_time_since_infection)
ve2 <- mutate(ve2, mid_point= if_else(is.na(ve2$mid_point1), as.numeric(ve2$average_time_since_infection), as.numeric(ve2$mid_point1)))

#run the model using logit space
ve2$efficacy_mean <- as.numeric(as.character(ve2$efficacy_mean))
ve2$efficacy_upper <- as.numeric(as.character(ve2$efficacy_upper))
ve2$efficacy_lower <- as.numeric(as.character(ve2$efficacy_lower))
ve2 <- mutate(ve2, se = (efficacy_upper - efficacy_lower)/3.92)

#transform to logit
library(crosswalk002, lib.loc = "/ihme/code/mscm/Rv4/packages/")
logit <- delta_transform(mean = ve2$efficacy_mean, sd = ve2$se, transformation = "linear_to_logit")

names(logit) <- c("mean_logit", "sd_logit")
vacc_logit <- cbind(ve2, logit)
ve2 <- vacc_logit

#meta-regression spline models with separate models by vaccine and outcome
library(dplyr)
library(mrbrt002, lib.loc = "/ihme/code/mscm/Rv4/packages/")

#Previous infection
infec <- subset(ve2, sev_severity == "infection")
#symptomatic <- subset(ve2, sev_severity == "symptomatic")
#severe <- subset(ve2, sev_severity == "severe")

#INFECTION
#all subsequent variant data for Ancestral, Alpha, Delta (#Figure 3 Panel A)
inf <- infec[!(infec$variant_name == "Omicron BA.1"), ]
inf <- inf[!(inf$variant_name == "Omicron BA.2"), ]
inf <- inf[!(inf$variant_name == "Beta"), ]

#all subsequent variant = Omicron BA.1 (#Figure 3 Panel B)
#inf <- subset(infec, variant_name == "Omicron BA.1")

#all subsequent variant BA.2 (#Figure 3 Panel B1)
#inf <- subset(infec, prim_variant_name == "Ancestral" & variant_name == "Omicron BA.2")

#SYMPTOMATIC
#all subsequent variant data for Ancestral, Alpha, Delta (#Figure 3 Panel C)
#symp <- symptomatic[!(symptomatic$variant_name == "Omicron BA.1"), ]
#symp <- symp[!(symp$variant_name == "Omicron BA.1"), ]
#symp <- symp[!(symp$variant_name == "Omicron BA.2"), ]

#all subsequent variant = Omicron BA.1 (#Figure3 Panel D)
#symp <- subset(symptomatic, variant_name == "Omicron BA.1")

#SEVERE
#all subsequent variant data for Ancestral, Alpha, Delta (#Figure 3 Panel E)
#sev <- severe[!(severe$variant_name == "Omicron BA.1"), ]
#sev <- sev[!(sev$variant_name == "Omicron BA.2"), ]

#all subsequent variant = Omicron BA.1 (#Figure 3 Panel F)
#sev <- subset(severe, variant_name == "Omicron BA.1")

#Saving data points for each plot
#write.csv(inf, "FILE_PATH/past_infection.csv", row.names = F)
#write.csv(inf, "FILE_PATH/past_infection_o.csv", row.names = F)
#write.csv(inf, "FILE_PATH/past_infection_oba2.csv", row.names = F)
#write.csv(symp, "FILE_PATH/past_symptomatic.csv", row.names = F)
#write.csv(symp, "FILE_PATH/past_symptomatic_o.csv", row.names = F)
#write.csv(sev, "FILE_PATH/pinf_curves/past_severe.csv", row.names = F)
#write.csv(sev, "FILE_PATH/pinf_curves/past_severe_o.csv", row.names = F)

#3.5.1 - Setting priors and shape constraints on splines
ve6 <- MRData()
ve6$load_df(
  data = inf,  #change to inf, symp or sev according each curve
  col_obs = "mean_logit",    
  col_obs_se = "sd_logit",   
  col_covs = list("mid_point"),
  col_study_id = "study_id2" )

mod1 <- MRBRT(
  data = ve6,
  cov_models = list(
    LinearCovModel("intercept", use_re = TRUE),  #TRUE #FALSE: Reduce the variation of the random effect #use_re  = use random effect
    LinearCovModel(
      alt_cov = "mid_point",
      use_spline = TRUE,
      spline_knots = array(seq(0, 1, length.out = 8)),  #"length.out" change according each group of data points: 6 if mid_point goes up to
                                                        #~60 weeks after infection or 8 if it is ~80 weeks after infection
      spline_degree = 2L,
      spline_knots_type = 'domain', #use domain or frequency (according to where are more data points density)
      spline_r_linear = TRUE,
      spline_l_linear = FALSE,
      prior_spline_monotonicity = 'decreasing'
      # prior_spline_convexity = "convex"
      # prior_spline_maxder_gaussian = array(c(0, 0.01))
      # prior_spline_maxder_gaussian = rbind(c(0,0,0,0,-1), c(Inf,Inf,Inf,Inf,0.0001))
    )
  )
)

mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)

df_pred3 <- data.frame(mid_point = seq(0, 78, by = 0.1)) #"seq(0, 78" 78 need to change based on maximum mid_point
                                                         #for each group of variant and outcome

dat_pred3 <- MRData()
dat_pred3$load_df(
  data = df_pred3, 
  col_covs=list('mid_point')
)

df_pred3$pred5 <- mod1$predict(data = dat_pred3)

#3.2.3 â Uncertainty from fixed effects and between-study heterogeneity

n_samples3 <- 1000L
set.seed(1)
samples3 <- mod1$sample_soln(sample_size = n_samples3)

draws3 <- mod1$create_draws(
  data = dat_pred3,
  beta_samples = samples3[[1]],
  gamma_samples = samples3[[2]], 
  random_study = TRUE )

df_pred3$pred5 <- mod1$predict(dat_pred3)
df_pred3$pred_lo <- apply(draws3, 1, function(x) quantile(x, 0.025)) #, na.rm = T
df_pred3$pred_hi <- apply(draws3, 1, function(x) quantile(x, 0.975)) #, na.rm = T

#convert pfin logit to pinf pinf = (plogis(pred5)
df_pred3 <- mutate(df_pred3, pinf = (plogis(pred5))) 
df_pred3 <- mutate(df_pred3, pinf_lower = (plogis(pred_lo))) 
df_pred3 <- mutate(df_pred3, pinf_upper = (plogis(pred_hi))) 


#change to inf, symp or sev according each curve
#plot to see knots
inf <- mutate(inf, insesqua = (1/sqrt(se))/10)

with(inf, plot(mid_point, efficacy_mean, xlim = c(0, 90), ylim = c(0, 1), cex = insesqua)) 
with(df_pred3, lines(mid_point, pinf))

#visualize knot locations
for (k in mod1$cov_models[[2]]$spline_knots) abline(v = k, col = "gray")

groups1 <- unique(inf$study_id2)

for (grp in groups1) {
  df_tmp <- filter(inf, study_id2 == grp)
  with(arrange(df_tmp, mid_point), lines(mid_point, efficacy_mean, lty = 2, col = "gray"))
}

# where the knot locations are on the spline  
get_knots <- function(model, cov_model_name) {
  model$cov_models[[which(model$cov_model_names == cov_model_name)]]$spline$knots
}

get_knots(model = mod1, cov_model_name = "mid_point")

#saving curves estimates
#Ancestral, Alpha, Delta
#past_inf <- df_pred3
#write.csv(past_inf, "FILE_PATH/pinf_curves/past_inf.csv", row.names = F)

#past_symp <- df_pred3
#write.csv(past_symp, "FILE_PATH/pinf_curves/past_symp.csv", row.names = F)

#past_sev <- df_pred3
#write.csv(past_sev, "FILE_PATH/pinf_curves/past_sev.csv", row.names = F)

#Omicron
#past_inf_o <- df_pred3
#write.csv(past_inf_o, "FILE_PATH/pinf_curves/past_inf_o.csv", row.names = F)

#past_inf_oba2 <- df_pred3
#write.csv(past_inf_oba2, "FILE_PATH/pinf_curves/past_inf_oba2.csv", row.names = F)

#past_symp_o <- df_pred3
#write.csv(past_symp_o, "FILE_PATH/pinf_curves/past_symp_o.csv", row.names = F)

#past_sev_o <- df_pred3
#write.csv(past_sev_o, "FILE_PATH/pinf_curves/past_sev_o.csv", row.names = F)


#######plot waning prev_inf (Figure 3)
library(ggplot2)
library(formattable)
f1 <- "Times"

previnfplot <- ggplot() + 
  geom_line(data = df_pred3, mapping = aes(x = mid_point, y = pinf), size = 1) +
  geom_point(data = inf, mapping = aes(x = mid_point, y = efficacy_mean, color = study_id2, size = insesqua), #)+ 
             show.legend = TRUE) +
  geom_ribbon(data=df_pred3, aes(x=mid_point, ymin=pinf_lower, ymax=pinf_upper), alpha=0.4, fill="pink") +
  ylim(c(0, 1))+
  xlim(c(0, 90)) +
  theme_classic() +
  theme(text=element_text(size=18, family = f1)) +
  theme(legend.position = "bottom") +
  guides(linetype = FALSE) +
  labs(y="Previous infection protection", x = "Week after infection", size = "Inverse of variance") + #colour = "Primary and subsequent variant"
  geom_hline(aes(yintercept = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), linetype = "yellow")) +
  ggtitle("Panel A: Protection against Ancestral, Alpha, and Delta reinfection") +
  guides(color = "none")
print(previnfplot)

#Use different title according each plot
#Panel A: Protection against Ancestral, Alpha, and Delta reinfection
#Panel B: Protection against Omicron BA.1 reinfection
#Panel B1: Protection against Omicron BA.2 reinfection
#Panel C: Protection against Ancestral, Alpha, and Delta symptomatic disease
#Panel D: Protection against Omicron BA.1 symptomatic disease
#Panel E: Protection against Ancestral, Alpha, and Delta severe disease
#Panel F: Protection against Omicron BA.1 severe disease


