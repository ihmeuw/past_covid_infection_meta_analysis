###Baseline Previous infection immunity##

library(tidyverse)
library(reshape2)
library(data.table)

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

ve2 <- mutate(ve2, location_name = if_else(location_name == "Belgium, Brazil, Colombia, The Philippines and South Africa",
                                           "Multiple countries", ve2$location_name))

#total number of studies
length(unique(ve2$study_id))
length(unique(ve2$author))
length(unique(ve2$location_name))
table(ve2$location_name)

#Table 1 main text
#save a dataset only with sub-lineages values to include in the text
#sublinegaes <- subset(ve2, prim_variant_name == "Omicron BA.1" | prim_variant_name == "Omicron BA.2" &
#                      variant_name == "Omicron BA.1" | variant_name == "Omicron BA.2" | variant_name == "Omicron BA.4/BA.5")
#write.csv(sublinegaes, "FILE_PATH/pinf_sublinegaes.csv", row.names = F)

#keep only baseline studies
ve2 <- subset(ve2, pinf == 0)

#new variant_name variable (to differentiate studies with Omicron sub-lineages)
ve2$variants = paste(ve2$prim_variant_name, ve2$variant_name,  sep = " and ")

#number of studies and countries (only baseline previous infection studies)
length(unique(ve2$study_id))
length(unique(ve2$location_name))
table(ve2$location_name)
table(ve2$study_id)

#run the model using logit space
ve2$efficacy_upper <- as.numeric(as.character(ve2$efficacy_upper))
ve2$efficacy_lower <- as.numeric(as.character(ve2$efficacy_lower))
ve2 <- mutate(ve2, se = (efficacy_upper - efficacy_lower)/3.92)

#transform to logit
library(crosswalk002, lib.loc = "/ihme/code/mscm/Rv4/packages/")
logit <- delta_transform(mean = ve2$efficacy_mean, sd = ve2$se, transformation = "linear_to_logit")

names(logit) <- c("mean_logit", "sd_logit")
vacc_logit <- cbind(ve2, logit)
ve2 <- vacc_logit

#Subset - infection, symptomatic, severe
#infec <- subset(ve2, sev_severity == "infection")
#symptomatic <- subset(ve2, sev_severity == "symptomatic")
#severe <- subset(ve2, sev_severity == "severe")

##Forest plots (Figure S1)
###subset by variant
#inf <- subset(infec, variant_name == "Ancestral")
#inf <- subset(infec, variant_name == "Alpha")
#inf <- subset(infec, variant_name == "Beta")
#inf <- subset(infec, variant_name == "Delta")
#inf <- subset(infec, variant_name == "Omicron BA.1")

#symp <- subset(symptomatic, variant_name == "Ancestral")
#symp <- subset(symptomatic, variant_name == "Alpha")
#symp <- subset(symptomatic, variant_name == "Beta")
#symp <- subset(symptomatic, variant_name == "Delta")
#symp <- subset(symptomatic, variant_name == "Omicron BA.1")

#sev <- subset(severe, variant_name == "Ancestral")
#sev <- subset(severe, variant_name == "Alpha")
#sev <- subset(severe, variant_name == "Beta")
#sev <- subset(severe, variant_name == "Delta")
#sev <- subset(severe, variant_name == "Omicron BA.1")


#Sensitivity analysis (Figure S2)
#subset by quality_rating Fair, Good - Sensitivity analysis

#inf <- subset(inf, quality_rating == "Fair" | quality_rating == "Good")
#inf <- subset(inf, quality_rating == "Good")

#symp <- subset(symp, quality_rating == "Fair" | quality_rating == "Good")
#symp <- subset(symp, quality_rating == "Good")

#sev <- subset(sev, quality_rating == "Fair" | quality_rating == "Good")
#sev <- subset(sev, quality_rating == "Good")


#Sensitivity analysis (Figure S3)
#subset by adjusted/unadjusted results - Sensitivity analysis

#inf <- subset(inf, adjustment == 0)
#inf <- subset(inf, adjustment == 1)
#inf <- subset(inf, adjustment == 2)
#inf <- subset(inf, adjustment == 3)

#symp <- subset(symp, adjustment == 0)
#symp <- subset(symp, adjustment == 1)
#symp <- subset(symp, adjustment == 2)
#symp <- subset(symp, adjustment == 3)

#sev <- subset(v, adjustment == 0)
#sev <- subset(sev, adjustment == 1)
#sev <- subset(sev, adjustment == 2)
#sev <- subset(sev, adjustment == 3)


#funnel plot (Figure S4)
#publication bias Funnel plot

library(metafor)
library(meta)

#refline is the meta-analysis result

#infection
#funnel(x = inf$efficacy_mean, sei = inf$se, refline=0.849, xlim = c(0,1.5), ylim = c(0,0.3)) #ancestral
#funnel(x = inf$efficacy_mean, sei = inf$se, refline=0.90, xlim = c(0,1.5), ylim = c(0,0.3)) #alpha 
#funnel(x = inf$efficacy_mean, sei = inf$se, refline=0.857, xlim = c(0,1.5), ylim = c(0,0.3)) #beta
#funnel(x = inf$efficacy_mean, sei = inf$se, refline=0.82, xlim = c(0,1.5), ylim = c(0,0.3)) #delta
#funnel(x = inf$efficacy_mean, sei = inf$se, refline=0.453, xlim = c(0,1.5), ylim = c(0,0.3)) #omicron
#regtest(x = inf$efficacy_mean, sei = inf$se, model="rma", predictor="sei")

#simptomatic
#funnel(x = symp$efficacy_mean, sei = symp$se, refline=0.821, xlim = c(0,1.5), ylim = c(0,0.3)) #ancestral
#funnel(x = symp$efficacy_mean, sei = symp$se, refline=0.872, xlim = c(0,1.5), ylim = c(0,0.3)) #alpha
#funnel(x = symp$efficacy_mean, sei = symp$se, refline=0.854, xlim = c(0,1.5), ylim = c(0,0.3)) #beta
#funnel(x = symp$efficacy_mean, sei = symp$se, refline=0.85, xlim = c(0,1.5), ylim = c(0,0.3)) #delta
#funnel(x = symp$efficacy_mean, sei = symp$se, refline=0.44, xlim = c(0,1.5), ylim = c(0,0.3)) #omicron
#regtest(x = symp$efficacy_mean, sei = symp$se, model="rma", predictor="sei")

#severe
#funnel(x = sev$efficacy_mean, sei = sev$se, refline=0.781, xlim = c(0,1.5), ylim = c(0,0.3)) #ancestral
#funnel(x = sev$efficacy_mean, sei = sev$se, refline=0.796, xlim = c(0,1.5), ylim = c(0,0.3)) #alpha
#funnel(x = sev$efficacy_mean, sei = sev$se, refline=0.88, xlim = c(0,1.5), ylim = c(0,0.3)) #beta
#funnel(x = sev$efficacy_mean, sei = sev$se, refline=0.972, xlim = c(0,1.5), ylim = c(0,0.3)) #delta
#funnel(x = sev$efficacy_mean, sei = sev$se, refline=0.819, xlim = c(0,1.5), ylim = c(0,0.3)) #omicron
#regtest(x = sev$efficacy_mean, sei = sev$se, model="rma", predictor="sei")


#Run MR-BRT
library(dplyr)
library(mrbrt002, lib.loc = "/ihme/code/mscm/Rv4/packages/")
library(data.table)

dt_combined <- data.table()

#inf
#symp
#sev

#1. Fitting a standard mixed effects model
ve6 <- MRData()
ve6$load_df(
  data = inf,  
  col_obs = "mean_logit",    
  col_obs_se = "sd_logit",   
  #col_covs = as.list(new_covs),
  col_study_id = 'study_id' )

mod1 <- MRBRT(
  data = ve6,
  cov_models = list(
    LinearCovModel("intercept", use_re = TRUE)
    #LinearCovModel(as.list(new_covs) ) 
  ))

mod1$fit_model(inner_print_level = 5L, inner_max_iter = 1000L)

pred_data <- as.data.table(expand.grid("intercept"=c(1)))

dat_pred1  <- MRData()
dat_pred1 $load_df(
  data = pred_data,
  col_covs = list('intercept')
)

mod1$beta_soln
mod1$gamma_soln

#3.2.2 â Uncertainty from fixed effects only (using fit-refit)
n_samples2 <- 1000L 
set.seed(1)

samples2_fitrefit <- mod1$sample_soln(sample_size = n_samples2)

draws2_fitrefit <- mod1$create_draws(
  data = dat_pred1,
  beta_samples = samples2_fitrefit[[1]],
  gamma_samples = samples2_fitrefit[[2]],
  random_study = FALSE )

pred_data$pred2 <- mod1$predict(data = dat_pred1)
pred_data$pred2_lo <- apply(draws2_fitrefit, 1, function(x) quantile(x, 0.025)) #, na.rm = T
pred_data$pred2_hi <- apply(draws2_fitrefit, 1, function(x) quantile(x, 0.975)) #, na.rm = T

#summary(pred_data)
#summary
exp(pred_data$pred2)
exp(pred_data$pred2_lo )
exp(pred_data$pred2_hi)

#3.2.3 â Uncertainty from fixed effects and between-study heterogeneity

n_samples3 <- 1000L
set.seed(1)
samples3 <- mod1$sample_soln(sample_size = n_samples3)

draws3 <- mod1$create_draws(
  data = dat_pred1,
  beta_samples = samples3[[1]],
  gamma_samples = samples3[[2]], #matrix(rep(mod1$gamma_soln, n_samples3), ncol = 1) #samples3[[2]],
  random_study = TRUE )

pred_data$pred3 <- mod1$predict(dat_pred1)
pred_data$pred3_lo <- apply(draws3, 1, function(x) quantile(x, 0.025)) #, na.rm = T
pred_data$pred3_hi <- apply(draws3, 1, function(x) quantile(x, 0.975)) #, na.rm = T

#summary(pred_data)
#summary
exp(pred_data$pred3)
exp(pred_data$pred3_lo )
exp(pred_data$pred3_hi)

#spreadsheet 
dtfinal <- copy(pred_data)

dtfinal[, pred2_exp:= exp(pred2)]
dtfinal[, pred2_lo_exp:= exp(pred2_lo)]
dtfinal[, pred2_hi_exp:= exp(pred2_hi)]
dtfinal[, pred3_exp:= exp(pred3)]
dtfinal[, pred3_lo_exp:= exp(pred3_lo)]
dtfinal[, pred3_hi_exp:= exp(pred3_hi)]

dt_combined <- rbind(dt_combined, dtfinal)

#convert inf protection logit to past inf protection
ve3 <- mutate(dt_combined, pastinf = (plogis(pred3))) 
ve3 <- mutate(ve3, inf_lower = (plogis(pred3_lo))) 
ve3 <- mutate(ve3, inf_upper = (plogis(pred3_hi))) 

#transform to percent style
ve3 <- mutate(ve3, pastinf = pastinf*100)
ve3 <- mutate(ve3, inf_lower = inf_lower*100)
ve3 <- mutate(ve3, inf_upper = inf_upper*100)
#inf <- mutate(inf, efficacy_mean = efficacy_mean*100)
#inf <- mutate(inf, efficacy_lower = efficacy_lower*100)
#inf <- mutate(inf, efficacy_upper = efficacy_upper*100)
#symp <- mutate(symp, efficacy_mean = efficacy_mean*100)
#symp <- mutate(symp, efficacy_lower = efficacy_lower*100)
#symp <- mutate(symp, efficacy_upper = efficacy_upper*100)
#sev <- mutate(sev, efficacy_mean = efficacy_mean*100)
#sev <- mutate(sev, efficacy_lower = efficacy_lower*100)
#sev <- mutate(sev, efficacy_upper = efficacy_upper*100)

#Forest plot previous infection_preventing infection_symptomatic_severe disease

#inf
#symp
#sev

gg1 <- fp_dat <- data.frame(
  mean = c(NA, NA, inf$efficacy_mean, NA, ve3$pastinf),
  lower = c(NA, NA, inf$efficacy_lower, NA, ve3$inf_lower),
  upper = c(NA, NA, inf$efficacy_upper, NA, ve3$inf_upper)
)

fp_text <- cbind(
  c("", "Author", as.character(inf$author), NA, "Mean estimate"),
  c("", "Country", as.character(inf$location_name), NA, NA),
  #c("", "Variant", as.character(inf$variants), NA, NA),
  c("", "Protection (95% UI)", paste0(
    format(round(inf$efficacy_mean, digits = 1), nsmall = 1), " (",
    format(round(inf$efficacy_lower, digits = 1), nsmall = 1), "-",
    format(round(inf$efficacy_upper, digits = 1), nsmall = 1), ")"
  ), NA, paste0(
    format(round(ve3[1, "pastinf"], digits = 1), nsmall = 1), " (",
    format(round(ve3[1, "inf_lower"], digits = 1), nsmall = 1), "-",
    format(round(ve3[1, "inf_upper"], digits = 1), nsmall = 1), ")"
  )))

#library(forestplot, lib.loc = "/home/j/temp/reed/prog/R_libraries/")
library(forestplot)

forestplot(fp_text, 
           fp_dat,# new_page = TRUE,
           is.summary=c(rep(TRUE,2), rep(FALSE, nrow(inf)),TRUE),
           #clip=c(0, 100),
           xticks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
           #xlog=FALSE, 
           #txt_gp = fpTxtGp(cex=1.2),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))

print(gg1)

