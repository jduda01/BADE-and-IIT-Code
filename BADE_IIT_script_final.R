##########Load everything#########
packages <- c("readxl","lme4","lmerTest", "ggplot2","tidyverse",
              "openxlsx","mediation","plotrix","ppcor","moments",
              "sjPlot","polycor","psych","simr", "effectsize")
lapply(packages, library, character.only=TRUE)


#Gender differences
p_bingen = filter(p, gender_ch=="female" | gender_ch=="male")

t.test(pq_tot ~ gender_ch, data = p_bingen)
t.test(bslAve ~ gender_ch, data = p_bingen)
t.test(bdiSumWAve ~ gender_ch, data = p_bingen)

cohens_d(pq_tot ~ gender_ch, data = p_bingen)
cohens_d(bslAve ~ gender_ch, data = p_bingen)
cohens_d(bdiSumWAve ~ gender_ch, data = p_bingen)

#Mean center
p$bdi_ctr = p$bdiSumWAve - mean(p$bdiSumWAve,na.rm=TRUE)
p$pq_ctr = p$pq_tot - mean(p$pq_tot, na.rm=TRUE)
p$nib_ctr = p$nib_ave - mean(p$nib_ave, na.rm=TRUE)
p$pib_ctr = p$pib_ave - mean(p$pib_ave, na.rm=TRUE)

p$fact1_ctr = p$fact1 - mean(p$fact1, na.rm =TRUE)
p$fact2_ctr = p$fact2 - mean(p$fact2, na.rm =TRUE)
p$fact3_ctr = p$fact3 - mean(p$fact3, na.rm =TRUE)
p$fact4_ctr = p$fact4 - mean(p$fact4, na.rm =TRUE)


###############  Main Mixed-Effects Models ################
p_flex_long = dplyr::select(p, subj, pq_ctr, bdi_ctr, pos_update_ave,	neg_update_ave, nib_ctr, pib_ctr) %>%
  gather("flex_type", "flexibility", contains("ave"))
p_flex_long$valence = -0.5*(p_flex_long$flex_type == "neg_update_ave") + 0.5*(p_flex_long$flex_type == "pos_update_ave")

me_main1 = lmer(flexibility ~ pq_ctr*valence +  bdi_ctr*valence + nib_ctr + pib_ctr + (1|subj), data=p_flex_long)
summary(me_main1)

#Partial ETA
eta_squared(me_main1, partial = TRUE, generalized = FALSE,
            ci = 0.90, alternative = "two.sided")

#Post hoc models (standardized betas)
  p$pq_z = scale(p$pq_tot)
  p$bdi_z = scale(p$bdiSumWAve)
  p$nib_z = scale(p$nib_ave)
  p$pib_z = scale(p$pib_ave)
  p$neg_update_z = scale(p$neg_update_ave)
  p$pos_update_z = scale(p$pos_update_ave)

  lm_neg_flex_main1_z = lm(neg_update_z ~ pq_z + bdi_z  + nib_z + pib_z, data = p)
  summary(lm_neg_flex_main1_z)
  
  lm_pos_flex_main1_z = lm(pos_update_z ~ pq_z + bdi_z  + nib_z + pib_z, data = p)
  summary(lm_pos_flex_main1_z)
        

###############Factor Analysis Mixed-Effects Models################
p_flex_fa_long = dplyr::select(p, subj, fact1_ctr, fact2_ctr, fact3_ctr, fact4_ctr,
                               pos_update_ave,	neg_update_ave_JD, nib_ctr, pib_ctr) %>%
  gather("flex_type", "flexibility", contains("ave"))

p_flex_fa_long$valence = -0.5*(p_flex_fa_long$flex_type == "neg_update_ave_JD") + 0.5*(p_flex_fa_long$flex_type == "pos_update_ave")

#Factor 1
me_fact1 = lmer(flexibility ~ fact1_ctr*valence + nib_ctr + pib_ctr + (1|subj), data = p_flex_fa_long)
summary(me_fact1)

  #Partial ETA Sq
  eta_squared(me_fact1, partial = TRUE, generalized = FALSE,
              ci = 0.90, alternative = "two.sided")

#Factor 2
me_fact2 = lmer(flexibility ~ fact2_ctr*valence + nib_ctr + pib_ctr + (1|subj), data = p_flex_fa_long)
summary(me_fact2)

    #Partial ETA Sq
    eta_squared(me_fact2, partial = TRUE, generalized = FALSE,
                ci = 0.90, alternative = "two.sided")

#Factor 3
me_fact3 = lmer(flexibility ~ fact3_ctr*valence + nib_ctr + pib_ctr + (1|subj), data = p_flex_fa_long)
summary(me_fact3)

  #Partial ETA Sq
  eta_squared(me_fact3, partial = TRUE, generalized = FALSE,
              ci = 0.90, alternative = "two.sided")

#Factor 4
me_fact4 = lmer(flexibility ~ fact4_ctr*valence + nib_ctr + pib_ctr + (1|subj), data = p_flex_fa_long)
summary(me_fact4)

  #Partial ETA Sq
  eta_squared(me_fact4, partial = TRUE, generalized = FALSE,
              ci = 0.90, alternative = "two.sided")


######################IIT ANALYSES######################
#Gender differences
p$gender_ch = ifelse(p$gender==1, "Male", "Female")

t.test(best_a ~ gender_ch, data = p)
t.test(best_b ~ gender, data = p)
t.test(BDISumWAve ~ gender_ch, data = p)

cohens_d(best_a ~ gender_ch, data = p)
cohens_d(best_b ~ gender_ch, data = p)
cohens_d(BDISumWAve ~ gender, data = p)
    
#Mean center
bdiMean = mean(p$BDISumWAve, na.rm=TRUE) #checked code
p$bdi_ctr <- p$BDISumWAve - bdiMean

best_a_mean = mean(p$best_a, na.rm=TRUE)#checked code
p$best_a_ctr = p$best_a - best_a_mean

best_b_mean = mean(p$best_b, na.rm=TRUE)#checked code
p$best_b_ctr = p$best_b - best_b_mean

p$nib_ctr = p$NIB - mean(p$NIB, na.rm = TRUE)#checked code
p$pib_ctr = p$PIB - mean(p$PIB, na.rm = TRUE)

#Convert to long form
p_iit_long = dplyr::select(p, id, N_flex, P_flex, bdi_ctr, best_a_ctr, best_b_ctr, nib_ctr, pib_ctr) %>%
  gather("flex_type", "flex", ends_with("_flex"))
p_iit_long$Valence = -0.5*(p_iit_long$flex_type == "N_flex") + 0.5*(p_iit_long$flex_type == "P_flex")
   
#BEST A
me_iit_besta = lmer(flex  ~ best_a_ctr*Valence + bdi_ctr*Valence + nib_ctr + pib_ctr + (1|id), data= p_iit_long)
summary(me_iit_besta) 
    
#Partial ETA Sq
eta_squared(me_iit_besta, partial = TRUE, generalized = FALSE,
            ci = 0.90, alternative = "two.sided")
      
#BEST B
me_iit_bestb = lmer(flex  ~ best_b_ctr*Valence + bdi_ctr*Valence + nib_ctr + pib_ctr + (1|id), data= p_iit_long)
summary(me_iit_bestb) 
    
eta_squared(me_iit_bestb, partial = TRUE, generalized = FALSE,
            ci = 0.90, alternative = "two.sided")
    
    