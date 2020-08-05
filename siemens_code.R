# Code for 'Spawning salmon density influences fruit production of salmonberry (Rubus spectabilis)'

#### Data preparations ####
# Load packages
library(ggplot2)
library(dplyr)
library(lme4)
library(car)
library(MuMIn)
library(AICcmodavg)
library(DHARMa)

# Read in data
siemens_dataset <- read.csv("siemens_dataset", stringsAsFactors = FALSE, strip.white = TRUE)
data <- siemens_dataset

# Prepare data for use 
data_brix <- data %>% 
  filter(fruit_quality == "R") %>% 
  drop_na(brix_value) 

data_weight <- data %>% 
  filter(fruit_quality == "R") %>% 
  drop_na(weight)

data_seed <- data %>% 
  filter(fruit_quality == "R")

data_fruit <- data %>% 
  #select (-c(reach, quadrat))
  dplyr::select(stream, shrub, fruit_count, salmon_density, pink_density, chum_density, pink_density_scaled, chum_density_scaled, salmon_density_scaled, relative_moisture_scaled, relative_soil_moisture, stem_density_scaled, stem_density, canopy_cover_scaled, canopy_cover, distance_upstream_scaled, distance_upstream, distance_from_stream_scaled, distance_from_stream, northness_scaled, northness, eastness_scaled, eastness) %>% 
  distinct()

#### Fruit count global model ####

# Fruit count model workflow notes:
# 1) Tried a Linear Mixed Effects Model (lmer): residual plot and qq plot show assumption violations
# 2) Then tried a General Linear Mixed Effects Model: data was over dispersed and uniformity issue was present
# 3) Then tried a General Linear Mixed Effects Model with a negative binomial distribution: model fits -> Used this model, see below

# Fruit Count Global Model: GLMM-Negative Binomial 
fruit_count_global <- glmer.nb(fruit_count ~ salmon_density_scaled + relative_moisture_scaled + stem_density_scaled + canopy_cover_scaled + distance_upstream_scaled + distance_from_stream_scaled + northness_scaled + eastness_scaled + (1|stream), data = data_fruit, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))

isSingular(fruit_count_global) # FALSE - Model is not on the boundary/singular

summary(fruit_count_global) 

# Check VIFs (variance inflation factors): Must be under 3 (Zuur, Ieno & Elphick 2010)
car::vif(fruit_count_global) # No values over 2

# Check Global Model: See Dharma diagnostics below since GLMM
hist(summary(fruit_count_global)$residuals) 
plot(x = fitted(fruit_count_global), y = resid(fruit_count_global), main = 'Residual Plot - Fruit Count', xlab = 'Fitted Value', ylab = 'Residual', las = 1, pch = 16) 
abline(a = 0, b = 0, lty = 2, col = 'red') 
qqnorm(resid(fruit_count_global), pch = 16, main = 'QQ Plot - Fruit Count', xlab = 'Theoretical Quantiles', ylab = 'Sample Quantiles')
qqline(resid(fruit_count_global)) 

# Analyze all possible combinations of variables to calculate meaningful Relative Variable importance (RVIs)
fruit_count_allsubsets <- dredge(fruit_count_global)
average_fruit_count <- model.avg(fruit_count_allsubsets, fit = TRUE) 
summary(average_fruit_count) # Equal proportions of each fixed effect 

# Check Relative Variable Importance (RVI)
MuMIn::importance(average_fruit_count) 

# GLMM MODEL DIAGNOSTICS: USE DHARMA PACKAGE 
# Simulate residuals
simulation_output_fruit_count <- simulateResiduals(fittedModel = fruit_count_global, n = 166)

# Residual check on simulated residuals 
plot(simulation_output_fruit_count) 
qqnorm(resid(fruit_count_global))
qqline(resid(fruit_count_global))

# Check individual predictors
par(mfrow = c(1,1))
plotResiduals(data_fruit$salmon_density_scaled, simulation_output_fruit_count$scaledResiduals)
plotResiduals(data_fruit$relative_moisture_scaled, simulation_output_fruit_count$scaledResiduals)
plotResiduals(data_fruit$stem_density_scaled, simulation_output_fruit_count$scaledResiduals)
plotResiduals(data_fruit$canopy_cover_scaled, simulation_output_fruit_count$scaledResiduals)
plotResiduals(data_fruit$distance_upstream_scaled, simulation_output_fruit_count$scaledResiduals)
plotResiduals(data_fruit$distance_from_stream_scaled, simulation_output_fruit_count$scaledResiduals)
plotResiduals(data_fruit$northness_scaled, simulation_output_fruit_count$scaledResiduals)
plotResiduals(data_fruit$eastness_scaled, simulation_output_fruit_count$scaledResiduals)

ggplot(data.frame(salmon_density_scaled=data_fruit$salmon_density_scaled,pearson=residuals(fruit_count_global, type="pearson")), aes(x=salmon_density_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(relative_moisture_scaled=data_fruit$relative_moisture_scaled,pearson=residuals(fruit_count_global, type="pearson")), aes(x=relative_moisture_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(stem_density_scaled=data_fruit$stem_density_scaled,pearson=residuals(fruit_count_global, type="pearson")), aes(x=stem_density_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(canopy_cover_scaled=data_fruit$canopy_cover_scaled,pearson=residuals(fruit_count_global, type="pearson")), aes(x=canopy_cover_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(distance_upstream_scaled=data_fruit$distance_upstream_scaled,pearson=residuals(fruit_count_global, type="pearson")), aes(x=distance_upstream_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(distance_from_stream_scaled=data_fruit$distance_from_stream_scaled,pearson=residuals(fruit_count_global, type="pearson")), aes(x=distance_from_stream_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(northness_scaled=data_fruit$northness_scaled,pearson=residuals(fruit_count_global, type="pearson")), aes(x=northness_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(eastness_scaled=data_fruit$eastness_scaled,pearson=residuals(fruit_count_global, type="pearson")), aes(x=eastness_scaled,y=pearson)) + geom_point() + theme_bw()

# Goodness of fit test
testUniformity(simulationOutput = simulation_output_fruit_count) # No uniformity issue

# Check for Outliers
testOutliers(simulationOutput = simulation_output_fruit_count) # Two potential outliers - Not an issue (See 'Simulate residuals' code above and the DHARMa Vignette)

# Check Dispersion
testDispersion(simulation_output_fruit_count) # No over/underdispersion detected


#### Fruit count model with chum and pink as separate variables ####

# Fruit count model workflow notes:
# 1) Tried a Linear Mixed Effects Model (lmer): residual plot and qq plot show assuption violations
# 2) Then tried a General Linear Mixed Effects Model: data was over dispersed and uniformity issue was present
# 3) Then tried a General Linear Mixed Effects Model with a negative binomial distribution: model fits -> Used this model, see below

# Fruit count model with chum and pink separate: GLMM - negative binomial
fruit_count_both <- glmer.nb(fruit_count ~ chum_density_scaled + pink_density_scaled + relative_moisture_scaled + stem_density_scaled + canopy_cover_scaled + distance_upstream_scaled + distance_from_stream_scaled + northness_scaled + eastness_scaled + (1|stream), data = data_fruit, na.action = "na.fail", control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))) 

isSingular(fruit_count_both) # FALSE - Model is not on the boundary/singular
summary(fruit_count_both) 

# Check VIFs (variance inflation factors): 
car::vif(fruit_count_both) # No values over 2

# Check model: Use Dharma package diagnostics below because it is a GLMM
hist(summary(fruit_count_both)$residuals) 
plot(x = fitted(fruit_count_both), y = resid(fruit_count_both), main = 'Residual Plot - Fruit Count', xlab = 'Fitted Value', ylab = 'Residual', las = 1, pch = 16) 
abline(a = 0, b = 0, lty = 2, col = 'red') 
qqnorm(resid(fruit_count_both), pch = 16, main = 'QQ Plot - Fruit Count', xlab = 'Theoretical Quantiles', ylab = 'Sample Quantiles')
qqline(resid(fruit_count_both)) 

# Analyze all possible combinations of variables to calculate meaningful Relative Variable importance (RVIs)
fruit_count_both_allsubsets <- dredge(fruit_count_both)
average_fruit_count_both <- model.avg(fruit_count_both_allsubsets, fit = TRUE) 
summary(average_fruit_count_both) # Equal proportions of each fixed effect 

# Check relative variable importance (RVI)
MuMIn::importance(average_fruit_count_both) 

# GLMM MODEL DIAGNOSTICS: USE DHARMA PACKAGE 
# Simulate residuals
simulation_output_fruit_count_both <- simulateResiduals(fittedModel = fruit_count_both, n = 166) 

# Residual Check on simulated residuals 
plot(simulation_output_fruit_count_both) 
qqnorm(resid(fruit_count_both))
qqline(resid(fruit_count_both))

#Check individual predictors
par(mfrow = c(1,1))
plotResiduals(data_fruit$chum_density_scaled, simulation_output_fruit_count_both$scaledResiduals)
plotResiduals(data_fruit$pink_density_scaled, simulation_output_fruit_count_both$scaledResiduals)
plotResiduals(data_fruit$relative_moisture_scaled, simulation_output_fruit_count_both$scaledResiduals)
plotResiduals(data_fruit$stem_density_scaled, simulation_output_fruit_count_both$scaledResiduals)
plotResiduals(data_fruit$canopy_cover_scaled, simulation_output_fruit_count_both$scaledResiduals)
plotResiduals(data_fruit$distance_upstream_scaled, simulation_output_fruit_count_both$scaledResiduals)
plotResiduals(data_fruit$distance_from_stream_scaled, simulation_output_fruit_count_both$scaledResiduals)
plotResiduals(data_fruit$northness_scaled, simulation_output_fruit_count_both$scaledResiduals)
plotResiduals(data_fruit$eastness_scaled, simulation_output_fruit_count_both$scaledResiduals)

ggplot(data.frame(chum_density_scaled=data_fruit$chum_density_scaled,pearson=residuals(fruit_count_both, type="pearson")), aes(x=chum_density_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(pink_density_scaled=data_fruit$pink_density_scaled,pearson=residuals(fruit_count_both, type="pearson")), aes(x=pink_density_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(relative_moisture_scaled=data_fruit$relative_moisture_scaled,pearson=residuals(fruit_count_both, type="pearson")), aes(x=relative_moisture_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(stem_density_scaled=data_fruit$stem_density_scaled,pearson=residuals(fruit_count_both, type="pearson")), aes(x=stem_density_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(canopy_cover_scaled=data_fruit$canopy_cover_scaled,pearson=residuals(fruit_count_both, type="pearson")), aes(x=canopy_cover_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(distance_upstream_scaled=data_fruit$distance_upstream_scaled,pearson=residuals(fruit_count_both, type="pearson")), aes(x=distance_upstream_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(distance_from_stream_scaled=data_fruit$distance_from_stream_scaled,pearson=residuals(fruit_count_both, type="pearson")), aes(x=distance_from_stream_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(northness_scaled=data_fruit$northness_scaled,pearson=residuals(fruit_count_both, type="pearson")), aes(x=northness_scaled,y=pearson)) + geom_point() + theme_bw()

ggplot(data.frame(eastness_scaled=data_fruit$eastness_scaled,pearson=residuals(fruit_count_both, type="pearson")), aes(x=eastness_scaled,y=pearson)) + geom_point() + theme_bw()

# Goodness of fit test
testUniformity(simulationOutput = simulation_output_fruit_count_both) # No Uniformity Issue

# Test outliers
testOutliers(simulation_output_fruit_count_both) # Two potential outliers - Not an issue (See 'Simulate residuals' code above and the DHARMa Vignette)

# Check Dispersion
testDispersion(simulation_output_fruit_count_both) # No over/underdispersion detected

#### Seed count global model ####

# Seed count global model: Linear Mixed Model (LMM)
seed_count_global <- lmer(total_seed_count ~ salmon_density_scaled + relative_moisture_scaled + stem_density_scaled + canopy_cover_scaled + distance_upstream_scaled + distance_from_stream_scaled + northness_scaled + eastness_scaled +  (1|stream/shrub), data = data_seed, REML = FALSE, na.action = na.fail, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))

summary(seed_count_global) 

# Check VIFs (variance inflation factors): 
car::vif(seed_count_global) # No values over 2

# Global model diagnostics: assumptions appear to be met
hist(summary(seed_count_global)$residuals)  # Residuals normally distributed
plot(x = fitted(seed_count_global), y = resid(seed_count_global), main = 'Residual Plot - Seed Count', xlab = 'Fitted Value', ylab = 'Residual', las = 1, pch = 16) 
abline(a = 0, b = 0, lty = 2, col = 'red') # Ok
qqnorm(resid(seed_count_global), pch = 16, main = 'QQ Plot - Seed Count', xlab = 'Theoretical Quantiles', ylab = 'Sample Quantiles') 
qqline(resid(seed_count_global)) # Ok 

# Analyze all possible combinations of variables to calculate meaningful Relative Variable importance (RVIs)
seed_count_allsubsets <- dredge(seed_count_global)
average_seed_count <- model.avg(seed_count_allsubsets, fit = TRUE) 
summary(average_seed_count) # Equal proportions of each fixed effect 

# Check Relative Variable Importance (RVI)
MuMIn::importance(average_seed_count) 

#### Seed count model with chum and pink as separate variables ####

# Seed count model with chum and pink separate: Liner Mixed Model (LMM)
seed_count_both <- lmer(total_seed_count ~ chum_density_scaled + pink_density_scaled + relative_moisture_scaled + stem_density_scaled + canopy_cover_scaled + distance_upstream_scaled + distance_from_stream_scaled + northness_scaled + eastness_scaled + (1|stream/shrub), data = data_seed, REML = FALSE, na.action = na.fail, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))

summary(seed_count_both) 

# Check VIFs (variance inflation factors): 
car::vif(seed_count_both) # No values over 2

# Model diagnostics: assumptions seem to be met
hist(summary(seed_count_both)$residuals) # Residuals normally distributed
plot(x = fitted(seed_count_both), y = resid(seed_count_both), main = 'Residual Plot - Seed Count', xlab = 'Fitted Value', ylab = 'Residual', las = 1, pch = 16) 
abline(a = 0, b = 0, lty = 2, col = 'red') # Ok
qqnorm(resid(seed_count_both), pch = 16, main = 'QQ Plot - Seed Count', xlab = 'Theoretical Quantiles', ylab = 'Sample Quantiles') 
qqline(resid(seed_count_both)) # Ok

# Analyze all possible combinations of variables to calculate meaningful Relative Variable importance (RVIs)
seed_count_both_allsubsets <- dredge(seed_count_both)
average_seed_count_both <- model.avg(seed_count_both_allsubsets, fit = TRUE) 
summary(average_seed_count_both) # Equal proportions of each fixed effect 

# Check Relative Variable Importance (RVI)
MuMIn::importance(average_seed_count_both) 


#### Fruit size global model ####

# Fruit size global model: Linear Mixed Model (LMM)
weight_global <- lmer(weight ~ salmon_density_scaled + relative_moisture_scaled + stem_density_scaled + canopy_cover_scaled + distance_upstream_scaled + distance_from_stream_scaled + northness_scaled + eastness_scaled + (1|stream/shrub), data = data_weight, REML = FALSE, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))) 

summary(weight_global)

# Check VIFs (variance inflation factors): 
car::vif(weight_global) # No values over 2

# Global model diagnostics: Assumptions seem to be met
hist(summary(weight_global)$residuals) # Residuals normally distributed
plot(x = fitted(weight_global), y = resid(weight_global), main = 'Residual Plot - Weight',
     xlab = 'Fitted Value', ylab = 'Residual', las = 1, pch = 16)  
abline(a = 0, b = 0, lty = 2, col = 'red') # Ok
qqnorm(resid(weight_global), pch = 16, main = 'QQ Plot - Weight', xlab = 'Theoretical Quantiles', ylab = 'Sample Quantiles') 
qqline(resid(weight_global)) # Ok

# Analyze all possible combinations of variables to calculate meaningful Relative Variable importance (RVIs)
weight_allsubsets <- dredge(weight_global)
average_weight <- model.avg(weight_allsubsets, fit = TRUE) 
summary(average_weight) # Equal proportions of each fixed effect 

# Check Relative Variable Importance (RVI)
MuMIn::importance(average_weight) 


#### Fruit size model with chum and pink as separate variables #### 

# Fruit size model with chum and pink separate: Liner Mixed Model (LMM)
weight_both <- lmer(weight ~ chum_density_scaled + pink_density_scaled + relative_moisture_scaled + stem_density_scaled + canopy_cover_scaled + distance_upstream_scaled + distance_from_stream_scaled + northness_scaled + eastness_scaled + (1|stream/shrub), data = data_weight, REML = FALSE, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6))) 

summary(weight_both)

# Check VIFs (variance inflation factors): 
car::vif(weight_both) # No values over 2

# Global model diagnostics: Assumptions seem to be met
hist(summary(weight_both)$residuals) # Residuals normally distributed
plot(x = fitted(weight_both), y = resid(weight_both), main = 'Residual Plot - Weight',
     xlab = 'Fitted Value', ylab = 'Residual', las = 1, pch = 16) 
abline(a = 0, b = 0, lty = 2, col = 'red') # Ok
qqnorm(resid(weight_both), pch = 16, main = 'QQ Plot - Weight', xlab = 'Theoretical Quantiles', ylab = 'Sample Quantiles') 
qqline(resid(weight_both)) # Ok

# Analyze all possible combinations of variables to calculate meaningful Relative Variable importance (RVIs)
weight_both_allsubsets <- dredge(weight_both)
average_weight_both <- model.avg(weight_both_allsubsets, fit = TRUE) 
summary(average_weight_both) # Equal proportions of each fixed effect 

# Check Relative Variable Importance (RVI)
MuMIn::importance(average_weight_both) 

#### Brix value (estimated sugar content) global model ####

# Brix value global model: Linear Mixed Model (LMM)
brix_value_global <- lmer(brix_value ~ salmon_density_scaled + relative_moisture_scaled + stem_density_scaled + canopy_cover_scaled + distance_upstream_scaled + distance_from_stream_scaled + northness_scaled + eastness_scaled + (1|stream/shrub), data = data_brix, REML = FALSE, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))
summary(brix_value_global)

# Check VIFs (variance inflation factors): 
car::vif(brix_value_global) # No values over 2

# Global Model Diagnostics: Assumptions seem to be met
hist(summary(brix_value_global)$residuals)  # Residuals normally distributed
plot(x = fitted(brix_value_global), y = resid(brix_value_global), main = 'Residual Plot - Brix Value', xlab = 'Fitted Value', ylab = 'Residual', las = 1, pch = 16) 
abline(a = 0, b = 0, lty = 2, col = 'red') # Ok
qqnorm(resid(brix_value_global), pch = 16, main = 'QQ Plot - Brix Value', xlab = 'Theoretical Quantiles', ylab = 'Sample Quantiles') 
qqline(resid(brix_value_global)) # Ok

# Analyze all possible combinations of variables to calculate meaningful Relative Variable importance (RVIs)
brix_value_allsubsets <- dredge(brix_value_global)
average_brix_value <- model.avg(brix_value_allsubsets, fit = TRUE) 
summary(average_brix_value) # Equal proportions of each fixed effect 

# Check Relative Variable Importance (RVI)
MuMIn::importance (average_brix_value) 


#### Brix value (estimated sugar content) model with chum and pink as separate variables ####

# Brix value model with chum and pink separate: Liner Mixed Model (LMM)
brix_value_both <- lmer(brix_value ~ chum_density_scaled + pink_density_scaled + relative_moisture_scaled + stem_density_scaled + canopy_cover_scaled + distance_upstream_scaled + distance_from_stream_scaled + northness_scaled + eastness_scaled + (1|stream/shrub), data = data_brix, REML = FALSE, na.action = "na.fail", control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e6)))

summary(brix_value_both)

# Check VIFs (variance inflation factors): 
car::vif(brix_value_both) # No values over 2

# Global model diagnostics: Assumptions seem to be met
hist(summary(brix_value_both)$residuals) # Residuals normally distributed
plot(x = fitted(brix_value_both), y = resid(brix_value_both), main = 'Residual Plot - Brix Value',
     xlab = 'Fitted Value', ylab = 'Residual', las = 1, pch = 16) 
abline(a = 0, b = 0, lty = 2, col = 'red') # Ok
qqnorm(resid(brix_value_both), pch = 16, main = 'QQ Plot - Brix Value', xlab = 'Theoretical Quantiles', ylab = 'Sample Quantiles') 
qqline(resid(brix_value_both)) # Ok

# Analyze all possible combinations of variables to calculate meaningful Relative Variable importance (RVIs)
brix_value_both_allsubsets <- dredge(brix_value_both)
average_brix_value_both <- model.avg(brix_value_both_allsubsets, fit = TRUE) 
summary(average_brix_value_both)  # Equal proportions of each fixed effect 

# Check Relative Variable Importance (RVI)
MuMIn::importance (average_brix_value_both)


#### Fruit count predictions - Global model ####

# Fruit count global model
# Average stream:
avg_stream_salmon <- expand.grid(salmon_density_scaled = mean(data_fruit$salmon_density_scaled),
                                     relative_moisture_scaled = mean(data_fruit$relative_moisture_scaled),
                                     stem_density_scaled = mean(data_fruit$stem_density_scaled),
                                     canopy_cover_scaled = mean(data_fruit$canopy_cover_scaled),
                                     distance_upstream_scaled = mean(data_fruit$distance_upstream_scaled),
                                     distance_from_stream_scaled = mean(data_fruit$distance_from_stream_scaled), 
                                     northness_scaled = mean(data_fruit$northness_scaled),
                                     eastness_scaled = mean(data_fruit$eastness_scaled))

pred_avg_stream_salmon <- predict(average_fruit_count, avg_stream_salmon, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_avg_stream_salmon$fit # 15.9 fruits per bush for an average stream
pred_avg_stream_salmon$se.fit*1.96 # +/- 4.5 fruits

mean(data_fruit$fruit_count) # mean is ~ 20 so this makes sense.

# More salmon:
more_salmon <- expand.grid(salmon_density_scaled = mean(data_fruit$salmon_density_scaled) + sd(data_fruit$salmon_density_scaled),
                               relative_moisture_scaled = mean(data_fruit$relative_moisture_scaled),
                               stem_density_scaled = mean(data_fruit$stem_density_scaled),
                               canopy_cover_scaled = mean(data_fruit$canopy_cover_scaled),
                               distance_upstream_scaled = mean(data_fruit$distance_upstream_scaled),
                               distance_from_stream_scaled = mean(data_fruit$distance_from_stream_scaled), 
                               northness_scaled = mean(data_fruit$northness_scaled),
                               eastness_scaled = mean(data_fruit$eastness_scaled))

pred_more_salmon <- predict(average_fruit_count, more_salmon, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_more_salmon

pred_more_salmon$fit # 19.6 fruits per bush for a stream with 1 sd more salmon
pred_more_salmon$se.fit*1.96 # +/- 9.6 fruits

sd(data_fruit$salmon_density) # 1 sd salmon = 8.8 kg/m stream 

#### Fruit count predictions - Chum and pink as separate variables ####

# Fruit count model with chum and pink as separate variables:
# Average stream:
avg_stream <- expand.grid(chum_density_scaled = mean(data_fruit$chum_density_scaled),
                              pink_density_scaled = mean(data_fruit$pink_density_scaled),
                              relative_moisture_scaled = mean(data_fruit$relative_moisture_scaled),
                              stem_density_scaled = mean(data_fruit$stem_density_scaled),
                              canopy_cover_scaled = mean(data_fruit$canopy_cover_scaled),
                              distance_upstream_scaled = mean(data_fruit$distance_upstream_scaled),
                              distance_from_stream_scaled = mean(data_fruit$distance_from_stream_scaled), 
                              northness_scaled = mean(data_fruit$northness_scaled),
                              eastness_scaled = mean(data_fruit$eastness_scaled))

pred_avg_stream <- predict(average_fruit_count_both, avg_stream, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_avg_stream$fit # 16.1 fruits per bush for an average stream
pred_avg_stream$se.fit*1.96 # +/- 4.4 fruits

mean(data_fruit$fruit_count)# mean is ~20 so this makes sense.

# More chum:
more_chum <- expand.grid(chum_density_scaled = mean(data_fruit$chum_density_scaled) + sd(data_fruit$chum_density_scaled),
                             pink_density_scaled = mean(data_fruit$pink_density_scaled),
                             relative_moisture_scaled = mean(data_fruit$relative_moisture_scaled),
                             stem_density_scaled = mean(data_fruit$stem_density_scaled),
                             canopy_cover_scaled = mean(data_fruit$canopy_cover_scaled),
                             distance_upstream_scaled = mean(data_fruit$distance_upstream_scaled),
                             distance_from_stream_scaled = mean(data_fruit$distance_from_stream_scaled), 
                             northness_scaled = mean(data_fruit$northness_scaled),
                             eastness_scaled = mean(data_fruit$eastness_scaled))

pred_more_chum <- predict(average_fruit_count_both, more_chum, type = "response", se.fit = TRUE, re.form = NA, full = T)

pred_more_chum

pred_more_chum$fit # 21.1 fruits per bush for a stream with 1 sd more chum
pred_more_chum$se.fit*1.96 # +/- 10.8 fruits

sd(data_fruit$chum_density) # 1 sd chum = 6.8 kg/m stream

