## Power calculation and simulation for 475 ASE locations
library(tidyverse)
library(lubridate)
library(readxl)
library(here)
library(MASS)
library(fastDummies)
library(fmesher)
library(inlabru)
library(Matrix)
library(INLA)

# Fit a Negative Binomial model to generate simulated data for 475 locations
pre_data <- read_excel("pre_pandemic_data_by 300 ASE location.xlsx")
nb_model <- glm.nb(Collisions ~ Q1 + Q2 + Q3 + T + Pandemic, data = pre_data)
summary(nb_model)

#Create dataframe with Time and Location variables
set.seed(123)  # For reproducibility
n_locations <- 475
n_rows_per_location <- 24
total_rows <- n_locations * n_rows_per_location

df <- data.frame(
  Location = rep(1:n_locations, each = n_rows_per_location),
  Time = rep(-19:4, times = n_locations)
)

# Function to generate dummy variables Q1 to Q4
generate_quarters <- function(time_series) {
  start_q1 <- sample(-19:-16, 1)  
  q1 <- ifelse(((time_series - start_q1) %% 4 == 0) & (time_series >= start_q1), 1, 0)
  q2 <- ifelse(((time_series - start_q1 - 1) %% 4 == 0) & (time_series >= start_q1), 1, 0)
  q3 <- ifelse(((time_series - start_q1 - 2) %% 4 == 0) & (time_series >= start_q1), 1, 0)
  q4 <- ifelse(((time_series - start_q1 - 3) %% 4 == 0) & (time_series >= start_q1), 1, 0)
  return(data.frame(Q1 = q1, Q2 = q2, Q3 = q3, Q4 = q4))
}

# Apply the function to each location
quarters <- df %>%
  group_by(Location) %>%
  do(generate_quarters(.$Time))

# Bind the generated quarter variables back to the original data frame
df.1 <- cbind(df, quarters)
df.1 <- df.1[ , !duplicated(colnames(df.1))]


# Function to generate pandemic variable for a single location
generate_pandemic <- function(time_series) {
  start_pandemic <- sample(-10:0, 1)  # Random start time between -10 and 0 (assuming that for most recent ASE, the pandemic occured 3.5 years ago)
  pandemic <- ifelse(time_series >= start_pandemic, 1, 0)
  return(pandemic)
}

# Apply the function to each location
df.2 <- df.1 %>%
  group_by(Location) %>%
  mutate(pandemic = generate_pandemic(Time)) %>%
  ungroup()

# Define six cities
cities <- c("Toronto", "Montreal", "Guelph", "Saskatoon", "Calgary", "Surrey")

# Assign locations : Total 475
num_locations <- c(260, 53, 27, 47, 53, 35)  

set.seed(123)  

location_data <- data.frame(
  location_id = 1:sum(num_locations),
  city = rep(cities, times = num_locations)
)

city_effects <- runif(6, min = 0, max = 0.1)

# Assign city effects to each location
location_data$city_effect <- city_effects[match(location_data$city, cities)]

random_values <- data.frame(ID = unique(df.2$Location), 
                            random_value = runif(length(unique(df.2$Location))))
df.2$ID <- df.2$Location

df.2 <- df.2 %>% left_join(random_values, by = "ID")

# Simulate new pre-period data using the model coefficients
Intercept <- as.numeric(coef(nb_model)["(Intercept)"])
X1 <- as.numeric(coef(nb_model)[2])
X2 <- as.numeric(coef(nb_model)[3])
X3 <- as.numeric(coef(nb_model)[4])
X4 <- as.numeric(coef(nb_model)[5])
X5 <- as.numeric(coef(nb_model)[6])

## Simulate new pre-data for 475 sites
# N= 475 sites x 24 pre-pandemic Qs = 11400 rows.
# time  1 to 4 will serve for the ASE-PERIOD
set.seed(123)  # For reproducibility
# df.2$Collisions <- rpois(nrow(df.2), exp(Intercept +df.2$random_value + X1*df.2$Q1 + 
#                                    X2*df.2$Q2 + 
#                                    X3*df.2$Q3 + 
#                                    X4*df.2$Time + 
#                                    X5*df.2$pandemic))
#BRICE JULY 2025
# HERE IS THE CODE TO ESTIMATE THETA WITH THE ORIGINAL TORONTO DATA
# Fit a Negative Binomial model
nb_model <- glm.nb(Collisions ~ Q1 + Q2 + Q3 + T + Pandemic, data = pre_data)

# View the estimated theta (dispersion parameter)
theta_est <- nb_model$theta
print(theta_est)
# Dispersion parameter (theta)  
theta <- theta_est  #

# Calculate the mean as before
mu <- exp(Intercept + df.2$random_value + 
            X1 * df.2$Q1 + 
            X2 * df.2$Q2 + 
            X3 * df.2$Q3 + 
            X4 * df.2$Time + 
            X5 * df.2$pandemic)

# Generate Negative Binomial data
df.2$Collisions <- rnbinom(n = nrow(df.2), size = theta, mu = mu)

write.csv(df.2,"simulated_data_2025Aug15.csv", row.names = FALSE)

df.2$Secular <- df.2$Time
df.2 <- dummy_cols(df.2,  select_columns = "ID")
str(df.2)
str(df.2, list.len=ncol(df.2))

# Calculate predicted counts

Sim_pre.2 <- subset(df.2, Time <1)

comp.a <- Collisions ~ Time + Q1 + Q2 + Q3 + pandemic + ID_1 + ID_2 + ID_3 + ID_4 + 
  ID_5 + ID_6 + ID_7 + ID_8 + ID_9 + ID_10 +  ID_11 + ID_12 + ID_13 + ID_14 + 
  ID_15 + ID_16 + ID_17 + ID_18 + ID_19 + ID_20 + ID_21 + ID_22 + ID_23 + 
  ID_24 + ID_25 + ID_26 + ID_27 + ID_28 + ID_29 + ID_30 + ID_31 + ID_32 + 
  ID_33 + ID_34 + ID_35 + ID_36 + ID_37 + ID_38 + ID_39 + ID_40 + ID_41 + 
  ID_42 + ID_43 + ID_44 + ID_45 + ID_46 + ID_47 + ID_48 + ID_49 + ID_50 + 
  ID_51 + ID_52 + ID_53 + ID_54 + ID_55 + ID_56 + ID_57 + ID_58 + ID_59 + ID_60 + 
  ID_61 + ID_62 + ID_63 + ID_64 + ID_65 + ID_66 + ID_67 + ID_68 + ID_69 + 
  ID_70 + ID_71 + ID_72 + ID_73 + ID_74 + ID_75 + ID_76 + ID_77 + ID_78 + 
  ID_79 + ID_80 + ID_81 + ID_82 + ID_83 + ID_84 + ID_85 + ID_86 + ID_87 + 
  ID_88 + ID_89 + ID_90 + ID_91 + ID_92 + ID_93 + ID_94 + ID_95 + ID_96 + 
  ID_97 + ID_98 + ID_99 + ID_100 + ID_101 + ID_102 + ID_103 + ID_104 + ID_105 + 
  ID_106 + ID_107 + ID_108 + ID_109 + ID_110 + ID_111 + ID_112 + ID_113 + 
  ID_114 + ID_115 + ID_116 + ID_117 + ID_118 + ID_119 + ID_120 + ID_121 + 
  ID_122 + ID_123 + ID_124 + ID_125 + ID_126 + ID_127 + ID_128 + ID_129 + 
  ID_130 + ID_131 + ID_132 + ID_133 + ID_134 + ID_135 + ID_136 + ID_137 + 
  ID_138 + ID_139 + ID_140 + ID_141 + ID_142 + ID_143 + ID_144 + ID_145 + 
  ID_146 + ID_147 + ID_148 + ID_149 + ID_150 + ID_151 + ID_152 + ID_153 + 
  ID_154 + ID_155 + ID_156 + ID_157 + ID_158 + ID_159 + ID_160 + ID_161 + 
  ID_162 + ID_163 + ID_164 + ID_165 + ID_166 + ID_167 + ID_168 + 
  ID_169 + ID_170 + ID_171 + ID_172 + ID_173 + ID_174 + ID_175 + 
  ID_176 + ID_177 + ID_178 + ID_179 + ID_180 + ID_181 + ID_182 + 
  ID_183 + ID_184 + ID_185 + ID_186 + ID_187 + ID_188 + ID_189 + 
  ID_190 + ID_191 + ID_192 + ID_193 + ID_194 + ID_195 + ID_196 + 
  ID_197 + ID_198 + ID_199 + ID_200 + ID_201 + ID_202 + ID_203 + 
  ID_204 + ID_205 + ID_206 + ID_207 + ID_208 + ID_209 + ID_210 + 
  ID_211 + ID_212 + ID_213 + ID_214 + ID_215 + ID_216 + ID_217 + 
  ID_218 + ID_219 + ID_220 + ID_221 + ID_222 + ID_223 + ID_224 + 
  ID_225 + ID_226 + ID_227 + ID_228 + ID_229 + ID_230 + ID_231 + 
  ID_232 + ID_233 + ID_234 + ID_235 + ID_236 + ID_237 + ID_238 + 
  ID_239 + ID_240 + ID_241 + ID_242 + ID_243 + ID_244 + ID_245 + 
  ID_246 + ID_247 + ID_248 + ID_249 + ID_250 + ID_251 + ID_252 + 
  ID_253 + ID_254 + ID_255 + ID_256 + ID_257 + ID_258 + ID_259 + 
  ID_260 + ID_261 + ID_262 + ID_263 + ID_264 + ID_265 + ID_266 + 
  ID_267 + ID_268 + ID_269 + ID_270 + ID_271 + ID_272 + ID_273 + 
  ID_274 + ID_275 + ID_276 + ID_277 + ID_278 + ID_279 + ID_280 + 
  ID_281 + ID_282 + ID_283 + ID_284 + ID_285 + ID_286 + ID_287 + 
  ID_288 + ID_289 + ID_290 + ID_291 + ID_292 + ID_293 + ID_294 + 
  ID_295 + ID_296 + ID_297 + ID_298 + ID_299 + ID_300 + ID_301 + 
  ID_302 + ID_303 + ID_304 + ID_305 + ID_306 + ID_307 + ID_308 + 
  ID_309 + ID_310 + ID_311 + ID_312 + ID_313 + ID_314 + ID_315 + 
  ID_316 + ID_317 + ID_318 + ID_319 + ID_320 + ID_321 + ID_322 + 
  ID_323 + ID_324 + ID_325 + ID_326 + ID_327 + ID_328 + ID_329 + 
  ID_330 + ID_331 + ID_332 + ID_333 + ID_334 + ID_335 + ID_336 + 
  ID_337 + ID_338 + ID_339 + ID_340 + ID_341 + ID_342 + ID_343 + ID_344 + ID_345 + ID_346 + ID_347 + ID_348 + ID_349 + ID_350 + ID_351 + ID_352 + ID_353 + ID_354 + ID_355 + ID_356 + ID_357 + ID_358 + ID_359 + ID_360 + ID_361 + ID_362 + ID_363 + ID_364 + ID_365 + ID_366 + ID_367 + ID_368 + ID_369 + ID_370 + ID_371 + ID_372 + ID_373 + ID_374 + ID_375 + ID_376 + ID_377 + ID_378 + ID_379 + ID_380 + ID_381 + ID_382 + ID_383 + ID_384 + ID_385 + ID_386 + ID_387 + ID_388 + ID_389 + ID_390 + ID_391 + ID_392 + ID_393 + ID_394 + ID_395 + ID_396 + ID_397 + ID_398 + ID_399 + ID_400 + ID_401 + ID_402 + ID_403 + ID_404 + ID_405 + ID_406 + ID_407 + ID_408 + ID_409 + ID_410 + ID_411 + ID_412 + ID_413 + ID_414 + ID_415 + ID_416 + ID_417 + ID_418 + ID_419 + ID_420 + ID_421 + ID_422 + ID_423 + ID_424 + ID_425 + ID_426 + ID_427 + ID_428 + ID_429 + ID_430 + ID_431 + ID_432 + ID_433 + ID_434 + ID_435 + ID_436 + ID_437 + ID_438 + ID_439 + ID_440 + ID_441 + ID_442 + ID_443 + ID_444 + ID_445 + ID_446 + ID_447 + ID_448 + ID_449 + ID_450 + ID_451 + ID_452 + ID_453 + ID_454 + ID_455 + ID_456 + ID_457 + ID_458 + ID_459 + ID_460 + ID_461 + ID_462 + ID_463 + ID_464 + ID_465 + ID_466 + ID_467 + ID_468 + ID_469 + ID_470 + ID_471 + ID_472 + ID_473 + ID_474 + name(Secular, model = "ar", order=1)

fit2b <- bru(comp.a, formula = Collisions ~ Intercept + Time + Q1 + Q2 + Q3 + 
               pandemic + ID_1 + ID_2 + ID_3 + ID_4 + 
               ID_5 + ID_6 + ID_7 + ID_8 + ID_9 + ID_10 +  ID_11 + ID_12 + ID_13 + ID_14 + 
               ID_15 + ID_16 + ID_17 + ID_18 + ID_19 + ID_20 + ID_21 + ID_22 + ID_23 + 
               ID_24 + ID_25 + ID_26 + ID_27 + ID_28 + ID_29 + ID_30 + ID_31 + ID_32 + 
               ID_33 + ID_34 + ID_35 + ID_36 + ID_37 + ID_38 + ID_39 + ID_40 + ID_41 + 
               ID_42 + ID_43 + ID_44 + ID_45 + ID_46 + ID_47 + ID_48 + ID_49 + ID_50 + 
               ID_51 + ID_52 + ID_53 + ID_54 + ID_55 + ID_56 + ID_57 + ID_58 + ID_59 + ID_60 + 
               ID_61 + ID_62 + ID_63 + ID_64 + ID_65 + ID_66 + ID_67 + ID_68 + ID_69 + 
               ID_70 + ID_71 + ID_72 + ID_73 + ID_74 + ID_75 + ID_76 + ID_77 + ID_78 + 
               ID_79 + ID_80 + ID_81 + ID_82 + ID_83 + ID_84 + ID_85 + ID_86 + ID_87 + 
               ID_88 + ID_89 + ID_90 + ID_91 + ID_92 + ID_93 + ID_94 + ID_95 + ID_96 + 
               ID_97 + ID_98 + ID_99 + ID_100 + ID_101 + ID_102 + ID_103 + ID_104 + ID_105 + 
               ID_106 + ID_107 + ID_108 + ID_109 + ID_110 + ID_111 + ID_112 + ID_113 + 
               ID_114 + ID_115 + ID_116 + ID_117 + ID_118 + ID_119 + ID_120 + ID_121 + 
               ID_122 + ID_123 + ID_124 + ID_125 + ID_126 + ID_127 + ID_128 + ID_129 + 
               ID_130 + ID_131 + ID_132 + ID_133 + ID_134 + ID_135 + ID_136 + ID_137 + 
               ID_138 + ID_139 + ID_140 + ID_141 + ID_142 + ID_143 + ID_144 + ID_145 + 
               ID_146 + ID_147 + ID_148 + ID_149 + ID_150 + ID_151 + ID_152 + ID_153 + 
               ID_154 + ID_155 + ID_156 + ID_157 + ID_158 + ID_159 + ID_160 + ID_161 + 
               ID_162 + ID_163 + ID_164 + ID_165 + ID_166 + ID_167 + ID_168 + 
               ID_169 + ID_170 + ID_171 + ID_172 + ID_173 + ID_174 + ID_175 + 
               ID_176 + ID_177 + ID_178 + ID_179 + ID_180 + ID_181 + ID_182 + 
               ID_183 + ID_184 + ID_185 + ID_186 + ID_187 + ID_188 + ID_189 + 
               ID_190 + ID_191 + ID_192 + ID_193 + ID_194 + ID_195 + ID_196 + 
               ID_197 + ID_198 + ID_199 + ID_200 + ID_201 + ID_202 + ID_203 + 
               ID_204 + ID_205 + ID_206 + ID_207 + ID_208 + ID_209 + ID_210 + 
               ID_211 + ID_212 + ID_213 + ID_214 + ID_215 + ID_216 + ID_217 + 
               ID_218 + ID_219 + ID_220 + ID_221 + ID_222 + ID_223 + ID_224 + 
               ID_225 + ID_226 + ID_227 + ID_228 + ID_229 + ID_230 + ID_231 + 
               ID_232 + ID_233 + ID_234 + ID_235 + ID_236 + ID_237 + ID_238 + 
               ID_239 + ID_240 + ID_241 + ID_242 + ID_243 + ID_244 + ID_245 + 
               ID_246 + ID_247 + ID_248 + ID_249 + ID_250 + ID_251 + ID_252 + 
               ID_253 + ID_254 + ID_255 + ID_256 + ID_257 + ID_258 + ID_259 + 
               ID_260 + ID_261 + ID_262 + ID_263 + ID_264 + ID_265 + ID_266 + 
               ID_267 + ID_268 + ID_269 + ID_270 + ID_271 + ID_272 + ID_273 + 
               ID_274 + ID_275 + ID_276 + ID_277 + ID_278 + ID_279 + ID_280 + 
               ID_281 + ID_282 + ID_283 + ID_284 + ID_285 + ID_286 + ID_287 + 
               ID_288 + ID_289 + ID_290 + ID_291 + ID_292 + ID_293 + ID_294 + 
               ID_295 + ID_296 + ID_297 + ID_298 + ID_299 + ID_300 + ID_301 + 
               ID_302 + ID_303 + ID_304 + ID_305 + ID_306 + ID_307 + ID_308 + 
               ID_309 + ID_310 + ID_311 + ID_312 + ID_313 + ID_314 + ID_315 + 
               ID_316 + ID_317 + ID_318 + ID_319 + ID_320 + ID_321 + ID_322 + 
               ID_323 + ID_324 + ID_325 + ID_326 + ID_327 + ID_328 + ID_329 + 
               ID_330 + ID_331 + ID_332 + ID_333 + ID_334 + ID_335 + ID_336 + 
               ID_337 + ID_338 + ID_339 + ID_340 + ID_341 + ID_342 + ID_343 + ID_344 + ID_345 + ID_346 + ID_347 + ID_348 + ID_349 + ID_350 + ID_351 + ID_352 + ID_353 + ID_354 + ID_355 + ID_356 + ID_357 + ID_358 + ID_359 + ID_360 + ID_361 + ID_362 + ID_363 + ID_364 + ID_365 + ID_366 + ID_367 + ID_368 + ID_369 + ID_370 + ID_371 + ID_372 + ID_373 + ID_374 + ID_375 + ID_376 + ID_377 + ID_378 + ID_379 + ID_380 + ID_381 + ID_382 + ID_383 + ID_384 + ID_385 + ID_386 + ID_387 + ID_388 + ID_389 + ID_390 + ID_391 + ID_392 + ID_393 + ID_394 + ID_395 + ID_396 + ID_397 + ID_398 + ID_399 + ID_400 + ID_401 + ID_402 + ID_403 + ID_404 + ID_405 + ID_406 + ID_407 + ID_408 + ID_409 + ID_410 + ID_411 + ID_412 + ID_413 + ID_414 + ID_415 + ID_416 + ID_417 + ID_418 + ID_419 + ID_420 + ID_421 + ID_422 + ID_423 + ID_424 + ID_425 + ID_426 + ID_427 + ID_428 + ID_429 + ID_430 + ID_431 + ID_432 + ID_433 + ID_434 + ID_435 + ID_436 + ID_437 + ID_438 + ID_439 + ID_440 + ID_441 + ID_442 + ID_443 + ID_444 + ID_445 + ID_446 + ID_447 + ID_448 + ID_449 + ID_450 + ID_451 + ID_452 + ID_453 + ID_454 + ID_455 + ID_456 + ID_457 + ID_458 + ID_459 + ID_460 + ID_461 + ID_462 + ID_463 + ID_464 + ID_465 + ID_466 + ID_467 + ID_468 + ID_469 + ID_470 + ID_471 + ID_472 + ID_473 + ID_474 +  name,
             family = "nbinomial", data = Sim_pre.2)

summary(fit2b)

comp.t2 <- Collisions ~ Time+  Q1 + Q2+ Q3 + pandemic + factor(ID) + f(Secular, model = "ar",  order=1) 
fit2 <- inla(comp.t2, data=Sim_pre.2,
             family="nbinomial",  
             control.compute=list(dic=TRUE, config = TRUE, cpo=TRUE), 
             control.predictor=list(compute=FALSE),
             verbose = F)
summary(fit2) #TO COMPARE, the results must very similar with fit2b

### POST PERIOD
#ASE.PERIOD <- subset(df.2, Time >0, select =-c(Collisions))
#Note that I'm removing collisions now
ASE.PERIOD2 <- df.2; ASE.PERIOD2$Y <- ASE.PERIOD2$Collisions
ASE.PERIOD2 <-  subset(ASE.PERIOD2, select =-c(Collisions))

Count.Pred <- predict(fit2b, seed=123, ASE.PERIOD2, formula = ~ exp(Intercept + Time + Q1 + Q2 + Q3 + pandemic + ID_1 + ID_2 + ID_3 + ID_4 + 
                                                                      ID_5 + ID_6 + ID_7 + ID_8 + ID_9 + ID_10 +  ID_11 + ID_12 + ID_13 + ID_14 + 
                                                                      ID_15 + ID_16 + ID_17 + ID_18 +   ID_19 + ID_20 + ID_21 + ID_22 + ID_23 + 
                                                                      ID_24 + ID_25 + ID_26 +   ID_27 + ID_28 + ID_29 + ID_30 + ID_31 + ID_32 + 
                                                                      ID_33 + ID_34 +   ID_35 + ID_36 + ID_37 + ID_38 + ID_39 + ID_40 + ID_41 + 
                                                                      ID_42 + ID_43 + ID_44 + ID_45 + ID_46 + ID_47 + ID_48 + ID_49 + ID_50 + 
                                                                      ID_51 + ID_52 + ID_53 + ID_54 + ID_55 + ID_56 + ID_57 + ID_58 + ID_59 + ID_60 + 
                                                                      ID_61 + ID_62 + ID_63 + ID_64 + ID_65 + ID_66 + ID_67 + ID_68 + ID_69 + 
                                                                      ID_70 + ID_71 + ID_72 + ID_73 + ID_74 + ID_75 + ID_76 + ID_77 + ID_78 + 
                                                                      ID_79 + ID_80 + ID_81 + ID_82 + ID_83 + ID_84 + ID_85 + ID_86 + ID_87 + 
                                                                      ID_88 + ID_89 + ID_90 + ID_91 + ID_92 + ID_93 + ID_94 + ID_95 + ID_96 + 
                                                                      ID_97 + ID_98 + ID_99 + ID_100 + ID_101 + ID_102 + ID_103 + ID_104 + ID_105 + 
                                                                      ID_106 + ID_107 + ID_108 + ID_109 + ID_110 + ID_111 + ID_112 + ID_113 + 
                                                                      ID_114 + ID_115 + ID_116 + ID_117 + ID_118 + ID_119 + ID_120 + ID_121 + 
                                                                      ID_122 + ID_123 + ID_124 + ID_125 + ID_126 +  ID_127 + ID_128 + ID_129 + 
                                                                      ID_130 + ID_131 + ID_132 + ID_133 + ID_134 + ID_135 + ID_136 + ID_137 + 
                                                                      ID_138 + ID_139 + ID_140 + ID_141 + ID_142 + ID_143 + ID_144 + ID_145 + 
                                                                      ID_146 + ID_147 + ID_148 + ID_149 + ID_150 + ID_151 + ID_152 + ID_153 + 
                                                                      ID_154 + ID_155 + ID_156 + ID_157 + ID_158 + ID_159 + ID_160 + ID_161 + 
                                                                      ID_162 + ID_163 + ID_164 + ID_165 + ID_166 + ID_167 + ID_168 + 
                                                                      ID_169 + ID_170 + ID_171 + ID_172 + ID_173 + ID_174 + ID_175 + 
                                                                      ID_176 + ID_177 + ID_178 + ID_179 + ID_180 + ID_181 + ID_182 + 
                                                                      ID_183 + ID_184 + ID_185 + ID_186 + ID_187 + ID_188 + ID_189 + 
                                                                      ID_190 + ID_191 + ID_192 + ID_193 + ID_194 + ID_195 + ID_196 + 
                                                                      ID_197 + ID_198 + ID_199 + ID_200 + ID_201 + ID_202 + ID_203 + 
                                                                      ID_204 + ID_205 + ID_206 + ID_207 + ID_208 + ID_209 + ID_210 + 
                                                                      ID_211 + ID_212 + ID_213 + ID_214 + ID_215 + ID_216 + ID_217 + 
                                                                      ID_218 + ID_219 + ID_220 + ID_221 + ID_222 + ID_223 + ID_224 + 
                                                                      ID_225 + ID_226 + ID_227 + ID_228 + ID_229 + ID_230 + ID_231 + 
                                                                      ID_232 + ID_233 + ID_234 + ID_235 + ID_236 + ID_237 + ID_238 + 
                                                                      ID_239 + ID_240 + ID_241 + ID_242 + ID_243 + ID_244 + ID_245 + 
                                                                      ID_246 + ID_247 + ID_248 + ID_249 + ID_250 + ID_251 + ID_252 + 
                                                                      ID_253 + ID_254 + ID_255 + ID_256 + ID_257 + ID_258 + ID_259 + 
                                                                      ID_260 + ID_261 + ID_262 + ID_263 + ID_264 + ID_265 + ID_266 + 
                                                                      ID_267 + ID_268 + ID_269 + ID_270 + ID_271 + ID_272 + ID_273 + 
                                                                      ID_274 + ID_275 + ID_276 + ID_277 + ID_278 + ID_279 + ID_280 + 
                                                                      ID_281 + ID_282 + ID_283 + ID_284 + ID_285 + ID_286 + ID_287 + 
                                                                      ID_288 + ID_289 + ID_290 + ID_291 + ID_292 + ID_293 + ID_294 + 
                                                                      ID_295 + ID_296 + ID_297 + ID_298 + ID_299 + ID_300 + ID_301 + 
                                                                      ID_302 + ID_303 + ID_304 + ID_305 + ID_306 + ID_307 + ID_308 + 
                                                                      ID_309 + ID_310 + ID_311 + ID_312 + ID_313 + ID_314 + ID_315 + 
                                                                      ID_316 + ID_317 + ID_318 + ID_319 + ID_320 + ID_321 + ID_322 + 
                                                                      ID_323 + ID_324 + ID_325 + ID_326 + ID_327 + ID_328 + ID_329 + 
                                                                      ID_330 + ID_331 + ID_332 + ID_333 + ID_334 + ID_335 + ID_336 + 
                                                                      ID_337 + ID_338 + ID_339 + ID_340 + ID_341 + ID_342 + ID_343 + ID_344 + ID_345 + ID_346 + ID_347 + ID_348 + ID_349 + ID_350 + ID_351 + ID_352 + ID_353 + ID_354 + ID_355 + ID_356 + ID_357 + ID_358 + ID_359 + ID_360 + ID_361 + ID_362 + ID_363 + ID_364 + ID_365 + ID_366 + ID_367 + ID_368 + ID_369 + ID_370 + ID_371 + ID_372 + ID_373 + ID_374 + ID_375 + ID_376 + ID_377 + ID_378 + ID_379 + ID_380 + ID_381 + ID_382 + ID_383 + ID_384 + ID_385 + ID_386 + ID_387 + ID_388 + ID_389 + ID_390 + ID_391 + ID_392 + ID_393 + ID_394 + ID_395 + ID_396 + ID_397 + ID_398 + ID_399 + ID_400 + ID_401 + ID_402 + ID_403 + ID_404 + ID_405 + ID_406 + ID_407 + ID_408 + ID_409 + ID_410 + ID_411 + ID_412 + ID_413 + ID_414 + ID_415 + ID_416 + ID_417 + ID_418 + ID_419 + ID_420 + ID_421 + ID_422 + ID_423 + ID_424 + ID_425 + ID_426 + ID_427 + ID_428 + ID_429 + ID_430 + ID_431 + ID_432 + ID_433 + ID_434 + ID_435 + ID_436 + ID_437 + ID_438 + ID_439 + ID_440 + ID_441 + ID_442 + ID_443 + ID_444 + ID_445 + ID_446 + ID_447 + ID_448 + ID_449 + ID_450 + ID_451 + ID_452 + ID_453 + ID_454 + ID_455 + ID_456 + ID_457 + ID_458 + ID_459 + ID_460 + ID_461 + ID_462 + ID_463 + ID_464 + ID_465 + ID_466 + ID_467 + ID_468 + ID_469 + ID_470 + ID_471 + ID_472 + ID_473 + ID_474 + name_eval(Secular)), n.samples=1000)  

Verif <- subset(Count.Pred, Time <1, select =c(ID, Time, Y, mean, sd, mean.mc_std_err, sd.mc_std_err, q0.025, q0.975, median)) # GOOD!

columns_to_sum <- c("Y", "mean", "q0.025", "median", "q0.975")
Verif2 <- Verif %>%
  group_by(Time) %>%
  summarise(across(all_of(columns_to_sum), sum))
Verif2$N <- 475

Post.1 <- subset(Count.Pred, Time ==1, select =c(ID, Time, mean, sd, mean.mc_std_err, sd.mc_std_err, q0.025, q0.975, median))

columns_to_sum <- c("mean", "q0.025", "median", "q0.975")
Post.1.Sum <- Post.1 %>%
  summarise(across(all_of(columns_to_sum), sum, na.rm = TRUE)) #THIS IS THE RESULT WE WANT
Post.1.Sum$N <- 475; Post.1.Sum$Time <- 1;

#For post 2 we can sample 275 locations, first we select TIME =2
Post.2 <- subset(Count.Pred, Time ==2, select =c(ID, Time, mean, sd, mean.mc_std_err, sd.mc_std_err, q0.025, q0.975, median))
#second we sample, set.seed For reproducibility
set.seed(123)
Post.2a <- Post.2[sample(nrow(Post.2), 401), ]
Post.2.Sum <- Post.2a %>%
  summarise(across(all_of(columns_to_sum), sum, na.rm = TRUE)) #THIS IS THE RESULT WE WANT
Post.2.Sum$N <- 401; Post.2.Sum$Time <- 2; 

#For post 3 we can sample 125 locations, first we select TIME =3
Post.3 <- subset(Count.Pred, Time==3, select =c(ID, Time, mean, sd, mean.mc_std_err, sd.mc_std_err, q0.025, q0.975, median))
#second we sample, set.seed For reproducibility
set.seed(123)
Post.3a <- Post.3[Post.3$ID %in% sample(nrow(Post.2), 201), ] # Note that I'm sampling from the 200 previous locations (Post.2)
Post.3.Sum <- Post.3a %>%
  summarise(across(all_of(columns_to_sum), sum, na.rm = TRUE)) #THIS IS THE RESULT WE WANT
Post.3.Sum$N <- 201; Post.3.Sum$Time <- 3

#For post 4 we can sample 50 locations, first we select TIME =4
Post.4 <- subset(Count.Pred, Time ==4, select =c(ID, Time, mean, sd, mean.mc_std_err, sd.mc_std_err, q0.025, q0.975, median))
#second we sample, set.seed For reproducibility
set.seed(123)
Post.4a <- Post.4[Post.4$ID %in% sample(nrow(Post.3), 118), ] # Note that I'm sampling from the 100 previous locations (Post.3)
Post.4.Sum <- Post.4a %>%
  summarise(across(all_of(columns_to_sum), sum, na.rm = TRUE)) #THIS IS THE RESULT WE WANT
Post.4.Sum$N <- 118;  Post.4.Sum$Time <- 4

RESULTS <- rbind.data.frame(Post.1.Sum, Post.2.Sum, Post.3.Sum, Post.4.Sum)

RESULTS$Lower_DIFF <- RESULTS$q0.025 - RESULTS$median
RESULTS$Upper_DIFF <- RESULTS$q0.975 - RESULTS$median

RESULTS$Lower_PERC <- RESULTS$q0.025 / RESULTS$median
RESULTS$Upper_PERC <- RESULTS$q0.975 / RESULTS$median
write.csv(RESULTS,"RESULTS.csv", row.names = FALSE) 

# Find the columns in Verif2 that are not in RESULTS
missing_columns <- setdiff(names(Verif2), names(RESULTS))
# Add these columns to RESULTS with NA values
RESULTS[missing_columns] <- NA
# Find the columns in RESULTS that are not in Verif2
missing_columns2 <- setdiff(names(RESULTS), names(Verif2))
# Add these columns to Verif2 with NA values
Verif2[missing_columns2] <- NA

# Now rbind the two datasets
Results_2 <- rbind(Verif2, RESULTS)

names(Results_2)[names(Results_2) == "Y"] <- "Collisions"
Results_2$CI_ratio <- Results_2$q0.975 / Results_2$q0.025