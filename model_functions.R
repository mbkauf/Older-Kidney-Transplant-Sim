`" 
*-------------------------------------------------------------------------------
* Simulation model functions
* model_functions.R
* Date Created: 01/27/2022
* Last Updated: 04/26/2022
* Matt Kaufmann, Stanford University
* Note: Modified to run on remote server
*-------------------------------------------------------------------------------
"`

library(fitdistrplus)
library(MASS)
library(copula)
library(psych)
library(dplyr)
library(survival)
library(Hmisc)

## Function to build simulation dataset
build.data <- function(seed, n) {
  ### Function Inputs
  # seed: used to replicate results
  # n: sample size used for simulation
  ### Create bootstrapped data set
  set.seed(seed)
  num_patients <- n
  
  ## Create the bootstrapped dataset and keep the names of the each variable
  boot_data <<- data[sample(nrow(data), num_patients, replace=T),]
  sim_data <- boot_data[,-1]
  rownames(sim_data) <- boot_data[,1]
  
  ## For analysis, we need to store time to event variables
  event             <<- boot_data[,"event"]
  start_time        <<- boot_data[,"start_time"]
  months_to_event   <<- boot_data[,"months_to_event"]
  event_cd          <<- boot_data[,"event_cd"]
  dgf               <<- boot_data[,"dgf"]
  graft_loss        <<- boot_data[,"graft_loss"]
  months_to_gl      <<- boot_data[,"months_to_gl"]  
  gs_death          <<- boot_data[,"gs_death"]
  months_to_gs_mort <<- boot_data[,"months_to_gs_mort"]
  gl_death          <<- boot_data[,"gl_death"]
  months_to_gl_mort <<- boot_data[,"months_to_gl_mort"]
  kdpi              <<- boot_data[,"kdpi"]
  tx_outcome        <<- boot_data[,"tx_outcome"]
  
  ## Remove variables not to be used in the analysis
  boot_data <<- boot_data[,-c(16:26, 37, 38)]
  
  ## For analysis by race/ethnicity
  data_list_white    <<- subset(boot_data, 
                                boot_data[,3]==0 & boot_data[,4]==0 & boot_data[,5]==0)
  data_list_black    <<- subset(boot_data, 
                                boot_data[,3]==1 & boot_data[,4]==0 & boot_data[,5]==0)
  data_list_hispanic <<- subset(boot_data, 
                                boot_data[,3]==0 & boot_data[,4]==1 & boot_data[,5]==0)
  data_list_other    <<- subset(boot_data, 
                                boot_data[,3]==0 & boot_data[,4]==0 & boot_data[,5]==1)
  
  ## Sorted bootstrapped data by race/ethnicity
  boot_data          <<- rbind(data_list_white, data_list_black, 
                               data_list_hispanic, data_list_other)
  
  ## Remove race variables
  data_list_white    <<- data_list_white[,-c(3, 4, 5)]
  data_list_black    <<- data_list_black[,-c(3, 4, 5)]
  data_list_hispanic <<- data_list_hispanic[,-c(3, 4, 5)]
  data_list_other    <<- data_list_other[,-c(3, 4, 5)]
  
  ## Cleaning up observed dataset for analysis 
  obs_df    <<- as.data.frame(data)
  obs_df    <<- obs_df %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0, 
                         ifelse(black==1 & hispanic==0 & other_race==0, 1, 
                                ifelse(black==0 & hispanic==1 & 
                                         other_race==0, 2, 3)))) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black", 
                                    "Hispanic", "Other")))
}
 
## Function to run waitlist simulation
list.simulation <- function(seed, n, 
                            data = sample_df, 
                            ddtx_coef = ddtx_coef_w, ddtx_p = ddtx_p_w, 
                            ldtx_coef = ldtx_coef_g, ldtx_gamma = ldtx_gamma_g,
                            mort_coef = mort_coef_w, mort_shape = mort_p_w,
                            remove_coef = remove_coef_w, remove_p = remove_p_w,
                            ddtx_rate = 0) {
  ### Function Inputs
  # seed: used to replicate results
  # n: sample size used for simulation
  # t: number of months to run simulation over
  ### Create bootstrapped data set
  set.seed(seed)
  num_patients <- n
  
  ### Set up wait list equations
  # Deceased Donor Tx (AFT Weibull)
  ddtx_lambda <- exp(-exp(ddtx_p) * data %*% ddtx_coef)
  
  # Living Donor Tx (PH Gompertz)
  ldtx_lambda <- exp(data %*% ldtx_coef)
  
  # Wait list mortality (AFT Weibull)
  mort_lambda <- exp(-exp(mort_shape) * data %*% mort_coef)

  # Other wait list removal (AFT Weibull)
  remove_lambda <- exp(-exp(remove_p) * data %*% remove_coef)
  
  ### Initialize vectors and matrices
  event_results   <- rep(0, num_patients)
  
  ddtx_results    <- rep(0, num_patients)
  
  ldtx_results    <- rep(0, num_patients)
  
  mort_results    <- rep(0, num_patients)
  
  remove_results  <- rep(0, num_patients)
  
  m2event         <- rep(0, num_patients)
  events          <- rep(0, num_patients)
  
  m2age100        <- (rep(35, num_patients) - data[,1]) * 12

  # Generate vector of months until candidate is 65
  # m2_65           <- boot_data[,1]
  # m2_65           <- replace(m2_65, m2_65>0, 0)
  # u65             <- replace(m2_65, m2_65<0, 1)
  # m2_65           <- (m2_65+0.5)*-12
  # m2_65           <- m2_65 * u65
  
  ### Run Simulation
  i <- 0
  while (all(event_results == 1) == FALSE) {
    i = i + 1
    ### Deceased Donor Tx?
    # Generate vector of ddtx probabilities
    ddtx_haz  <- (exp(ddtx_p) * ddtx_lambda * i^(exp(ddtx_p) - 1))*(1+ddtx_rate)
    ddtx_prob <- 1 - exp(-ddtx_haz)
    
    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_ddtx_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y              <- as.integer(rand_ddtx_surv < ddtx_prob)  # vector of those who had event in period i
    new_ddtx       <- as.integer(y==1 & pmax(y, ddtx_results) > ddtx_results)  # vector of new events
    
    # Check check that they are still at risk for event
    ddtx_results <- pmax(y, ddtx_results)
    match        <- as.integer(ddtx_results == pmax(ldtx_results, mort_results, 
                                                    remove_results) & 
                                 pmax(ldtx_results, mort_results, remove_results)==1)
    ddtx_results <- ddtx_results - match  # remove those who were not at risk for ddtx
    new_ddtx     <- new_ddtx - match  # final vector of new ddtx events
    
    ### Living Donor Tx?
    # Generate vector of ldtx probabilities
    ldtx_prob <- 1 - exp(-(ldtx_lambda * exp(ldtx_gamma*i)))  # Graft loss probability vector
    
    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_ldtx_surv <- runif(num_patients, min = 0, max = 1)  # random number vector
    x              <- as.integer(rand_ldtx_surv < ldtx_prob)  # vector of possible events
    new_ldtx       <- as.integer(x==1 & pmax(x, ldtx_results) > ldtx_results)
    
    # Check check that they are still at risk for event
    ldtx_results <- pmax(x, ldtx_results)
    match        <- as.integer(ldtx_results == pmax(ddtx_results, mort_results, 
                                                    remove_results) & 
                                 pmax(ddtx_results, mort_results, remove_results)==1)
    ldtx_results <- ldtx_results - match
    new_ldtx     <- new_ldtx - match    

    ### Wait list mortality?
    mort_prob <- 1 - exp(-(exp(mort_shape) * mort_lambda * i^(exp(mort_shape) - 1)))

    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_mort_surv <- runif(num_patients, min = 0, max = 1)  # random number vector
    x              <- as.integer(rand_mort_surv < mort_prob)  # vector of possible events
    age100         <- as.integer(i >= m2age100)      # force death at age 100
    x              <- pmax(x, age100)   # new vector of possible events
    new_mort       <- as.integer(x==1 & pmax(x, mort_results) > mort_results)
    
    # Check check that they are still at risk for event
    mort_results <- pmax(x, mort_results)
    match        <- as.integer(mort_results == pmax(ddtx_results, ldtx_results, 
                                                    remove_results) & 
                                 pmax(ddtx_results, ldtx_results, remove_results)==1)
    mort_results <- mort_results - match
    new_mort     <- new_mort - match    
    
    ### Other List Removal?
    # Generate vector of list removal probabilities
    remove_prob <- 1 - exp(-(exp(remove_p) * remove_lambda * i^(exp(remove_p) - 1)))
    
    # Determine who has event by using random uniform draws and determine who
    # had a new event this period
    rand_remove_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y                <- as.integer(rand_remove_surv < remove_prob)  # vector of those who had event in period i
    new_remove       <- as.integer(y==1 & pmax(y, remove_results) > remove_results)
    
    # Check check that they are still at risk for event
    remove_results <- pmax(y, remove_results)
    match          <- as.integer(remove_results == pmax(ddtx_results, 
                                                        ldtx_results, mort_results) 
                                 & pmax(ddtx_results, ldtx_results, mort_results)==1)
    remove_results <- remove_results - match
    new_remove     <- new_remove - match 
    
    # Create vector the tracks cumulative removal events, create patient trace
    # remove_trace[,i] <- remove_results
    
    ### Time to any event
    event_results <- pmax(ddtx_results, ldtx_results,
                          mort_results, remove_results)
    m2event       <- m2event + pmax(new_ddtx, new_ldtx, 
                                    new_mort, new_remove)*i 
  }
  
  ### Build dataset for analysis
  results_sim <- as.data.frame(cbind(data, event_results, ddtx_results,
                                     ldtx_results, mort_results, remove_results,
                                     m2event))
  results_sim <- results_sim %>%
    rename(
      event = event_results,
      ddtx = ddtx_results,
      ldtx = ldtx_results,
      mort = mort_results,
      remove = remove_results,
      months_to_event = m2event
    ) %>%
    mutate(months_to_event = ifelse(event==0, months, months_to_event)) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0, 
                         ifelse(black==1 & hispanic==0 & other_race==0, 1, 
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(sim = rep(1, num_patients)) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black", 
                                    "Hispanic", "Other"))) %>%
    mutate(sim = factor(sim, levels = c(0, 1), 
                        labels = c("Actual", "Simulated")))

  list_sim_df <<- results_sim
}

## Function to run post-tx simulation
run.simulation <- function(seed, data = list_sim_df, 
                           gl30_coef = gl30_coef_ml_noreg, 
                           dgf_coef = dgf_coef_ml_noreg,
                           mort30_coef = mort30_coef_ml_noreg, 
                           gl_coef = gl_coef_g_noreg,
                           gl_shape = gl_gamma_g_noreg, 
                           gs_mort_coef = gs_mort_coef_g_noreg,
                           gs_mort_shape = gs_mort_gamma_g_noreg, 
                           dial_mort_coef = dial_mort_coef_w_noreg, 
                           dial_mort_shape = dial_mort_lnp_w_noreg, 
                           gl_mort_coef = gl_mort_coef_logit_noreg) {
  ### Function Inputs
  # seed: used to replicate results
  # n: sample size used for simulation
  # t: number of months to run simulation over
  ### Create bootstrapped data set
  set.seed(seed)
  # months <- t
  
  ### Start first 30 days
  ### Make the following changes to the data used for the next part
  # keep only deceased donor transplants
  # update to age at transplant
  # update to time on dialysis prior to transplant
  # update list year to transplant year
  # draw from joint distribution of cold ischemia time and KDPI
  # reorder variables to match order for transplant outcomes
  invisible(rand_kdpi_cold(n = nrow(data)))
  data <- cbind(data, Z)
  names(data)[ncol(data)] <- "kdpi"
  names(data)[ncol(data)-1] <- "rec_cold_isch_tm"
  
  ddtx            <- data$ddtx
  ldtx            <- data$ldtx
  mort            <- data$mort
  remove          <- data$remove
  months_to_event <- data$months_to_event

  optn_reg2  <- data$optn_reg2
  optn_reg3  <- data$optn_reg3
  optn_reg4  <- data$optn_reg4
  optn_reg5  <- data$optn_reg5
  optn_reg6  <- data$optn_reg6
  optn_reg7  <- data$optn_reg7
  optn_reg8  <- data$optn_reg8
  optn_reg9  <- data$optn_reg9
  optn_reg10 <- data$optn_reg10
  optn_reg11 <- data$optn_reg11
    
  data <- data %>%
    mutate(rec_age_at_tx_c = can_age_at_listing_c + (months_to_event/12)) %>%
    mutate(years_dial = baseline_yrs_dial + (months_to_event/12)) %>%
    mutate(tx_year_c = list_year_c + (months_to_event/12)) %>%
    mutate(tx_year_c = ifelse(tx_year_c > 10, 10, tx_year_c)) %>%
    mutate(constant = 1) %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race, 
                  years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b, 
                  can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat, 
                  copd, pvd,	ang_cad, constant) 
  
  
  data <- as.matrix(data)
  num_patients <- nrow(data)
  
  ### Which 30-Day event? (Multinomial Logit)
  # Base Outcome: Graft Success
  gl30_exb   <- exp(data %*% gl30_coef)  
  dgf_exb    <- exp(data %*% dgf_coef)
  mort30_exb <- exp(data %*% mort30_coef)
  
  gl30_prob     <- gl30_exb / (1 + gl30_exb + dgf_exb + mort30_exb)
  dgf_prob      <- dgf_exb / (1 + gl30_exb + dgf_exb + mort30_exb)
  mort30_prob   <- mort30_exb / (1 + gl30_exb + dgf_exb + mort30_exb)
  success_prob  <- 1 / (1 + gl30_exb + dgf_exb + mort30_exb)
  
  gl30_range    <- gl30_prob
  dgf_range     <- gl30_prob + dgf_prob
  mort30_range  <- dgf_prob + mort30_prob
  
  dgf_results     <- rep(0, num_patients)
  mort30_results  <- rep(0, num_patients)
  success_results <- rep(0, num_patients)
  tx_outcome      <- rep(0, num_patients)
  
  rand_num_v  <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
  
  gl30_results    <- as.integer(rand_num_v <= gl30_range)
  gl_results      <- as.integer(rand_num_v <= gl30_range)
  dgf_results     <- as.integer((rand_num_v > gl30_range) & 
                                  (rand_num_v <= dgf_range))
  mort30_results  <- as.integer((rand_num_v > dgf_range) & 
                                  (rand_num_v <= mort30_range))
  success_results <- as.integer(rand_num_v > mort30_range)
  
  tx_outcome <- tx_outcome + gl30_results + (dgf_results*2) + (mort30_results*3)
    
  data_gl <- as.data.frame(cbind(data, dgf_results))
  

  data_gl_mort <- data_gl %>%
    mutate(m2gl = rep(0, nrow(data_gl))) %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race, 
                  years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b, 
                  can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat, copd, 
                  pvd, ang_cad, dgf_results, m2gl, constant)
  data_gs_mort <- data_gl %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race, 
                  years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
                  can_blood_o, rec_cold_isch_tm, diab_stat, copd, 
                  pvd, ang_cad, dgf_results, constant) # tx_year_c
  data_gl <- data_gl %>%
    dplyr::select(rec_age_at_tx_c,	male,	black, hispanic, other_race, 
                  years_dial, canhx_cpra, kdpi, can_blood_ab, can_blood_b,
                  can_blood_o, rec_cold_isch_tm, tx_year_c, diab_stat, copd, 
                  pvd, ang_cad, dgf_results, constant)
  
  data_gl_mort <- as.matrix(data_gl_mort)
  data_gl      <- as.matrix(data_gl)
  data_gs_mort <- as.matrix(data_gs_mort)
  
  # build.data(n=num_patients)
  ### Graft loss (PH Gompertz)
  gl_lambda <- exp(data_gl %*% gl_coef)
  
  ### Death w/ functioning graft (PH Gompertz)
  gs_mort_lambda <- exp(data_gs_mort %*% gs_mort_coef)
  
  ### Initialize vectors and matrices
  gs_mort_results <- rep(0, num_patients)
  gl_results      <- rep(0, num_patients)
  gl_mort_results <- rep(0, num_patients)
  
  m2gl      <- rep(0, num_patients)
  m2gs_mort <- rep(0, num_patients)
  m2gl_mort <- rep(0, num_patients)
  
  death    <- rep(0, num_patients)
  m2age100 <- (rep(35, num_patients) - data_gs_mort[,1]) * 12
  
  ### Run Simulation
  i <- 0
  while (all(death == 1) == FALSE) {
    i = i + 1
    ### First we will see who has died w/ functioning graft
    gs_mort_haz <- gs_mort_lambda * exp(gs_mort_shape*i) # vector of hazards
    gs_mort_prob <- 1 - exp(-gs_mort_haz)  # Graft loss probability vector
    rand_gs_mort_surv <- runif(num_patients, min = 0, max = 1)  # random number vector
    x <- as.integer(rand_gs_mort_surv < gs_mort_prob)  # vector of possible events
    age100         <- as.integer(i >= m2age100)      # force death at age 100
    x              <- pmax(x, age100)   # new vector of possible events
    new_gs_mort <- as.integer(x==1 & pmax(x, gs_mort_results) > gs_mort_results)
    gs_mort_results <- pmax(x, gs_mort_results)   # create vector that combines those with previous deaths to new deaths
    match <- as.integer(gs_mort_results == 1 & 
                          (gl_results==1 | mort30_results==1 | gl30_results==1))
    gs_mort_results <- gs_mort_results - match
    new_gs_mort <- new_gs_mort - match
    m2gs_mort <- m2gs_mort + new_gs_mort*i
    # gs_mort_trace[,i] <- gs_mort_results  
    
    ### Next we see who died after graft failure
    data_gl_mort[,ncol(data_gl_mort) - 1] <- m2gl # Update data set
    dial_mort_lambda <- exp(-exp(dial_mort_shape) * data_gl_mort %*% dial_mort_coef) # Update lambda vector
    dial_mort_haz <- exp(dial_mort_shape) * dial_mort_lambda * i^(exp(dial_mort_shape) - 1)  # Graft loss hazard vector
    dial_mort_prob <- 1 - exp(-dial_mort_haz)  # Graft loss probability vector
    rand_dial_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    z <- as.integer(rand_dial_surv < dial_mort_prob)  # vector of those who had event in period i
    age100         <- as.integer(i >= m2age100)      # force death at age 100
    Z              <- pmax(x, age100)   # new vector of possible events
    new_gl_mort <- as.integer(z==1 & pmax(z, gl_mort_results) > gl_mort_results)
    gl_mort_results <- pmax(z, gl_mort_results)  # create vector that combines those with previous graft loss to new graft loss
    gl_mort_results <- gl_results * gl_mort_results
    new_gl_mort <- gl_results * new_gl_mort
    m2gl_mort <- m2gl_mort + new_gl_mort*i 

    gl_haz <- gl_lambda * exp(gl_shape*i) # vector of hazards (Gompertz)
    gl_prob <- 1 - exp(-gl_haz)  # Graft loss probability vector
    
    rand_gl_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    y <- as.integer(rand_gl_surv < gl_prob)  # vector of those who had event in period i
    new_gl <- as.integer(y==1 & pmax(y, gl_results) > gl_results)
    gl_results <- pmax(y, gl_results)  # create vector that combines those with previous graft loss to new graft loss
    match <- as.integer(gl_results==1 & 
                          (gs_mort_results==1 | mort30_results==1))
    gl_results <- gl_results - match
    new_gl <- new_gl - match
    m2gl <- m2gl + new_gl*i
 
    
    ### Day of graft loss death
    gl_mort_prob <- exp(data_gl_mort %*% gl_mort_coef) / 
      (1 + exp(data_gl_mort %*% gl_mort_coef))
    rand_gl_mort_surv <- runif(num_patients, min = 0, max = 1)  # Vector of random probabilities from uniform
    a <- as.integer(rand_gl_mort_surv < gl_mort_prob)  # vector of those who had event in period i
    new_gl_mort <- as.integer(a==1 & pmax((new_gl *a), gl_mort_results) > gl_mort_results)
    gl_mort_results <- pmax(new_gl * a, gl_mort_results)
    m2gl_mort <- m2gl_mort + new_gl_mort*i 
    
    # Is anyone still alive?
    death = pmax(gs_mort_results, gl_mort_results, age100)
  }
  
  ### Build dataset for analysis
  list_vars_tx = c("graft_loss", "gs_death", "gl_death", "gl_time", 
                   "gs_death_time", "gl_death_time", "dgf",
                   "kdpi", "rec_cold_isch_tm")

  sim_df <<- as.data.frame(cbind(data_gl, ddtx, ldtx, mort, remove, 
                                 months_to_event, gl_results, gs_mort_results, 
                                 gl_mort_results, m2gl, m2gs_mort, m2gl_mort,
                                 optn_reg2, optn_reg3, optn_reg4, optn_reg5,
                                 optn_reg6, optn_reg7, optn_reg8, optn_reg9,
                                 optn_reg10, optn_reg11, tx_outcome))
  
  sim_df <<- sim_df %>%
    rename(
      graft_loss = gl_results,
      gs_death = gs_mort_results,
      gl_death = gl_mort_results,
      gl_time = m2gl,
      gs_death_time = m2gs_mort,
      gl_death_time = m2gl_mort,
      dgf = dgf_results
    ) %>%
    mutate_at(.vars = list_vars_tx, .funs = list(~ifelse(ddtx==1, ., NA))) %>%
    mutate(gl_time = ifelse(graft_loss==0 & gs_death==0,
                            i, gl_time)) %>%
    mutate(gl_time = ifelse(graft_loss==0 & gs_death==1,
                            gs_death_time, gl_time)) %>%
    mutate(gs_death_time = ifelse(gs_death==0 & graft_loss==0,
                                  i, gs_death_time)) %>%
    mutate(gs_death_time = ifelse(gs_death==0 & graft_loss==1,
                                  gl_time, gs_death_time)) %>%
    mutate(gl_death_time = gl_death_time - gl_time) %>%
    mutate(gl_death_time = ifelse(gl_death==0 & graft_loss==1,
                                  -gl_time+i, gl_death_time)) %>%
    mutate(gl_death_time = ifelse(gl_death_time < 0, NA, gl_death_time)) %>%
    mutate(race = ifelse(black==0 & hispanic==0 & other_race==0, 0,
                         ifelse(black==1 & hispanic==0 & other_race==0, 1,
                                ifelse(black==0 & hispanic==1 & other_race==0, 2, 3)))) %>%
    mutate(sim = rep(1, num_patients)) %>%
    mutate(race = factor(race, levels = c(0, 1, 2, 3),
                         labels = c("White", "Black",
                                    "Hispanic", "Other"))) %>%
    mutate(sim = factor(sim, levels = c(0, 1),
                        labels = c("Observed", "Simulated"))) %>%
    mutate(tx_outcome = factor(tx_outcome, levels = c(0, 1, 2, 3),
                         labels = c("GS", "GL",
                                    "DGF", "Mort")))
  
}

## Function to generate random values for KDPI and cold ischemia time
rand_kdpi_cold <- function(seed=11111, n, data=kdpi_cold, france=FALSE) {
  set.seed(seed)
  if (france==TRUE) {
    kdpi_shape1 = 1.045481
    kdpi_shape2 = 0.697808
  } else {
    kdpi_shape1 = 1.458402
    kdpi_shape2 = 1.033602
  }
  myCop <- normalCopula(param=0.130092, dim = 2, dispstr = "un")
  myMvd <- mvdc(copula=myCop, margins=c("gamma", "beta"),
                paramMargins=list(list(shape = 3.949532, 
                                       rate = 0.2128366),
                                  list(shape1 = kdpi_shape1, 
                                       shape2 = kdpi_shape2)))
  Z <<- rMvdc(n, myMvd)
  colnames(Z) <<- c("rec_cold_isch_tm", "kdpi")
  Z[,2] <<- Z[,2] * 100
  Z[,2] <<- round(Z[,2],0)
}

## Function to run code with output suppressed
hush=function(code){
  sink("NUL") # use /dev/null in UNIX
  tmp = code
  sink()
  return(tmp)
}

## Function to vary parameters for survival equations
chol_coef_surv <- function(cov_mat, v_coef, shape) {
  chol_mat <- chol(cov_mat)
  v_in <- qnorm(runif(nrow(cov_mat), 0, 1), 0, 1)  # vector of random inverse normal
  v_tz <- rep(0, nrow(cov_mat))
  
  for (i in 1:nrow(cov_mat)) {
    v_tz[i] <- v_in %*% chol_mat[i,]
  }
  shape_tz <- v_tz[nrow(cov_mat)]
  coef_tz <- v_tz[1: nrow(cov_mat)-1]
  coef_psa <<- v_coef + coef_tz
  shape_psa <<- shape + shape_tz
}

## Function to vary parameters for non-survival equations
chol_coef_logit <- function(cov_mat, v_coef) {
  chol_mat <- chol(cov_mat)
  v_in <- qnorm(runif(nrow(cov_mat), 0, 1), 0, 1)  # vector of random inverse normal
  v_tz <- rep(0, nrow(cov_mat))
  
  for (i in 1:nrow(cov_mat)) {
    v_tz[i] <- v_in %*% chol_mat[i,]
  }
  coef_tz <- v_tz
  coef_psa <<- v_coef + coef_tz
}

## Function to get probabilistically varied coefficients for waitlist equations
psa_list <- function(ddtx_cov_mat = ddtx_cov, ddtx_v_coef = ddtx_coef_w, 
                     ddtx_shape = ddtx_p_w, ldtx_cov_mat = ldtx_cov, 
                     ldtx_v_coef = ldtx_coef_g, ldtx_shape = ldtx_gamma_g,
                     mort_cov_mat = mort_cov, mort_v_coef = mort_coef_w, 
                     mort_shape = mort_p_w, remove_cov_mat = remove_cov, 
                     remove_v_coef = remove_coef_w, remove_shape = remove_p_w) {
  chol_coef_surv(cov_mat = ddtx_cov_mat, v_coef = ddtx_v_coef, shape = ddtx_shape)
  ddtx_coef_psa <<- coef_psa
  ddtx_shape_psa <<- shape_psa
  
  chol_coef_surv(cov_mat = ldtx_cov_mat, v_coef = ldtx_v_coef, shape = ldtx_shape)
  ldtx_coef_psa <<- coef_psa
  ldtx_shape_psa <<- shape_psa
  
  chol_coef_surv(cov_mat = mort_cov_mat, v_coef = mort_v_coef, shape = mort_shape)
  mort_coef_psa <<- coef_psa
  mort_shape_psa <<- shape_psa
  
  chol_coef_surv(cov_mat = remove_cov_mat, v_coef = remove_v_coef, shape = remove_shape)
  remove_coef_psa <<- coef_psa
  remove_shape_psa <<- shape_psa
  
}

## Function to get probabilistically varied coefficients for post-tx equations
psa_posttx <- function(gl_cov_mat = gl_cov, 
                       gl_v_coef = gl_coef_g_noreg, 
                       gl_shape = gl_gamma_g_noreg, 
                       gs_mort_cov_mat = gs_mort_cov, 
                       gs_mort_v_coef = gs_mort_coef_g_noreg, 
                       gs_mort_shape = gs_mort_gamma_g_noreg,
                       dial_mort_cov_mat = dial_mort_cov,
                       dial_mort_v_coef = dial_mort_coef_w_noreg, 
                       dial_mort_shape = dial_mort_lnp_w_noreg,
                       gl_mort_cov_mat = gl_mort_cov,
                       gl_mort_v_coef = gl_mort_coef_logit_noreg,
                       mlogit_cov_mat = mlogit_cov,
                       mlogit_v_coef = rbind(gl30_coef_ml_noreg, 
                                             dgf_coef_ml_noreg, 
                                             mort30_coef_ml_noreg)) {
  chol_coef_surv(cov_mat = gl_cov_mat, v_coef = gl_v_coef, shape = gl_shape)
  gl_coef_psa <<- coef_psa
  gl_shape_psa <<- shape_psa
  
  chol_coef_surv(cov_mat = gs_mort_cov_mat, v_coef = gs_mort_v_coef, 
                 shape = gs_mort_shape)
  gs_mort_coef_psa <<- coef_psa
  gs_mort_shape_psa <<- shape_psa
  
  chol_coef_surv(cov_mat = dial_mort_cov_mat, v_coef = dial_mort_v_coef, 
                 shape = dial_mort_shape)
  dial_mort_coef_psa <<- coef_psa
  dial_mort_shape_psa <<- shape_psa
  
  chol_coef_logit(cov_mat = gl_mort_cov_mat, v_coef = gl_mort_v_coef)
  gl_mort_coef_psa <<- coef_psa
  
  chol_coef_logit(cov_mat = mlogit_cov_mat, v_coef = mlogit_v_coef)
  mlogit_coef_psa <<- coef_psa
  gl30_coef_psa   <<- coef_psa[1:18]
  dgf_coef_psa    <<- coef_psa[19:36]
  mort30_coef_psa <<- coef_psa[37:54]
}

## Function to clean data between SIR runs
clean_SIR <- function() {
    sim_df <<- sim_df %>%
    dplyr::select(-c(constant, rec_cold_isch_tm)) %>%
    rename(months_to_gl = gl_time) %>%
    rename(months_to_gs_mort = gs_death_time) %>%
    rename(months_to_gl_mort = gl_death_time) %>%
    mutate(can_age_at_listing_c = rec_age_at_tx_c - (months_to_event/12)) %>%
    mutate(baseline_yrs_dial = years_dial - (months_to_event/12)) %>%
    mutate(list_year_c = tx_year_c - (months_to_event/12)) %>%
    mutate(start_time = ifelse(can_age_at_listing_c >=65, 0, 
                               abs(can_age_at_listing_c)-0.5)) %>%
    mutate(months_to_event = ifelse(months_to_event < 0, 0, months_to_event)) %>%
    mutate(ddtx = ifelse(ddtx==1 & months_to_event>240 & sim=="Simulated",
                         0, ddtx)) %>%
    mutate(ldtx = ifelse(ldtx==1 & months_to_event>240 & sim=="Simulated",
                         0, ldtx)) %>%
    mutate(mort = ifelse(mort==1 & months_to_event>240 & sim=="Simulated",
                         0, mort)) %>%
    mutate(remove = ifelse(remove==1 & months_to_event>240 & sim=="Simulated",
                           0, remove)) %>%
    mutate(months_to_event = ifelse(months_to_event>240 & sim=="Simulated",
                                    240, months_to_event)) %>%
    mutate(graft_loss = ifelse(graft_loss==1 & months_to_gl>120,
                               0, graft_loss)) %>%
    mutate(gs_death = ifelse(gs_death==1 & months_to_gs_mort>120,
                             0, gs_death)) %>%
    mutate(gl_death = ifelse(gl_death==1 & months_to_gl_mort>120,
                             0, gl_death)) %>%
    mutate(months_to_gl = ifelse(months_to_gl>120,
                                 120, months_to_gl)) %>%
    mutate(months_to_gs_mort = ifelse(months_to_gs_mort>120,
                                      120, months_to_gs_mort)) %>%
    mutate(months_to_gl_mort = ifelse(months_to_gl_mort>120,
                                      120, months_to_gl_mort)) %>%
    filter(!(ddtx==1 & rec_age_at_tx_c < 0))
}


##### Run Simulation ######
sim_time <- proc.time()  # start timer
list.simulation(n=100000, seed = 11111)
run.simulation(seed = 11111)
proc.time() - sim_time

##### Run Simulation on best SIR fit ######
best_list <- as.integer(which.max(score_list))
best_post <- as.integer(which.max(score_post))
best_mlog <- as.integer(which.max(score_mlogit))
best_glmort <- as.integer(which.max(score_gl_mort))

sim_time <- proc.time()  # start timer
list.simulation(n=100000, seed = 11111,
                ddtx_coef = ddtx_coef_mat[,best_list], ddtx_p = ddtx_shape_v[best_list],
                ldtx_coef = ldtx_coef_mat[,best_list], ldtx_gamma = ldtx_shape_v[best_list],
                mort_coef = mort_coef_mat[,best_list], mort_shape = mort_shape_v[best_list],
                remove_coef = remove_coef_mat[,best_list], remove_p = remove_shape_v[best_list])
run.simulation(seed = 11111, gl30_coef = mlogit_coef_mat[1:18, best_mlog], 
               dgf_coef = mlogit_coef_mat[19:36, best_mlog], 
               mort30_coef = mlogit_coef_mat[37:54, best_mlog], 
               gl_coef = gl_coef_mat[,best_post], gl_shape = gl_shape_v[best_post], 
               gs_mort_coef = gs_mort_coef_mat[,best_post], 
               gs_mort_shape = gs_mort_shape_v[best_post], 
               dial_mort_coef = dial_mort_coef_mat[,best_glmort], 
               dial_mort_shape = dial_mort_shape_v[best_glmort], 
               gl_mort_coef = gl_mort_coef_mat[,best_glmort])
proc.time() - sim_time

list.simulation(n=100000, seed = 11111,
                ddtx_coef = ddtx_coef_mat[,as.integer(best_list)])
##### Generate simulated data #####
data <- as.data.frame(data)
test_sim_df <- data %>% 
  select(-c(event, start_time, months_to_event, event_cd, dgf, graft_loss, 
            months_to_gl, gs_death, months_to_gs_mort, gl_death, 
            months_to_gl_mort, kdpi, tx_outcome, no_miss))
sim_dataset(n=100000)

##### Run Wait List Simulation 5 Times ####
y  <- NULL
list.simulation(n=100000, seed = 12345)
tmp <- list_sim_df
tmp$iter <- 1
y <- rbind(y, tmp)

list.simulation(n=100000, seed = 45643)
tmp <- list_sim_df
tmp$iter <- 2
y <- rbind(y, tmp)

list.simulation(n=100000, seed = 94562)
tmp <- list_sim_df
tmp$iter <- 3
y <- rbind(y, tmp)

list.simulation(n=100000, seed = 13865)
tmp <- list_sim_df
tmp$iter <- 4
y <- rbind(y, tmp)

list.simulation(n=100000, seed = 78946)
tmp <- list_sim_df
tmp$iter <- 5
y <- rbind(y, tmp)



