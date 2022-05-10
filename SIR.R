## Load packages
library(doParallel)
library(doRNG)
library(foreach)


## Set-up parameters for SIR
n <- 10000
seed_num <- 12345
in_path <- "//tsclient/Documents/Simulation Model/Model Inputs/"

## Vectors of survival times
surv_times_15 <- seq(from = 60, to = 180, by = 60)
surv_times_10 <- seq(from = 60, to = 120, by = 60)
surv_times_gl <- seq(from = 3, to = 18, by = 3)

## Initialize waitlist matrices and vectors
ddtx_coef_mat   <- matrix(data=NA, nrow=nrow(ddtx_coef_w), ncol=n)
ldtx_coef_mat   <- matrix(data=NA, nrow=nrow(ldtx_coef_g), ncol=n)
mort_coef_mat   <- matrix(data=NA, nrow=nrow(mort_coef_w), ncol=n)
remove_coef_mat <- matrix(data=NA, nrow=nrow(remove_coef_w), ncol=n)

ddtx_shape_v    <- rep(NA, n)
ldtx_shape_v    <- rep(NA, n)
mort_shape_v    <- rep(NA, n)
remove_shape_v  <- rep(NA, n)

surv_sim_mat    <- matrix(data=NA, nrow=60, ncol=n)

ddtx_surv_sim_mat    <- matrix(data=NA, nrow=length(surv_times_15), ncol=n)
ldtx_surv_sim_mat    <- matrix(data=NA, nrow=length(surv_times_15), ncol=n)
mort_surv_sim_mat    <- matrix(data=NA, nrow=length(surv_times_15), ncol=n)
remove_surv_sim_mat  <- matrix(data=NA, nrow=length(surv_times_15), ncol=n)

ddtx_surv_sim_mat_race    <- matrix(data=NA, nrow=length(surv_times_10)*4, ncol=n)
ldtx_surv_sim_mat_race    <- matrix(data=NA, nrow=length(surv_times_10)*4, ncol=n)
mort_surv_sim_mat_race    <- matrix(data=NA, nrow=length(surv_times_10)*4, ncol=n)
remove_surv_sim_mat_race  <- matrix(data=NA, nrow=length(surv_times_10)*4, ncol=n)

## Initialize post-transplant matrices and vectors
mlogit_coef_mat    <- matrix(data=NA, nrow=nrow(rbind(gl30_coef_ml_noreg, 
                                                   dgf_coef_ml_noreg, 
                                                   mort30_coef_ml_noreg)), ncol=n)
gs_mort_coef_mat   <- matrix(data=NA, nrow=nrow(gs_mort_coef_g_noreg), ncol=n)
gl_coef_mat        <- matrix(data=NA, nrow=nrow(gl_coef_g_noreg), ncol=n)
dial_mort_coef_mat <- matrix(data=NA, nrow=nrow(dial_mort_coef_w_noreg), ncol=n)
gl_mort_coef_mat   <- matrix(data=NA, nrow=nrow(gl_mort_coef_logit_noreg), ncol=n)

gs_surv_sim_mat      <- matrix(data=NA, nrow=length(surv_times_10), ncol=n)
gl_surv_sim_mat      <- matrix(data=NA, nrow=length(surv_times_10), ncol=n)
gl_mort_surv_sim_mat <- matrix(data=NA, nrow=length(surv_times_gl), ncol=n)

gs_surv_sim_mat_race      <- matrix(data=NA, nrow=length(surv_times_10)*4, ncol=n)
gl_surv_sim_mat_race      <- matrix(data=NA, nrow=length(surv_times_10)*4, ncol=n)
gl_mort_surv_sim_mat_race <- matrix(data=NA, nrow=length(surv_times_gl)*4, ncol=n)

mlogit_sim_mat <- matrix(data=NA, nrow=4, ncol=n)

gs_mort_shape_v   <- rep(NA, n)
gl_shape_v        <- rep(NA, n)
dial_mort_shape_v <- rep(NA, n)

for (i in 1:n) {
  set.seed(seed_num) 
  psa_list()
  psa_posttx()
  
  ## Store coefficients for waitlist functions
  ddtx_coef_mat[,i]   <- ddtx_coef_psa
  ldtx_coef_mat[,i]   <- ldtx_coef_psa
  mort_coef_mat[,i]   <- mort_coef_psa
  remove_coef_mat[,i] <- remove_coef_psa
  
  ddtx_shape_v[i]     <- ddtx_shape_psa
  ldtx_shape_v[i]     <- ldtx_shape_psa
  mort_shape_v[i]     <- mort_shape_psa
  remove_shape_v[i]   <- remove_shape_psa
  
  ## Store coefficients for post-transplant functions
  mlogit_coef_mat[,i]    <- mlogit_coef_psa
  gs_mort_coef_mat[,i]   <- gs_mort_coef_psa
  gl_coef_mat[,i]        <- gl_coef_psa
  dial_mort_coef_mat[,i] <- dial_mort_coef_psa
  gl_mort_coef_mat[,i]   <- gl_mort_coef_psa
  
  gs_mort_shape_v[i]   <- gs_mort_shape_psa
  gl_shape_v[i]        <- gl_shape_psa
  dial_mort_shape_v[i] <- dial_mort_shape_psa
  
  seed_num = seed_num + i
}

comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}

parallel::detectCores()
n.cores <- 16

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "PSOCK"
)

doParallel::registerDoParallel(cl = my.cluster)
doRNG::registerDoRNG(seed = seed_num)


out <- foreach(big_loop = 1:n, .combine='comb', .multicombine=TRUE,
               .packages = c("MASS", "copula", "psych", "dplyr", "survival")) %dopar% {
                 
                 # set.seed(seed_num) 
                 psa_list()
                 psa_posttx()
                 
                 ## Store coefficients for waitlist functions
                 ddtx_coef   <- t(ddtx_coef_psa)
                 ldtx_coef   <- t(ldtx_coef_psa)
                 mort_coef   <- t(mort_coef_psa)
                 remove_coef <- t(remove_coef_psa)

                 ddtx_shape     <- ddtx_shape_psa
                 ldtx_shape     <- ldtx_shape_psa
                 mort_shape     <- mort_shape_psa
                 remove_shape   <- remove_shape_psa

                 ## Store coefficients for post-transplant functions
                 mlogit_coef    <- t(mlogit_coef_psa)
                 gs_mort_coef   <- t(gs_mort_coef_psa)
                 gl_coef        <- t(gl_coef_psa)
                 dial_mort_coef <- t(dial_mort_coef_psa)
                 gl_mort_coef   <- t(gl_mort_coef_psa)

                 gs_mort_shape   <- gs_mort_shape_psa
                 gl_shape        <- gl_shape_psa
                 dial_mort_shape <- dial_mort_shape_psa
                 
  list.simulation(n=100000, seed = 11111,
                  ddtx_coef = ddtx_coef_psa, ddtx_p = ddtx_shape_psa,
                  ldtx_coef = ldtx_coef_psa, ldtx_gamma = ldtx_shape_psa,
                  mort_coef = mort_coef_psa, mort_shape = mort_shape_psa,
                  remove_coef = remove_coef_psa, remove_p = remove_shape_psa)
  run.simulation(seed = 11111, gl30_coef = mlogit_coef_psa[1:18], 
                 dgf_coef = mlogit_coef_psa[19:36], mort30_coef = mlogit_coef_psa[37:54], 
                 gl_coef = gl_coef_psa, gl_shape = gl_shape_psa, 
                 gs_mort_coef = gs_mort_coef_psa, gs_mort_shape = gs_mort_shape_psa, 
                 dial_mort_coef = dial_mort_coef_psa, dial_mort_shape = dial_mort_shape_psa, 
                 gl_mort_coef = gl_mort_coef_psa)
  clean_SIR()
  
  ### Full Cohort Waitlist Outcomes
  fit_ddtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ 1,
                         data = sim_df)  # Deceased Donor Tx
  fit_ldtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ 1,
                         data = sim_df)  # Living Donor Tx
  fit_mort    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ 1,
                         data = sim_df)  # Waitlist Mortality
  fit_remove  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ 1,
                         data = sim_df)  # Other Removal
  
  ### Waitlist Outcomes By Race/Ethnicity
  fit_ddtx_race    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ 1 + race,
                              data = sim_df)  # Deceased Donor Tx
  fit_ldtx_race    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ 1 + race,
                              data = sim_df)  # Living Donor Tx
  fit_mort_race    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ 1 + race,
                              data = sim_df)  # Waitlist Mortality
  fit_remove_race  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ 1 + race,
                              data = sim_df)  # Other Removal
  

  ### 30-Day Outcome Table
  mlogit_t <- prop.table(table(sim_df$tx_outcome))
  
  ### Store Waitlist Survival Times
  ddtx_surv   <- summary(fit_ddtx, times = surv_times_15, extend=TRUE)
  ldtx_surv   <- summary(fit_ldtx, times = surv_times_15, extend=TRUE)
  mort_surv   <- summary(fit_mort, times = surv_times_15, extend=TRUE)
  remove_surv <- summary(fit_remove, times = surv_times_15, extend=TRUE)

  ddtx_surv_race   <- summary(fit_ddtx_race, times = surv_times_10, extend=TRUE)
  ldtx_surv_race   <- summary(fit_ldtx_race, times = surv_times_10, extend=TRUE)
  mort_surv_race   <- summary(fit_mort_race, times = surv_times_10, extend=TRUE)
  remove_surv_race <- summary(fit_remove_race, times = surv_times_10, extend=TRUE)
  
  ddtx_surv   <- as.list(ddtx_surv[["surv"]])
  ldtx_surv   <- as.list(ldtx_surv[["surv"]])
  mort_surv   <- as.list(mort_surv[["surv"]])
  remove_surv <- as.list(remove_surv[["surv"]])
  
  ddtx_surv_race   <- as.list(ddtx_surv_race[["surv"]])
  ldtx_surv_race   <- as.list(ldtx_surv_race[["surv"]])
  mort_surv_race   <- as.list(mort_surv_race[["surv"]])
  remove_surv_race <- as.list(remove_surv_race[["surv"]])
  
  ddtx_surv_sim   <- as.numeric(ddtx_surv[1:3])
  ldtx_surv_sim   <- as.numeric(ldtx_surv[1:3])
  mort_surv_sim   <- as.numeric(mort_surv[1:3])
  remove_surv_sim <- as.numeric(remove_surv[1:3])
  
  ddtx_surv_race_sim   <- as.numeric(ddtx_surv_race[1:8])
  ldtx_surv_race_sim   <- as.numeric(ldtx_surv_race[1:8])
  mort_surv_race_sim   <- as.numeric(mort_surv_race[1:8])
  remove_surv_race_sim <- as.numeric(remove_surv_race[1:8])
  
  ### Store Post-Transplant Outcomes
  # Multinomial logit 30-day outcomes
  mlogit_sim <- mlogit_t

  ### Graft Loss
  fit_gl            <- survfit(Surv(time = months_to_gl, graft_loss) ~ 1,
                               data = sim_df)  # Graft loss
  fit_gl_race       <- survfit(Surv(time = months_to_gl, graft_loss) ~ 1 + race,
                               data = sim_df)  # Graft loss
  gl_surv           <- summary(fit_gl, times = surv_times_10, extend=TRUE)
  gl_surv_race      <- summary(fit_gl_race, times = surv_times_10, extend=TRUE)
  
  gl_surv           <- as.list(gl_surv[["surv"]])
  gl_surv_race      <- as.list(gl_surv_race[["surv"]])
  
  gl_surv_sim        <- as.numeric(gl_surv[1:2])
  gl_surv_race_sim   <- as.numeric(gl_surv_race[1:8])
  
  ### Death with Function
  fit_gs_death <- survfit(Surv(time = months_to_gs_mort, gs_death) ~ 1, 
                          data = sim_df)  # Death w/ function
  fit_gs_death_race <- survfit(Surv(time = months_to_gs_mort, gs_death) ~ 1 + race, 
                               data = sim_df)  # Death w/ function
  
  gs_mort_surv      <- summary(fit_gs_death, times = surv_times_10, extend=TRUE)
  gs_mort_surv_race <- summary(fit_gs_death_race, times = surv_times_10, extend=TRUE)
  
  gs_mort_surv      <- as.list(gs_mort_surv[["surv"]])
  gs_mort_surv_race <- as.list(gs_mort_surv_race[["surv"]])
  
  gs_surv_sim        <- as.numeric(gs_mort_surv[1:2])
  gs_surv_race_sim   <- as.numeric(gs_mort_surv_race[1:8])
  
  ### Death after graft loss
  if(all(is.na(sim_df$months_to_gl_mort))==FALSE) {
    fit_gl_death <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ 1, 
                            data = sim_df)  # Death after graft loss
    fit_gl_death_race <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ 1 + race, 
                                 data = sim_df)  # Death after graft loss
    gl_mort_surv      <- summary(fit_gl_death, times = surv_times_gl, extend=TRUE)
    gl_mort_surv_race <- summary(fit_gl_death_race, times = surv_times_gl, extend=TRUE)
    
    gl_mort_surv      <- as.list(gl_mort_surv[["surv"]])
    gl_mort_surv_race <- as.list(gl_mort_surv_race[["surv"]])
    
    if (length(gl_mort_surv)==6) {
      gl_mort_surv_sim   <- as.numeric(gl_mort_surv[1:6])
    } else {
      gl_mort_surv_sim   <- rep(NA, 6)
    }
    
    if (length(gl_mort_surv_race)==24) {
      gl_mort_surv_race_sim   <- as.numeric(gl_mort_surv_race[1:24])
    } else {
      gl_mort_surv_race_sim   <- rep(NA, 24)
    }
  } else {
    gl_mort_surv_sim   <- rep(NA, 6)
    gl_mort_surv_race_sim   <- rep(NA, 24)
  }

  full_ddtx_sim   <- c(ddtx_surv_sim, ddtx_surv_race_sim)
  full_ldtx_sim   <- c(ldtx_surv_sim, ldtx_surv_race_sim)
  full_mort_sim   <- c(mort_surv_sim, mort_surv_race_sim)
  full_remove_sim <- c(remove_surv_sim, remove_surv_race_sim)

  full_gs_mort_sim   <- c(gs_surv_sim, gs_surv_race_sim)
  full_gl_sim        <- c(gl_surv_sim, gl_surv_race_sim)
  full_gl_mort_sim   <- c(gl_mort_surv_sim, gl_mort_surv_race_sim)

  # seed_num = seed_num + i
  
  list(full_ddtx_sim, full_ldtx_sim, full_mort_sim, full_remove_sim, 
       full_gs_mort_sim, full_gl_sim, full_gl_mort_sim, mlogit_sim,
       ddtx_coef, ldtx_coef, mort_coef, remove_coef, 
       ddtx_shape, ldtx_shape, mort_shape, remove_shape,
       mlogit_coef, gs_mort_coef, gl_coef, dial_mort_coef, gl_mort_coef,
       gs_mort_shape, gl_shape, dial_mort_shape)
}

parallel::stopCluster(cl = my.cluster)

for (i in 1:n) {
  list.simulation(n=100000, seed = 11111,
                  ddtx_coef = ddtx_coef_mat[,i], ddtx_p = ddtx_shape_v[i],
                  ldtx_coef = ldtx_coef_mat[,i], ldtx_gamma = ldtx_shape_v[i],
                  mort_coef = mort_coef_mat[,i], mort_shape = mort_shape_v[i],
                  remove_coef = remove_coef_mat[,i], remove_p = remove_shape_v[i])
  run.simulation(seed = 11111, gl30_coef = mlogit_coef_mat[1:18,i], 
                 dgf_coef = mlogit_coef_mat[19:36,i], mort30_coef = mlogit_coef_mat[37:54,i], 
                 gl_coef = gl_coef_mat[,i], gl_shape = gl_shape_v[i], 
                 gs_mort_coef = gs_mort_coef_mat[,i], 
                 gs_mort_shape = gs_mort_shape_v[i], 
                 dial_mort_coef = dial_mort_coef_mat[,i], 
                 dial_mort_shape = dial_mort_shape_v[i], 
                 gl_mort_coef = gl_mort_coef_mat[,i])
  clean_SIR()
  
  ### Full Cohort Waitlist Outcomes
  fit_ddtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ 1,
                         data = sim_df)  # Deceased Donor Tx
  fit_ldtx    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ 1,
                         data = sim_df)  # Living Donor Tx
  fit_mort    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ 1,
                         data = sim_df)  # Waitlist Mortality
  fit_remove  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ 1,
                         data = sim_df)  # Other Removal
  ### Waitlist Outcomes By Race/Ethnicity
  fit_ddtx_race    <- survfit(Surv(time = start_time, time2 = months_to_event, ddtx) ~ 1 + race,
                              data = sim_df)  # Deceased Donor Tx
  fit_ldtx_race    <- survfit(Surv(time = start_time, time2 = months_to_event, ldtx) ~ 1 + race,
                              data = sim_df)  # Living Donor Tx
  fit_mort_race    <- survfit(Surv(time = start_time, time2 = months_to_event, mort) ~ 1 + race,
                              data = sim_df)  # Waitlist Mortality
  fit_remove_race  <- survfit(Surv(time = start_time, time2 = months_to_event, remove) ~ 1 + race,
                              data = sim_df)  # Other Removal
  
  ### Full Cohort Post-Transplant Outcomes
  fit_gl       <- survfit(Surv(time = months_to_gl, graft_loss) ~ 1,
                          data = sim_df)  # Graft loss
  fit_gs_death <- survfit(Surv(time = months_to_gs_mort, gs_death) ~ 1, 
                          data = sim_df)  # Death w/ function
  fit_gl_death <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ 1, 
                          data = sim_df)  # Death after graft loss
  ### Post-Transplant Outcomes by Race/Ethnicity
  fit_gl_race       <- survfit(Surv(time = months_to_gl, graft_loss) ~ 1 + race,
                               data = sim_df)  # Graft loss
  fit_gs_death_race <- survfit(Surv(time = months_to_gs_mort, gs_death) ~ 1 + race, 
                               data = sim_df)  # Death w/ function
  fit_gl_death_race <- survfit(Surv(time = months_to_gl_mort, gl_death) ~ 1 + race, 
                               data = sim_df)  # Death after graft loss
  ### 30-Day Outcome Table
  mlogit_t <- prop.table(table(sim_df$tx_outcome))
  
  ### Store Waitlist Survival Times
  ddtx_surv   <- summary(fit_ddtx, times = surv_times_15, extend=TRUE)
  ldtx_surv   <- summary(fit_ldtx, times = surv_times_15, extend=TRUE)
  mort_surv   <- summary(fit_mort, times = surv_times_15, extend=TRUE)
  remove_surv <- summary(fit_remove, times = surv_times_15, extend=TRUE)
  
  ddtx_surv_race   <- summary(fit_ddtx_race, times = surv_times_10, extend=TRUE)
  ldtx_surv_race   <- summary(fit_ldtx_race, times = surv_times_10, extend=TRUE)
  mort_surv_race   <- summary(fit_mort_race, times = surv_times_10, extend=TRUE)
  remove_surv_race <- summary(fit_remove_race, times = surv_times_10, extend=TRUE)
  
  ddtx_surv   <- as.list(ddtx_surv[["surv"]])
  ldtx_surv   <- as.list(ldtx_surv[["surv"]])
  mort_surv   <- as.list(mort_surv[["surv"]])
  remove_surv <- as.list(remove_surv[["surv"]])
  
  ddtx_surv_race   <- as.list(ddtx_surv_race[["surv"]])
  ldtx_surv_race   <- as.list(ldtx_surv_race[["surv"]])
  mort_surv_race   <- as.list(mort_surv_race[["surv"]])
  remove_surv_race <- as.list(remove_surv_race[["surv"]])
  
  ddtx_surv_sim_mat[,i]   <- as.numeric(ddtx_surv[1:3])
  ldtx_surv_sim_mat[,i]   <- as.numeric(ldtx_surv[1:3])
  mort_surv_sim_mat[,i]   <- as.numeric(mort_surv[1:3])
  remove_surv_sim_mat[,i] <- as.numeric(remove_surv[1:3])
  
  ddtx_surv_sim_mat_race[,i]   <- as.numeric(ddtx_surv_race[1:8])
  ldtx_surv_sim_mat_race[,i]   <- as.numeric(ldtx_surv_race[1:8])
  mort_surv_sim_mat_race[,i]   <- as.numeric(mort_surv_race[1:8])
  remove_surv_sim_mat_race[,i] <- as.numeric(remove_surv_race[1:8])
  
  ### Store Post-Transplant Outcomes
  # Multinomial logit 30-day outcomes
  mlogit_sim_mat[,i] <- mlogit_t
  
  # Survival outcomes
  gs_mort_surv      <- summary(fit_gs_death, times = surv_times_10, extend=TRUE)
  gl_surv           <- summary(fit_gl, times = surv_times_10, extend=TRUE)
  gl_mort_surv      <- summary(fit_gl_death, times = surv_times_gl, extend=TRUE)
  
  gs_mort_surv_race <- summary(fit_gs_death_race, times = surv_times_10, extend=TRUE)
  gl_surv_race      <- summary(fit_gl_race, times = surv_times_10, extend=TRUE)
  gl_mort_surv_race <- summary(fit_gl_death_race, times = surv_times_gl, extend=TRUE)
  
  gs_mort_surv      <- as.list(gs_mort_surv[["surv"]])
  gl_surv           <- as.list(gl_surv[["surv"]])
  gl_mort_surv      <- as.list(gl_mort_surv[["surv"]])
  
  gs_mort_surv_race <- as.list(gs_mort_surv_race[["surv"]])
  gl_surv_race      <- as.list(gl_surv_race[["surv"]])
  gl_mort_surv_race <- as.list(gl_mort_surv_race[["surv"]])
  
  gs_surv_sim_mat[,i]        <- as.numeric(gs_mort_surv[1:2])
  gl_surv_sim_mat[,i]        <- as.numeric(gl_surv[1:2])
  
  gs_surv_sim_mat_race[,i]        <- as.numeric(gs_mort_surv_race[1:8])
  gl_surv_sim_mat_race[,i]        <- as.numeric(gl_surv_race[1:8])
  
  if (length(gl_mort_surv)==6) {
    gl_mort_surv_sim_mat[,i]   <- as.numeric(gl_mort_surv[1:6])
  }
  
  if (length(gl_mort_surv_race)==24) {
    gl_mort_surv_sim_mat_race[,i]   <- as.numeric(gl_mort_surv_race[1:24])
  }
}

full_ddtx_sim_mat   <- rbind(ddtx_surv_sim_mat, ddtx_surv_sim_mat_race)
full_ldtx_sim_mat   <- rbind(ldtx_surv_sim_mat, ldtx_surv_sim_mat_race)
full_mort_sim_mat   <- rbind(mort_surv_sim_mat, mort_surv_sim_mat_race)
full_remove_sim_mat <- rbind(remove_surv_sim_mat, remove_surv_sim_mat_race)
full_list_sim_mat   <- rbind(full_ddtx_sim_mat, full_ldtx_sim_mat,
                             full_mort_sim_mat, full_remove_sim_mat)

full_ddtx_obs_mat   <- rbind(ddtx_surv_obs_mat, ddtx_surv_obs_mat_race)
full_ldtx_obs_mat   <- rbind(ldtx_surv_obs_mat, ldtx_surv_obs_mat_race)
full_mort_obs_mat   <- rbind(mort_surv_obs_mat, mort_surv_obs_mat_race)
full_remove_obs_mat <- rbind(remove_surv_obs_mat, remove_surv_obs_mat_race)
full_list_obs_mat   <- rbind(full_ddtx_obs_mat, full_ldtx_obs_mat,
                             full_mort_obs_mat, full_remove_obs_mat)

full_gs_mort_sim_mat   <- rbind(gs_surv_sim_mat, gs_surv_sim_mat_race)
full_gl_sim_mat        <- rbind(gl_surv_sim_mat, gl_surv_sim_mat_race)
full_gl_mort_sim_mat   <- rbind(gl_mort_surv_sim_mat, gl_mort_surv_sim_mat_race)
full_post_sim_mat      <- rbind(full_gs_mort_sim_mat, full_gl_sim_mat)

full_gs_mort_obs_mat   <- rbind(gs_surv_obs_mat, gs_surv_obs_mat_race)
full_gl_obs_mat        <- rbind(gl_surv_obs_mat, gl_surv_obs_mat_race)
full_gl_mort_obs_mat   <- rbind(gl_mort_surv_obs_mat, gl_mort_surv_obs_mat_race)
full_post_obs_mat      <- rbind(full_gs_mort_obs_mat, full_gl_obs_mat)

sse_ddtx    <- colSums((full_ddtx_obs_mat - full_ddtx_sim_mat))^2
sse_ldtx    <- colSums((full_ldtx_obs_mat - full_ldtx_sim_mat))^2
sse_mort    <- colSums((full_mort_obs_mat - full_mort_sim_mat))^2
sse_remove  <- colSums((full_remove_obs_mat - full_remove_sim_mat))^2
sse_list    <- colSums((full_list_obs_mat - full_list_sim_mat))^2
sse_mlogit  <- colSums((mlogit_obs_mat - mlogit_sim_mat))^2

sse_gs       <- colSums((full_gs_mort_obs_mat - full_gs_mort_sim_mat))^2
sse_gl       <- colSums((full_gl_obs_mat - full_gl_sim_mat))^2
sse_gl_mort  <- colSums((full_gl_mort_obs_mat - full_gl_mort_sim_mat), na.rm = TRUE)^2
sse_post     <- colSums((full_post_obs_mat - full_post_sim_mat))^2

score_ddtx   <- 1/(sse_ddtx)
score_ddtx   <- score_ddtx/sum(score_ddtx) 
score_ldtx   <- 1/(sse_ldtx)
score_ldtx   <- score_ldtx/sum(score_ldtx) 
score_mort   <- 1/(sse_mort)
score_mort   <- score_mort/sum(score_mort)
score_remove <- 1/(sse_remove)
score_remove <- score_remove/sum(score_remove)
score_list   <- 1/(sse_list)
score_list   <- score_list/sum(score_list)
score_mlogit <- 1/(sse_mlogit)
score_mlogit <- score_mlogit/sum(score_mlogit)

score_gs      <- 1/(sse_gs)
score_gs      <- score_gs/sum(score_gs) 
score_gl      <- 1/(sse_gl)
score_gl      <- score_gl/sum(score_gl) 
score_gl_mort <- 1/(sse_gl_mort)
score_gl_mort[!is.finite(score_gl_mort)] <- NA
score_gl_mort <- score_gl_mort/sum(score_gl_mort, na.rm = TRUE)
score_post      <- 1/(sse_post)
score_post      <- score_post/sum(score_post)

na_gl_mort <- gl_mort_surv_obs_mat_race[1,]
na_gl_mort[!is.na(na_gl_mort)] <- 1
na_gl_mort[is.na(na_gl_mort)] <- 0

score_gl_mort <- score_gl_mort*na_gl_mort
score_gl_mort[is.na(score_gl_mort)] <- 0


### Create Dataframe for Analysis
survival_yr  <- "surv_yr_"
survival_mo  <- "surv_mo_"

year_surv_15     <- surv_times_15/12
year_surv_10     <- surv_times_10/12
year_surv_9      <- surv_times_9/12

surv_vars_15 <- lapply(year_surv_15, function(x) paste0(survival_yr, x))
surv_vars_10 <- lapply(year_surv_10, function(x) paste0(survival_yr, x))
surv_vars_9  <- lapply(year_surv_9, function(x) paste0(survival_yr, x))
surv_vars_gl <- lapply(surv_times_gl, function(x) paste0(survival_mo, x))

### Create DDTX Analysis
ddtx_sim <- as.data.frame(t(ddtx_surv_obs_mat))
colnames(ddtx_sim) <- surv_vars_15
ddtx_sim$sim <- 1
ddtx_sim$www <- score_list

ddtx_obs <- as.data.frame(t(ddtx_surv_sim_mat))
colnames(ddtx_obs) <- surv_vars_15
ddtx_obs$sim <- 0
ddtx_obs$www <- 0.005

ddtx_df <- rbind(ddtx_sim, ddtx_obs)


ddtx5 <- ggplot(ddtx_df, aes(x = sim, y=surv_yr_5, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0.2, 0.65) + scale_x_discrete(limits=c(0, 1))
ddtx10 <- ggplot(ddtx_df, aes(x = sim, y=surv_yr_10, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0.2, 0.65) + scale_x_discrete(limits=c(0, 1))
ddtx15 <- ggplot(ddtx_df, aes(x = sim, y=surv_yr_15, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0.2, 0.65) + scale_x_discrete(limits=c(0, 1))
plots_ddtx <- grid.arrange(ddtx5, ddtx10, ddtx15, nrow=1, top = "Deceased Donor Tx")

### Create LDTX Analysis
ldtx_sim <- as.data.frame(t(ldtx_surv_obs_mat))
colnames(ldtx_sim) <- surv_vars_15
ldtx_sim$sim <- 1
ldtx_sim$www <- score_list

ldtx_obs <- as.data.frame(t(ldtx_surv_sim_mat))
colnames(ldtx_obs) <- surv_vars_15
ldtx_obs$sim <- 0
ldtx_obs$www <- 0.005

ldtx_df <- rbind(ldtx_sim, ldtx_obs)


ldtx5 <- ggplot(ldtx_df, aes(x = sim, y=surv_yr_5, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0.6, 0.9) + scale_x_discrete(limits=c(0, 1))
ldtx10 <- ggplot(ldtx_df, aes(x = sim, y=surv_yr_10, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0.6, 0.9) + scale_x_discrete(limits=c(0, 1))
ldtx15 <- ggplot(ldtx_df, aes(x = sim, y=surv_yr_15, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0.6, 0.9) + scale_x_discrete(limits=c(0, 1))
plots_ldtx <- grid.arrange(ldtx5, ldtx10, ldtx15, nrow=1, top = "Living Donor Tx")

### Create Waitlist Mortality Analysis
mort_sim <- as.data.frame(t(mort_surv_obs_mat))
colnames(mort_sim) <- surv_vars_15
mort_sim$sim <- 1
mort_sim$www <- score_list

mort_obs <- as.data.frame(t(mort_surv_sim_mat))
colnames(mort_obs) <- surv_vars_15
mort_obs$sim <- 0
mort_obs$www <- 0.005

mort_df <- rbind(mort_sim, mort_obs)


mort5 <- ggplot(mort_df, aes(x = sim, y=surv_yr_5, weight=www, group=sim)) + 
  geom_boxplot() +  
  ylim(0.1, 0.8) + scale_x_discrete(limits=c(0, 1))
mort10 <- ggplot(mort_df, aes(x = sim, y=surv_yr_10, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0.1, 0.8) + scale_x_discrete(limits=c(0, 1))
mort15 <- ggplot(mort_df, aes(x = sim, y=surv_yr_15, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0.1, 0.8) + scale_x_discrete(limits=c(0, 1))
plots_mort <- grid.arrange(mort5, mort10, mort15, nrow=1, top = "Waitlist Mortality")

### Create Waitlist Remove Analysis
remove_sim <- as.data.frame(t(remove_surv_obs_mat))
colnames(remove_sim) <- surv_vars_15
remove_sim$sim <- 1
remove_sim$www <- score_list

remove_obs <- as.data.frame(t(remove_surv_sim_mat))
colnames(remove_obs) <- surv_vars_15
remove_obs$sim <- 0
remove_obs$www <- 0.005

remove_df <- rbind(remove_sim, remove_obs)


remove5 <- ggplot(remove_df, aes(x = sim, y=surv_yr_5, weight=www, group=sim)) + 
  geom_boxplot() +  
  ylim(0, 0.6) + scale_x_discrete(limits=c(0, 1))
remove10 <- ggplot(remove_df, aes(x = sim, y=surv_yr_10, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0, 0.6) + scale_x_discrete(limits=c(0, 1))
remove15 <- ggplot(remove_df, aes(x = sim, y=surv_yr_15, weight=www, group=sim)) + 
  geom_boxplot() + 
  ylim(0, 0.6) + scale_x_discrete(limits=c(0, 1))
plots_mort <- grid.arrange(remove5, remove10, remove15, nrow=1, top = "Waitlist Removal")

### Create Death with Function Analysis
gsmort_sim <- as.data.frame(t(gs_surv_obs_mat))
colnames(gsmort_sim) <- surv_vars_10
gsmort_sim$sim <- 1
gsmort_sim$www <- score_post

gsmort_obs <- as.data.frame(t(gs_surv_sim_mat))
colnames(gsmort_obs) <- surv_vars_10
gsmort_obs$sim <- 0
gsmort_obs$www <- 0.005

gsmort_df <- rbind(gsmort_sim, gsmort_obs)

gsmort5 <- ggplot(gsmort_df, aes(x = sim, y=surv_yr_5, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))
gsmort10 <- ggplot(gsmort_df, aes(x = sim, y=surv_yr_10, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))

plots_gsmort <- grid.arrange(gsmort5, gsmort10, nrow=1, top = "Death W/ Function")

### Create Graft Loss Analysis
gl_sim <- as.data.frame(t(gl_surv_obs_mat))
colnames(gl_sim) <- surv_vars_10
gl_sim$sim <- 1
gl_sim$www <- score_post

gl_obs <- as.data.frame(t(gl_surv_sim_mat))
colnames(gl_obs) <- surv_vars_10
gl_obs$sim <- 0
gl_obs$www <- 0.005

gl_df <- rbind(gl_sim, gl_obs)

gl5 <- ggplot(gl_df, aes(x = sim, y=surv_yr_5, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))
gl10 <- ggplot(gl_df, aes(x = sim, y=surv_yr_10, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))

plots_gl <- grid.arrange(gl5, gl10, nrow=1, top = "Graft Loss")

### Create Graft Loss Mortality Analysis
gl_mort_sim <- as.data.frame(t(gl_mort_surv_obs_mat))
colnames(gl_mort_sim) <- surv_vars_gl
gl_mort_sim$sim <- 1
gl_mort_sim$www <- score_gl_mort

gl_mort_obs <- as.data.frame(t(gl_mort_surv_sim_mat))
colnames(gl_mort_obs) <- surv_vars_gl
gl_mort_obs$sim <- 0
gl_mort_obs$www <- 0.005

gl_mort_df <- rbind(gl_mort_sim, gl_mort_obs)

glmort6 <- ggplot(gl_mort_df, aes(x = sim, y=surv_mo_6, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))
glmort12 <- ggplot(gl_mort_df, aes(x = sim, y=surv_mo_12, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))
glmort18 <- ggplot(gl_mort_df, aes(x = sim, y=surv_mo_18, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))

plots_gl <- grid.arrange(glmort6, glmort12, glmort18, nrow=1, top = "Graft Loss Mortality")

### Multinomial Logit Analysis ###
mlogit_sim <- as.data.frame(t(mlogit_sim_mat))
colnames(mlogit_sim) <- c("No_Comps", "Graft_Loss", "DGF", "Mortality")
mlogit_sim$sim <- 1
mlogit_sim$www <- score_mlogit

mlogit_obs <- as.data.frame(t(mlogit_obs_mat))
colnames(mlogit_obs) <- c("No_Comps", "Graft_Loss", "DGF", "Mortality")
mlogit_obs$sim <- 0
mlogit_obs$www <- 0.005

mlogit_df <- rbind(mlogit_sim, mlogit_obs)

mlogit1 <- ggplot(mlogit_df, aes(x = sim, y=No_Comps, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))
mlogit2 <- ggplot(mlogit_df, aes(x = sim, y=Graft_Loss, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))
mlogit3 <- ggplot(mlogit_df, aes(x = sim, y=DGF, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))
mlogit4 <- ggplot(mlogit_df, aes(x = sim, y=Mortality, weight=www, group=sim)) + 
  geom_boxplot() +  
  scale_x_discrete(limits=c(0, 1))

plots_mlogit <- grid.arrange(mlogit1, mlogit2, mlogit3, mlogit4, nrow=2, top = "30-Day Outcomes")
