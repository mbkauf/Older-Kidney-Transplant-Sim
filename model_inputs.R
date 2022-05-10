`" 
*-------------------------------------------------------------------------------
* Simulation Inputs
* model_inputs.R
* Date Created: 01/27/2022
* Last Updated: 04/25/2022
* Matt Kaufmann, Stanford University
* Note: Modified to run on remote server
*-------------------------------------------------------------------------------
"`

library(readxl)
library(haven)

##### file paths #####
in_path <- "//tsclient/Documents/Simulation Model/Model Inputs/"

##### Data set #####
sample_df       <- as.matrix(read.csv(file = paste0(in_path, "sample_df.csv"), 
                                 sep=","))

##### Equation Coefficients #####
### Wait List Equation Variables ###
# Deceased Donor Tx (AFT Weibull)
ddtx_coef_w   <- as.matrix(read_excel(paste0(in_path, "list_results_input.xlsx"),
                                      range = "B5:B30", col_names = FALSE))
ddtx_p_w      <- read_excel(paste0(in_path, "list_results_input.xlsx"), 
                            range = "B31", col_names = FALSE)[[1]]
ddtx_cov <- as.matrix(read_excel(paste0(in_path, "ddtx_results.xlsx"), 
                                 sheet = "cov", col_names = FALSE))

# Living Donor Tx (PH Gompertz)
ldtx_coef_g   <- as.matrix(read_excel(paste0(in_path, "list_results_input.xlsx"),
                                      range = "C5:C30", col_names = FALSE))
ldtx_gamma_g  <- read_excel(paste0(in_path, "list_results_input.xlsx"), 
                            range = "C31", col_names = FALSE)[[1]]
ldtx_cov <- as.matrix(read_excel(paste0(in_path, "ldtx_results.xlsx"), 
                                 sheet = "cov", col_names = FALSE))

# Wait list mortality (AFT Weibull)
mort_coef_w   <- as.matrix(read_excel(paste0(in_path, "list_results_input.xlsx"),
                                      range = "D5:D30", col_names = FALSE))
mort_p_w  <- read_excel(paste0(in_path, "list_results_input.xlsx"), 
                        range = "D31", col_names = FALSE)[[1]]
mort_cov <- as.matrix(read_excel(paste0(in_path, "mort_results.xlsx"), 
                                 sheet = "cov", col_names = FALSE))

# Other wait list removals (AFT Weibull)
remove_coef_w   <- as.matrix(read_excel(paste0(in_path, "list_results_input.xlsx"),
                                        range = "E5:E30", col_names = FALSE))
remove_p_w      <- read_excel(paste0(in_path, "list_results_input.xlsx"), 
                              range = "E31", col_names = FALSE)[[1]]
remove_cov <- as.matrix(read_excel(paste0(in_path, "remove_results.xlsx"), 
                                   sheet = "cov", col_names = FALSE))

# 30-Day Outcomes (w/o region)
# Base Outcome: Graft Success
gl30_coef_ml_noreg   <- as.matrix(read_excel(paste0(in_path, "tx_results_input_no_reg.xlsx"),
                                             range = "C5:C22", col_names = FALSE))
dgf_coef_ml_noreg    <- as.matrix(read_excel(paste0(in_path, "tx_results_input_no_reg.xlsx"),
                                             range = "D5:D22", col_names = FALSE))
mort30_coef_ml_noreg <- as.matrix(read_excel(paste0(in_path, "tx_results_input_no_reg.xlsx"),
                                             range = "E5:E22", col_names = FALSE))
mlogit_cov <- as.matrix(read_excel(paste0(in_path, "mlogit_results.xlsx"),
                                   sheet = "cov", col_names = FALSE))

# Post-Tx: Graft Loss (Gompertz) (w/o region)
gl_coef_g_noreg <- as.matrix(read_excel(paste0(in_path, "gl_results_input_no_reg.xlsx"),
                                        range = "B5:B23", col_names = FALSE))
gl_gamma_g_noreg <- read_excel(paste0(in_path, "gl_results_input_no_reg.xlsx"), 
                               range = "B24", col_names = FALSE)[[1]]
gl_cov <- as.matrix(read_excel(paste0(in_path, "gl_results.xlsx"), 
                               sheet = "cov", col_names = FALSE))

# Post-Tx: Graft Loss (Weibull) (w/o region)
gl_coef_w_noreg <- as.matrix(read_excel(paste0(in_path, "gl_results_input_no_reg.xlsx"),
                                        range = "C5:C23", col_names = FALSE))
gl_p_w_noreg <- read_excel(paste0(in_path, "gl_results_input_no_reg.xlsx"), 
                               range = "C24", col_names = FALSE)[[1]]

# Post-Tx: Death with Function (Gompertz) (w/o region)
gs_mort_coef_g_noreg <- as.matrix(read_excel(paste0(in_path, "gs_mort_results_input_no_reg.xlsx"),
                                             range = "B5:B22", col_names = FALSE))
gs_mort_gamma_g_noreg <- read_excel(paste0(in_path, "gs_mort_results_input_no_reg.xlsx"), 
                                    range = "B23", col_names = FALSE)[[1]]
gs_mort_cov <- as.matrix(read_excel(paste0(in_path, "gs_mort_results.xlsx"),
                                    sheet = "cov", col_names = FALSE))

# w/o region 
# Post-Tx: Death after graft loss (put months2gl variable last) (Weibull)
dial_mort_coef_w_noreg <- as.matrix(read_excel(paste0(in_path, "gl_mort_results_input_no_reg.xlsx"),
                                               range = "B5:B24", col_names = FALSE))
dial_mort_lnp_w_noreg <- read_excel(paste0(in_path, "gl_mort_results_input_no_reg.xlsx"), 
                                    range = "B25", col_names = FALSE)[[1]]
dial_mort_cov <- as.matrix(read_excel(paste0(in_path, "gl_mort_results.xlsx"),
                                      sheet = "cov", col_names = FALSE))

### Same day graft loss - mortality 
gl_mort_coef_logit_noreg <- as.matrix(read_excel(paste0(in_path, 
                                                        "gl_mort_results_logit_input_no_reg.xlsx"),
                                                 range = "B5:B24", col_names = FALSE))
gl_mort_cov <- as.matrix(read_excel(paste0(in_path, "gl_mort_results_logit.xlsx"),
                                    sheet = "cov", col_names = FALSE))

##### Policy Inputs #####
tx_bump <- ((1+0.50)^(1/12))-1  # Monthly Rate (50% Annual)

##### Costs #####
c_dialysis = 3639  # monthly dialysis cost
c_ddtx = 106675 # one-time deceased donor transplant cost
c_ldtx = 94000 # one-time deceased donor transplant cost
c_dgf = 3639 # cost of experiencing delayed graft function
c_graftloss = 94508 # one-time graft loss cost
c_gsdeath = 24228  # one-time death with function cost
c_wldeath = 64494  # one-time death on wait list
c_txcare = 2000 # monthly post-tx care costs

##### QALYs #####
u_dialysis = 0.71  # utility of dialysis
u_graft_func = 0.82  # utility with functioning graft
u_post_gl = 0.71  # utility post-graft loss
u_gl = 0.20 # disutility of graft loss event
u_death = 0  # utility after death

##### Discount Rate #####
r <- ((1+0.03)^(1/12))-1  # Monthly Discount Rate (3% Annual)


myfun=function(x){
  mysum <- 0
  for (i in 1:length(x)){
    mysum = mysum + x[i]/(1+r)^i
  }
  return(mysum)
}