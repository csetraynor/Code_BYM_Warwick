#----Load libraries---#

library(purrr)
library(httr)
library(readr)
library(survival)
library(rstan)
library(spBayesSurv)
library(pracma)
library(assertthat)
library(cgdsr)
suppressMessages(library(dplyr))
library(ggplot2)
require(ggfortify)
theme_set(theme_bw())
library(VIM)
library(scales)
library(caret)
library(ipred)
library(e1071)
library(MASS)
library(lattice)
library(impute)
library(qdap)
library(doMC)
registerDoMC(cores = 5)
###############################################
#Data obtantion
#get data from with MSKCC package 

mycgds = CGDS("http://www.cbioportal.org/public-portal/")

glioblastome_2013_id_sutdy = getCancerStudies(mycgds)[55,1]
glioblastome_2013_case_list = getCaseLists(mycgds, glioblastome_2013_id_sutdy)[2,1]
glioblastome_2013_clinical_data <-  getClinicalData(mycgds, glioblastome_2013_case_list)

glio_clin_dat <- tbl_df(glioblastome_2013_clinical_data %>% tibble::rownames_to_column("sample_id")) 

####################################################################
#Data Cleaning and Preprocessing

#convert to lower case
names(glio_clin_dat) <- tolower(names(glio_clin_dat)) 

#convert missig values
convert_blank_to_na <- function(x) {
  if(!purrr::is_character(x)){
    warning('input variate not character - return original')
    return(x)
  } else {
    ifelse(x == '', NA, x)
  }
}
glio_clin_dat <- glio_clin_dat %>%
  dplyr::mutate_all(funs(convert_blank_to_na))

glio_clin_dat %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

#filter unknown or negative survival times (os_monts < 0)

glio_clin_dat %>%
  filter(is.na(os_status) | os_status != '') %>%
  filter(os_months <= 0 | is.na(os_months)) %>%
  dplyr::select(os_status, os_months) %>%
  dplyr::glimpse()

#for now this observation will be remove from the analysis

glio_short_dat <- glio_clin_dat %>%
  filter(!is.na(os_status) & os_status != '') %>%
  filter(os_months > 0 & !is.na(os_months))

#Check 44 fewer obsrvations than original if working with dfs
assertthat::assert_that(nrow(glio_short_dat) == nrow(glio_clin_dat) - 44)
glio_clin_dat <- glio_short_dat
remove(glio_short_dat)
glio_clin_dat %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

#correlataion matrix
corrM <-   tbl_df(model.matrix(~x2012_methylation_class + 
                                 expression_subtype + g.cimp_methylation + 
                                 idh1_status +
                                 mgmt_status, data = glio_clin_dat))
corrM %>%
  mutate_all(funs(as.integer)) %>%
  sjPlot::sjp.corr(sort.corr = T, show.legend = T)

# table(glio_clin_dat$g.cimp_methylation, glio_clin_dat$idh1_status, useNA = "always")
# table(glio_clin_dat$g.cimp_methylation, glio_clin_dat$idh1_status,
#       glio_clin_dat$mgmt_status, useNA = "always" )
# table(glio_clin_dat$mgmt_status, useNA = "always" )
# table(glio_clin_dat$g.cimp_methylation, glio_clin_dat$idh1_status, useNA = "always")
# glio_clin_dat$idh1_status[is.na(glio_clin_dat$idh1_status) & glio_clin_dat$g.cimp_methylation == "G-CIMP"] <- "R132H"
# glio_clin_dat$idh1_status[is.na(glio_clin_dat$idh1_status) & glio_clin_dat$g.cimp_methylation == "non-G-CIMP"] <- "WT"

#Create dummy vars

Xdummies <- dummyVars(Surv(os_months, os_deceased) ~ age +
                        g.cimp_methylation + idh1_status +
                        mgmt_status, data =  glio_clin_dat %>%
                        mutate(os_deceased = (os_status == "DECEASED")))

X <- tbl_df(predict(Xdummies, newdata =  glio_clin_dat %>%
               mutate(os_deceased = (os_status == "DECEASED"))))

names(X)<- stringr::str_replace_all(names(X), "-", "")
names(X) <- tolower(names(X))

#impute using caret package method bagged trees. For each predictor in the data, a bagged tree is created using all of the other predictors in the training set. When a new sample has a missing predictor value, the bagged model is used to predict the value. 
preProc <- caret::preProcess(X[-1], method = c("bagImpute"))
X[,-1] <- predict(preProc, X[-1], na.action = na.pass)
X[-1] <- round(X[-1])
#standardise continuos covariates such as age
preProc <- caret::preProcess(X[1], method = c("center", "scale"))
X[1] <- predict(preProc, X[1])
histogram(X$age)

#Near Zero Variance Predictors

nzv <- caret::nearZeroVar(X, saveMetrics= TRUE)
nzv[nzv$nzv,]
#Prepare clinical covariates and  ensemble covariates
X <- X %>% mutate(g.noncimp_or_idh1_wt = (g.cimp_methylationnongcimp | idh1_statuswt)) %>% dplyr::select(age, g_noncimp_wt = g.noncimp_or_idh1_wt, mgmt_meth = mgmt_statusmethylated) 
nzv <- caret::nearZeroVar(X, saveMetrics= TRUE)
nzv[nzv$nzv,]

################################################################
## Create data partition K folds for Cross Validation
K = 10
folds <- caret::createFolds(glio_clin_dat$os_status,
                                           k = 10,
                                           list = TRUE)
str(folds)
data <- glio_clin_dat %>%
  dplyr::select(os_status, os_months, sample_id) %>%
  cbind(X)
split_up <- lapply(folds, function(ind, dat) dat[ind,], dat = data)
unlist(lapply(split_up, nrow))

for( i in 1:K){
  train_data <- do.call(rbind, split_up[c(-i)])
  test_data <- split_up[[i]]
  
  ##Run Stan
  stan_file <- "Weibull.stan"
  
  #Function gen_stan_data and gen_init from Weibull_clinical.R
  nChain <- 4
  wei_fullfit <- rstan::stan(stan_file,
                             data = gen_stan_data(train_data, '~age+ 
                                                  g_noncimp_wt+ 
                                                  mgmt_meth'),
                             cores = min(nChain, parallel::detectCores()),
                             seed = 7327,
                             chains = nChain,
                             iter = 2000,
                             init = gen_inits(M = 3),
                             control = list(adapt_delta = 0.99, max_treedepth = 10)
  )
  
  ######### Posterior predicitive checks ######
  
  pp_alpha <- rstan::extract(wei_fullfit,'alpha')$alpha
  pp_mu <- rstan::extract(wei_fullfit,'mu')$mu
  pp_beta <- rstan::extract(wei_fullfit, 'beta_bg')$beta_bg
  
  pp_beta <-  split(pp_beta, seq(nrow(pp_beta)))
  pp_alpha <-  split(pp_alpha, seq(nrow(pp_alpha)))
  pp_mu <-  split(pp_mu, seq(nrow(pp_mu)))
  
  X_test <- test_data %>%
    select(age, g_noncimp_wt, mgmt_meth)
  
  pp_newdata <- 
    purrr::pmap(list(pp_beta, pp_alpha, pp_mu),
                function(pp_beta, pp_alpha, pp_mu) 
                  {weibull_sim_data(alpha = pp_alpha,
                  mu = pp_mu,
                  n = n_distinct(test_data$sample_id),
                  beta = pp_beta,
                X = X_test)             } )
  
  pp_survdata <- 
    pp_newdata %>%
    map(~ mutate(., os_deceased = os_status == 'DECEASED')) %>%
    map(~ survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
    map(fortify)
  
  #Calculate Brier Score
  
  
  
}





