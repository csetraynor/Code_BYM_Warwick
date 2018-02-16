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
###############################################
#Data obtantion
#get data from with MSKCC package 

mycgds = CGDS("http://www.cbioportal.org/public-portal/")

glioblastome_2013_id_sutdy = getCancerStudies(mycgds)[55,1]
glioblastome_2013_case_list = getCaseLists(mycgds, glioblastome_2013_id_sutdy)[2,1]
glioblastome_2013_clinical_data <-  getClinicalData(mycgds, glioblastome_2013_case_list)

glioblastome_2008_id_sutdy = getCancerStudies(mycgds)[56,1]
glioblastome_2008_case_list = getCaseLists(mycgds, glioblastome_2008_id_sutdy)[2,1]
glioblastome_2008_clinical_data <-  getClinicalData(mycgds, glioblastome_2008_case_list)

#inspect dataframes
glimpse(glioblastome_2013_clinical_data)
glimpse(glioblastome_2008_clinical_data)

# glioblastome_2013_clinical_data <- glioblastome_2013_clinical_data %>% tibble::rownames_to_column("sample_id"); glioblastome_2008_clinical_data <- glioblastome_2008_clinical_data %>% tibble::rownames_to_column("sample_id") 
# 
# glio_clin_dat <- glioblastome_2013_clinical_data %>%
#   filter(!(sample_id %in%  glioblastome_2008_clinical_data$sample_id ))


glio_clin_dat <- tbl_df(glioblastome_2013_clinical_data %>% tibble::rownames_to_column("sample_id")) 
####################################################################
#Data Cleaning

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

#inspect resulting dataframe
glimpse(glio_clin_dat)

######################################################################
#Data Exploration
#Considering overall survival#

glio_clin_dat %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

#filter unknown or negative survival times (os_monts < 0)

glio_clin_dat %>%
  filter(is.na(os_status) | os_status != '') %>%
  filter(os_months <= 0 | is.na(os_months)) %>%
  select(os_status, os_months) %>%
  dplyr::glimpse()

#for now this observation will be remove from the analysis

glio_short_dat <- glio_clin_dat %>%
  filter(!is.na(os_status) & os_status != '') %>%
  filter(os_months > 0 & !is.na(os_months))

#Check 44 fewer obsrvations than original if working with dfs
assertthat::assert_that(nrow(glio_short_dat) == nrow(glio_clin_dat) - 44)
glio_clin_dat <- glio_short_dat

glio_clin_dat %>%
  VIM::aggr(prop = FALSE, combined = TRUE, numbers = TRUE, sortVars = TRUE, sortCombs = TRUE)

########## Distribution of event times  ######################

glio_clin_dat %>%
  ggplot(aes(x = os_months,
             group = os_status,
             colour = os_status,
             fill = os_status)) +
  geom_density(alpha = 0.5)

#KM curve

mle.surv <- survfit(Surv(os_months, os_deceased) ~ 1,
                    data = glio_clin_dat %>%
                      mutate(os_deceased = (os_status == "DECEASED")))
ggplot2::autoplot(mle.surv, conf.int = F) +
  ggtitle('KM survival for GGM 2013 Cohort')



############ Parametric Survival Model #####################

gen_stan_data <- function(data, formula = as.formula(~ 1)) {
  
  if (!inherits(formula, 'formula'))
    formula <- as.formula(formula)
  
  observed_data <- data %>%
    dplyr::filter(os_status == 'DECEASED')
  
  censored_data <- data %>%
    dplyr::filter(os_status != 'DECEASED')
  
  Xobs_bg <- observed_data %>%
    model.matrix(formula, data = .)
  
  Xcen_bg <- censored_data %>% 
    model.matrix(formula, data = . )
  
  assertthat::assert_that(ncol(Xcen_bg) == ncol(Xobs_bg))
  M <- ncol(Xcen_bg)
  
  if (M > 1) {
    if ("(Intercept)" %in% colnames(Xobs_bg))
      Xobs_bg <- array(Xobs_bg[,-1], dim = c(nrow(observed_data), M - 1))
    if ("(Intercept)" %in% colnames(Xcen_bg))
      Xcen_bg <- array(Xcen_bg[,-1], dim = c(nrow(censored_data), M - 1))
    assertthat::assert_that(ncol(Xcen_bg) == ncol(Xobs_bg))
    M <- ncol(Xcen_bg)
  }
  
  stan_data <- list(
    Nobs = nrow(observed_data),
    Ncen = nrow(censored_data),
    yobs = observed_data$os_months,
    ycen = censored_data$os_months,
    M_bg = M,
    Xcen_bg = array(Xcen_bg, dim = c(nrow(censored_data), M)),
    Xobs_bg = array(Xobs_bg, dim = c(nrow(observed_data), M))
  )
}


##update initial values
load("historic_cohort.Rdata")
gen_inits <- function(M) {
  function()
    list(
      alpha_raw = 0.01*rnorm(1),
      mu = rnorm(1),
      tau_s_bg_raw = 0.1*abs(rnorm(1)),
      tau_bg_raw = array(abs(rnorm(M)), dim = c(M)),
      beta_bg_raw = array(rnorm(M), dim = c(M))
    )
}


##Prepare for fit Clinical Model

M <- cor(train[sapply(train, function(x) !is.character(x))])
#corrplot(M, method = "ellipse",order = "hclust")
##Updated

corrM <-   tbl_df(model.matrix(~x2012_methylation_class + 
                             expression_subtype + g.cimp_methylation + 
                             idh1_status +
               mgmt_status, data = glio_clin_dat))
corrM %>%
  mutate_all(funs(as.integer)) %>%
  sjPlot::sjp.corr(sort.corr = T, show.legend = T)

## impute idh1_status based on g.cimp
table(glio_clin_dat$g.cimp_methylation, glio_clin_dat$idh1_status, useNA = "always")

glio_clin_dat$idh1_status[is.na(glio_clin_dat$idh1_status) & glio_clin_dat$g.cimp_methylation == "G-CIMP"] <- "R132C"
glio_clin_dat$idh1_status[is.na(glio_clin_dat$idh1_status) & glio_clin_dat$g.cimp_methylation == "non-G-CIMP"] <- "WT"

table(glio_clin_dat$g.cimp_methylation, glio_clin_dat$idh1_status,
      glio_clin_dat$mgmt_status, useNA = "always" )


#Prepare clinical covariates, just for now we will drop the remaining missing values
glio_clin_dat <- glio_clin_dat %>%
  filter(!is.na( mgmt_status)) %>%
  mutate(age_centered = scale(age),  #scaling age because is a continuos covariate
         g.cimp_or_idh1_r = (I(g.cimp_methylation == "G-CIMP") | I(idh1_status != "WT")))


##Run Stan
stan_file <- "Weibull.stan"
#open stan file
if (interactive())
  file.edit(stanfile)
testfit <- rstan::stan(stan_file,
                        data = gen_stan_data(glio_clin_dat, '~ age_centered + 
                                                              g.cimp_or_idh1_r + 
                                                              I( mgmt_status=="METHYLATED") '),
                        init = gen_inits(M = 3),
                        iter = 4,
                        chains = 1
)
nChain <- 4
wei_fullfit <- rstan::stan(stan_file,
                        data = gen_stan_data(glio_clin_dat, '~ age_centered + 
                                                              g.cimp_or_idh1_r+ 
                                                              I(mgmt_status=="METHYLATED") '),
                        cores = min(nChain, parallel::detectCores()),
                        seed = 7327,
                        chains = nChain,
                        iter = 2000,
                        init = gen_inits(M = 3),
                        control = list(adapt_delta = 0.99, max_treedepth = 10)
)


####################### Checking convergence ###################

print(wei_fullfit) #(Check Rhat close to 1)

rstan::traceplot(wei_fullfit, c('lp__', 'beta_bg'), ncol = 2) #Review traceplot for log-posterior

rstan::traceplot(wei_fullfit, c('alpha','mu'), ncol = 1)    #Review traceplot for parameters of interest

if(interactive())
  shinystan::launch_shinystan(wei_fullfit)        #Launch shiny stan

######################### Posterior predicitive checks ###################################
#Simulate time to event data

weibull_sim_data <- function(alpha, mu, n, beta, X) {
  
  beta <- as.vector(as.numeric(beta))
  X <- array(matrix(as.numeric(X)), dim = c(n, length(beta)))
  
  #prognostic index
  hazard_ratio = X %*% beta 
  
  
  
  t = lapply(hazard_ratio, function(x) 
    rweibull(n = 1, shape = alpha, scale = exp(-(mu + hazard_ratio)/alpha)))
  t = do.call(rbind, t)
  
  
  data <- data.frame(surv_months = t,
                     censor_months = rexp(n = n, rate = 1/100),
                     stringsAsFactors = F
  ) %>%
    dplyr::mutate(os_status = ifelse(surv_months < censor_months,
                                     'DECEASED', 'LIVING'
    ),
    os_months = ifelse(surv_months < censor_months,
                       surv_months, censor_months
    )
    )
  
  return(data)
}

test_alpha <- 0.8
test_mu <- -4
test_n <- 100
test_X = matrix(c(rnorm(100), sample(c(0,1), 100, replace = TRUE), ncol=2))
test_beta = c(0.5, 1)
weibull_sim_data(alpha = test_alpha, mu = test_mu, n = test_n, beta = test_beta, X = test_X)



#Censoring is "arbitrarily" rexp() , censoring is assumed to be noninformative.

######## Simulating data for each posterior draw #
pp_alpha <- rstan::extract(wei_fullfit,'alpha')$alpha
pp_mu <- rstan::extract(wei_fullfit,'mu')$mu
pp_beta <- rstan::extract(wei_fullfit, 'beta_bg')$beta_bg


# create list
pp_beta <-  split(pp_beta, seq(nrow(pp_beta)))
pp_alpha <-  split(pp_alpha, seq(nrow(pp_alpha)))
pp_mu <-  split(pp_mu, seq(nrow(pp_mu)))

X <- model.matrix(~ age_centered + 
                 I(g.cimp_methylation=="G-CIMP")+ 
                 I(mgmt_status=="METHYLATED"), data = glio_clin_dat)
M <- ncol(X)
if ("(Intercept)" %in% colnames(X))
  X <- array(X[,-1], dim = c(nrow(glio_clin_dat), M - 1))

pp_newdata <- 
  purrr::pmap(list(pp_beta, pp_alpha, pp_mu),
              function(pp_beta, pp_alpha, pp_mu) {weibull_sim_data(alpha = pp_alpha,
                                                                   mu = pp_mu,
                                                              n = n_distinct(glio_clin_dat$sample_id),
                                                                   beta = pp_beta,
                                                                   X = X)
              } )

###### Plot time to event in the posterior draws compare to actual time in dataset
ggplot(pp_newdata %>%
         bind_rows() %>%
         mutate(type = 'posterior predicted values') %>%
         bind_rows(glio_clin_dat %>% mutate(type = 'actual data'))
       , aes(x = os_months, group = os_status, colour = os_status, fill = os_status))+ geom_density(alpha = 0.5) + facet_wrap(~type, ncol = 1)

#### summarise posterior predictive draws

## cumulative survival rate at each draw from the posterior

pp_survdata <- 
  pp_newdata %>%
  map(~ mutate(., os_deceased = os_status == 'DECEASED')) %>%
  map(~ survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
  map(fortify)

## summarise cum survival for each unit time (month), summarised at 95% confidence interval
pp_survdata_agg <- 
  pp_survdata %>%
  map(~mutate(., time_group = floor(time))) %>%
  bind_rows() %>%
  group_by(time_group) %>%
  summarize(surv_mean = mean(surv),
            surv_p50 = median(surv),
            surv_lower = quantile(surv, probs = 0.025),
            surv_upper = quantile(surv, probs = 0.975)) %>%
  ungroup()

## km
kmcurve_data <-   fortify(
  survfit(
    Surv(os_months, os_deceased) ~ 1,
    data = glio_clin_dat %>%
      mutate(os_deceased = os_status == 'DECEASED')
  )) %>%
  mutate(lower =  surv,
         upper = surv)

ggplot(pp_survdata_agg %>%
         mutate(type = 'posterior predicted values') %>%
         rename(surv = surv_p50, lower = surv_lower, upper = surv_upper, time = time_group)
       %>%
         bind_rows(kmcurve_data %>% mutate(type = 'actual data')),
       aes(x = time, group = type, linetype = type)) +
  geom_line(aes(y = surv, colour = type)) + 
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  xlim(c(0, 200))+
  ggtitle("Glioblastome cohort")





############# Save as a function  ######################

pp_predict_surv <- function(pp_alpha, pp_mu, n, pp_beta, X, 
                            level = 0.9,
                            plot = F, data = NULL,
                            sim_data_fun = weibull_sim_data) {
  pp_newdata <- 
    purrr::pmap(list(pp_beta, pp_alpha, pp_mu),
                function(pp_beta, pp_alpha, pp_mu) {sim_data_fun(alpha = pp_alpha,
                                                                     mu = pp_mu,
                                                                     n = n,
                                                                     beta = pp_beta,
                                                                     X = X)
                } )
  
  pp_survdata <-
    pp_newdata %>%
    purrr::map(~ dplyr::mutate(., os_deceased = os_status == 'DECEASED')) %>%
    purrr::map(~ survival::survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
    purrr::map(fortify)
  
  ## compute quantiles given level 
  lower_p <- 0 + ((1 - level)/2)
  upper_p <- 1 - ((1 - level)/2)
  
  pp_survdata_agg <- 
    pp_survdata %>%
    purrr::map(~ dplyr::mutate(.,
                               time_group = floor(time))) %>%
    dplyr::bind_rows() %>%
    dplyr::group_by(time_group) %>%
    dplyr::summarize(surv_mean = mean(surv)
                     , surv_p50 = median(surv)
                     , surv_lower = quantile(surv,
                                             probs = lower_p)
                     , surv_upper = quantile(surv,
                                             probs = upper_p)
    ) %>%
    dplyr::ungroup()
  
  if (plot == FALSE) {
    return(pp_survdata_agg)
  } 
  
  ggplot_data <- pp_survdata_agg %>%
    dplyr::mutate(type = 'posterior predicted values') %>%
    dplyr::rename(surv = surv_p50,
                  lower = surv_lower,
                  upper = surv_upper, time = time_group)
  
  if (!is.null(data)){
    ggplot_data <- 
      ggplot_data %>% 
      bind_rows(
        fortify(
          survival::survfit(
            Surv(os_months, os_deceased) ~ 1, 
            data = data %>% 
              dplyr::mutate(
                os_deceased = os_status == 'DECEASED')
          )) %>%
          dplyr::mutate(lower = surv,
                        upper = surv, type = 'actual data')
      )}
  
  pl <- ggplot(ggplot_data,
               aes(x = time, group = type, linetype = type)) + 
    geom_line(aes(y = surv, colour = type)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2)
  
  pl 
}



################# Posterior predictive check ###################


pp_alpha <- rstan::extract(wei_fullfit,'alpha')$alpha
pp_mu <- rstan::extract(wei_fullfit,'mu')$mu
pp_beta <- rstan::extract(wei_fullfit, 'beta_bg')$beta_bg

pp_beta <-  split(pp_beta, seq(nrow(pp_beta)))
pp_alpha <-  split(pp_alpha, seq(nrow(pp_alpha)))
pp_mu <-  split(pp_mu, seq(nrow(pp_mu)))

X <- model.matrix(~ age_centered + 
                    I(g.cimp_methylation=="G-CIMP")+ 
                    I(mgmt_status=="METHYLATED"), data = glio_clin_dat)
M <- ncol(X)
if ("(Intercept)" %in% colnames(X))
  X <- array(X[,-1], dim = c(nrow(glio_clin_dat), M - 1))

pl <- pp_predict_surv(pp_alpha = pp_alpha,
                      pp_mu = pp_mu,
                      pp_beta = pp_beta,
                      X = X, 
                      n = nrow(glio_clin_dat),
                      data = glio_clin_dat,
                      plot = T
) 
pl + 
  xlim(NA, 250) +
  ggtitle('Posterior predictive checks for NULL weibull model\nfit to GBC 2008 historical cohort; showing 90% CI')
