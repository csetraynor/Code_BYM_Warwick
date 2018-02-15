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

glioblastome_2013_clinical_data <- glioblastome_2013_clinical_data %>% tibble::rownames_to_column("sample_id"); glioblastome_2008_clinical_data <- glioblastome_2008_clinical_data %>% tibble::rownames_to_column("sample_id") 

glio_clin_dat <- glioblastome_2013_clinical_data %>%
  filter(!(sample_id %in%  glioblastome_2008_clinical_data$sample_id ))



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
    M = M,
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
      mu = rnorm(n = 1, mean = hist_cohort$mean_mu, sd = hist_cohort$sd_mu),
      tau_s_bg_raw = 0.1*abs(rnorm(1)),
      tau_bg_raw = array(abs(rnorm(M)), dim = c(M)),
      beta_bg_raw = array(rnorm(M), dim = c(M))
    )
}


##Prepare for fit Clinical Model

standardise <- function(x) {
  x <- (x-mean(x)) / sd(x) 
return(x)
}

glio_clin_dat <- glio_clin_dat %>%
  filter(!is.na(g.cimp_methylation) & !is.na( mgmt_status)) %>%
  mutate(age_centered = standardise(age))


##Run Stan
stan_file <- "Weibull.stan"
#open stan file
if (interactive())
  file.edit(stanfile)
testfit <- rstan::stan(stan_file,
                        data = gen_stan_data(glio_clin_dat, '~ age_centered + 
                                                              I(g.cimp_methylation=="G-CIMP")+ 
                                                              I( mgmt_status=="METHYLATED") '),
                        init = gen_inits(M = 3),
                        iter = 4,
                        chains = 1
)
nChain <- 4
wei_fullfit <- rstan::stan(stan_file,
                        data = gen_stan_data(glio_clin_dat, '~ age_centered + 
                                                              I(g.cimp_methylation=="G-CIMP")+ 
                                             I(mgmt_status=="METHYLATED") '),
                        cores = min(nChain, parallel::detectCores()),
                        seed = 7327,
                        chains = nChain,
                        iter = 2000,
                        init = gen_inits(M = 3),
                        control = list(adapt_delta = 0.99, max_treedepth = 10)
)



