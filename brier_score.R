#Brier Score

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
              function(pp_beta, pp_alpha, pp_mu) 
                {weibull_sim_data(alpha = pp_alpha,
                    mu = pp_mu,
                n = n_distinct(glio_clin_dat$sample_id),
                    beta = pp_beta,
                   X = X)
              } )


  
pp_smod <- lapply(pp_newdata %>%
       map(~ mutate(., os_deceased = os_status == 'DECEASED')), function(x){
         with(x, Surv(os_months, os_deceased))
       })
     
pp_survdata <- lapply(pp_smod, function(x){
  survival::survfit(x ~ 1)
}) 
  
# integrated Brier score up to max(DLBCL$time)
sbrier(pp_smod, pp_survdata)

pp_brier <- purrr::map2(.x = pp_smod,
                          .y = pp_survdata,
                          .f = ~ipred::sbrier(obj = .x,
                                       pred = .y)) %>% unlist()

ggplot2::ggplot(pp_brier)



# integrated Brier score up to time=50
sbrier(smod, KM, btime=c(0, 50))

# Brier score for time=50
sbrier(smod, KM, btime=50)


pp_smod %>%
  map(~ survfit(Surv(os_months, os_deceased) ~ 1, data = .)) %>%
  map(fortify)




