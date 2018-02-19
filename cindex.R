ApparrentCindex  <- pec::cindex(list("Cox X1"=cox1,
                                     "Cox X2"=cox2,
                                     "Cox X1+X2"=cox12,
                                     "RSF"=rsf1),
                                formula=Surv(time,status)~X1+X2,
                                data=dat,
                                eval.times=seq(5,500,50))

library(rms)

Cindex <- purrr::map(.x = pred_KM,
                     .f = ~pec::cindex(object = .x,
                              formula =  Surv(os_months, os_deceased) ~ )


print(ApparrentCindex)
plot(ApparrentCindex)