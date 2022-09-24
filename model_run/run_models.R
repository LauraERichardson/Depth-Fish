library(tidyverse)
library(brms)

inv_logit <- function(z) 1/(1 + exp(-z))

# Finalized on 26th July 2022

# By LE Richardson, AJ Delargy and P Neubauer 

load('../intermed_data/Depth_study_fish_data.RData')

####
####---- Priors ----#####
####

# MacNeil et al. 2015 Nature: Resident fish biomass in absense of fishing averages 1,013 (963, 1469) kg ha-1 
# (posterior median (95% highest posterior density intervals)). So in g-m2: 101.3 (96.3, 146.9)
# Prior mean set to log of the expected intercept, with prior sd inflated to reflect uncertainty about translating ranges from MacNeil et al. 2015
prior_fish <- set_prior('normal(4.6,1)', class='Intercept')

#Trophic group proportions

#PISC: MacNeil et al. Extended Data Figure 6: Average reef fish functional group across a biomass gradient
#Prop: 0.2 mean, reflect uncertainty in proportions by drawing from beta distribution * totfish to derive prior
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,10,40)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,10,40)))
#2.95,0.35
# need to adjust for prior of absence in hurdle
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,10,40)*(1-inv_logit(rlogis(10000,-2,0.5)))))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,10,40)*(1-inv_logit(rlogis(10000,-2,0.5)))))
#2.77,0.4

# Prior mean set to log of the expected intercept, with prior sd inflated to reflect uncertainty about translating ranges from MacNeil et al. 2015
prior_PISC <- set_prior('normal(2.78,1)', class='Intercept')


#PLANK: MacNeil et al. Extended Data Figure 6: Average reef fish functional group across a biomass gradient
#Prop: 0.3 mean
# adjusted for hurdle prior
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,15,35)*(1-inv_logit(rlogis(10000,-2,0.5)))))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,15,35)*(1-inv_logit(rlogis(10000,-2,0.5)))))
#3.19,0.35
# Prior mean set to log of the expected intercept, with prior sd inflated to reflect uncertainty about translating ranges from MacNeil et al. 2015
prior_PLANK <- set_prior('normal(3.19,1)', class='Intercept')

#PRIM: 
#(a+b+c+d in MacNeil et al. Extended Data Figure 6: Average reef fish functional group across a biomass gradient)
#Prop: 0.8 mean
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,40,10)))
sd(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,40,10)))
#4.37, 0.22
# Prior mean set to log of the expected intercept, with prior sd inflated to reflect uncertainty about translating ranges from MacNeil et al. 2015
prior_PRIM <- set_prior('normal(4.37,1)', class='Intercept')

#SEC:
#Mean proportion of secondary consumers in study data:
PropSEC<-fish$SECONDARY/fish$TotFish
summary(PropSEC)
#Min.      1st Qu.   Median    Mean      3rd Qu.   Max. 
#0.0003287 0.1145833 0.2025365 0.2730447 0.3523987 1.0000000
#0.27
hist(rbeta(10000,27,73))
mean(log(rlnorm(10000, 4.6, 0.21)*rbeta(10000,27,73)))
#3.27793
prior_SEC <- set_prior('normal(3.27793,1)', class='Intercept')

def_priors = c(set_prior('normal(0,2)', class='b'), #regressions pars
               set_prior('cauchy(0,5)', class='sd'), #group effect sd
               set_prior('cauchy(0,10)', class='sds'), #wiggliness
               set_prior('cauchy(0,1)', class='sd', coef = 'DEPTH_c', group = 'ISLAND')) #slope of depth within island effect

####---- Model formula ----#####

form <- ~ DEPTH_c +
  POP_STATUS +
  SITE_SLOPE_400m_c +
  #DEPTH_c:POP_STATUS +
  s(DEPTH_c, by=POP_STATUS) +
  s(DEPTH_c, SITE_SLOPE_400m_c) +
  (1|OBS_YEAR) +
  (1|OBS_YEAR:ECOREGION) + 
  (1|OBS_YEAR:ISLAND) + 
  (1|ECOREGION)  + 
  (1|SITE) +
  (1|DIVER) +
  (1+DEPTH_c||ISLAND)


####---- Model settings ----#####

input_frame <- data.frame(
  SHORT = c('fish', ## shorts for naming
            'PISC',
            'PLANK',
            'PRIM',
            'SEC'),
  RESP  = c( "TotFish", ## model responses
            'PISCIVORE',
            'PLANKTIVORE',
            'PRIMARY',
            'SECONDARY'),
  HURDLE = c(F, T, T, F, F) ## Run hurdle?
)

run_depth_model <- function(input_frame, 
                            iter = 2500, 
                            warmup = 500){
  
  attach(input_frame, warn.conflicts = F) # note, input frame is assumed to be dim=c(1,p)
  
  pos_form <- formula(paste(c(RESP, form), collapse = " ")) # Formula for models with positive data only
  bin_form <- formula(paste(c('hu'       , form), collapse = " "))
  
  hu_form <- bf(pos_form, bin_form) # Formula for hurdle models with zeroes
  
  priors <- c(def_priors,
              get(paste0("prior_",SHORT)))
  
  dir.create(file.path("../model_output/",SHORT),showWarnings = F,recursive = T)
  
   if(HURDLE) {
    assign(SHORT, fish)
    priors <- c(priors,
                set_prior('logistic(-2,0.5)', class='Intercept', dpar = 'hu'),
                set_prior('normal(0,2)', class='b', dpar = 'hu'),
                set_prior('cauchy(0,5)', class='sd', dpar = 'hu'),
                set_prior('cauchy(0,10)', class='sds', dpar = 'hu'),
                set_prior('cauchy(0,1)', class='sd', coef = 'DEPTH_c', group = 'ISLAND', dpar = 'hu'))
  } else {
    assign(SHORT, filter(fish,!!sym(RESP)>0) %>% droplevels)
  }
  
  Gamma.brms.mod <<- brm(formula = if(HURDLE) hu_form else pos_form,
                         data = if(HURDLE) fish else get(SHORT), 
                         family = if(HURDLE) 'hurdle_gamma' else 'gamma',
                         prior = priors,
                         chains = 4,
                         thin = 4,
                         cores=4,
                         iter = get('iter'),
                         backend = 'cmdstanr',
                         file_refit = "always",
                         warmup = get('warmup'),
                         control = list(adapt_delta=0.99),
                         file = paste0("../model_output/",SHORT,"/",SHORT,"_mod"), 
                         refresh=10)
  
  png(filename = paste0("../model_output/",SHORT,"/model_diagnostics%d.png"), width = 900, height = 600)
  
  print({
    mcmc_plot(Gamma.brms.mod, type='trace')
  })
  print({
    pp = brms::pp_check(Gamma.brms.mod, type = 'ecdf_overlay', nsamples=min(iter-warmup, 200)) 
    pp + theme_bw()+ xlim(c(0,500)) + xlab('Biomass (g-m2)') + ylab("Cumulative probability")
  })
  print({
    brms::pp_check(Gamma.brms.mod, type = 'loo_pit_qq', nsamples=min(iter-warmup, 200)) 
  })
  print({
    brms::pp_check(Gamma.brms.mod, type = 'error_scatter_avg_vs_x', nsamples=min(iter-warmup, 200), x='DEPTH_c') 
  })
  dev.off()
  
}

parallel::mclapply(c(1,4,5,2,3),function(m) {
  run_depth_model(input_frame[m,], 
                  iter = 2500,
                  warmup=500)
  }, mc.cores=3, mc.preschedule = F)

save(input_frame, file = "../intermed_data/model_options.RData")

if(as.logical(Sys.getenv('CLOUD_RUN', F))){
  
  sapply(input_frame$SHORT, function(d) dir.create(file.path('/output', d)))
  
  # copy model files and diags to cloud bucket
  mod_files <- dir('../model_output/', recursive = T, full.names = T)
  sapply(mod_files, function(f) file.copy(f, file.path('/output', strsplit(f,'//')[[1]][2])))
  
  int_files <- dir('../intermed_data/', recursive = T, full.names = T)
  sapply(int_files, function(f) file.copy(f, file.path('/output', strsplit(f,'//')[[1]][2])))
  
}