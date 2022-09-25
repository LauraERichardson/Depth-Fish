require(tidyverse)

if(as.logical(Sys.getenv('CLOUD_RUN', F))){
  load('/input/PACIFIC-FISH-DEPTH-MODELS/Depth_study_fish_data.RData', v=T)
  load('/input/PACIFIC-FISH-DEPTH-MODELS/model_options.RData') 
  mod_dir <- '/input/PACIFIC-FISH-DEPTH-MODELS/'
} else {
  load('../intermed_data/Depth_study_fish_data.RData', v=T)
  load('../intermed_data/model_options.RData')
  mod_dir <- "../model_output/"
}

####---- load models ---####

# R names of models
models <- c("TotFish.Gamma.brms.full", 
            "PRIM.Gamma.brms.full",
            "PLANK.Gamma.brms.HU.full",
            "SEC.Gamma.brms.full",
            "PISC.Gamma.brms.HU.full")

dats <- c("fish","PRIM","PLANK","SEC","PISC")

ord <- sapply(input_frame$SHORT, grep, models, ignore.case = T)

lapply(1:length(models), function(idx){
  assign(models[ord[idx]], read_rds(paste0(mod_dir,
                                           input_frame$SHORT[idx],"/",
                                           input_frame$SHORT[idx],'_mod.rds')),
         envir = globalenv()
  )
  return()
})

with(input_frame,{
  sapply(1:length(HURDLE), function(h){
    if(HURDLE[h]) {
      assign(SHORT[h], fish, envir = globalenv())
    } else {
      assign(SHORT[h], filter(fish,!!sym(RESP[h])>0) %>% droplevels, envir = globalenv())
    }
  })
})
