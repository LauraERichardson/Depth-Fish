# Depth-Fish
## Sequence of files to run:

 ### Data and models 
 - prepare_data.R: loads and grooms data, saves for models and plots.
 - run_models.R: attention - runs all models, may over-write existing files!
 
 ### Figures:
 - pop_effect_figs.R: figures for population level effects (Figure 1, Figure S1)
 - depth_dens_plot.R: Figure 2
 - proba_of_increase_plot.R: Different options for plotting probability of increase in biomass (Figure 3A-C, Figure S2A-B)
 - slope_dens_plot.R: Figure 4A 
 - proba_of_increase_slope.R: Figure 4B
 - ternery_plot.R: ternary variance partitioning (Figure 5B, Table S10)
 
 - plot_opts.R: global plotting options
 
 ### Tables
 
 - Table 1.R: Table 1 Proba of increase by trophic group, depth bin and pop status
 - Summary tables: Table S4 Model coefficient estimate summaries; Table S5 Max probability of effect estimates; Table S6 Unadjusted Bayesian R2 values; Table S11 Proba of difference in variation among scales
 - TablesS7_S8.R: Table S7 Proba of greater absolute increase; Table S8 Proba of greater proportionate increase (%)
