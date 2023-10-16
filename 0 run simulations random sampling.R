source('simulation.R')
library(parallel)

# set up simulations
DEBUG = FALSE

if (DEBUG==TRUE)
{
  NUM_SAMPLES = 10
  NUM_CORES = 1  
} else
{
  NUM_SAMPLES = 20000
  NUM_CORES = 24
}


# set up functions

pick_n_less_than_sum_one <- function(n, max_scale_factor=1)
{
  points = runif(n=n, min=0, max=1)
  points_normalized = points / sum(points)

  scale_factor = runif(n=1, min=0, max=1)
  points_normalized = points_normalized * scale_factor * max_scale_factor
  
  return(points_normalized)
}




# pick parameters

sample_params <- function(num_samples,
                          dir_out,
                          include_apomixis=TRUE,
                          include_haploid_tetraploid=TRUE,
                          include_full_range=TRUE,
                          include_survival_variation=TRUE,
                          include_triploid_fertility=TRUE)
{
  prob_gametes_haploid = replicate(n=num_samples,pick_n_less_than_sum_one(n=2),simplify='matrix')
  prob_gametes_diploid = replicate(n=num_samples,pick_n_less_than_sum_one(n=3),simplify='matrix')
  prob_gametes_triploid = replicate(n=num_samples,pick_n_less_than_sum_one(n=4),simplify='matrix')
  prob_gametes_tetraploid = replicate(n=num_samples,pick_n_less_than_sum_one(n=5),simplify='matrix')
  
  df_out = data.frame(
    dir_out=dir_out,
    prob_survival_offspring_haploid = runif(n=num_samples, min=0, max=1) * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    prob_survival_offspring_diploid = runif(n=num_samples, min=0, max=1),
    prob_survival_offspring_triploid = runif(n=num_samples, min=0, max=1),
    prob_survival_offspring_tetraploid = runif(n=num_samples, min=0, max=1) * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    
    prob_gametes_haploid_0 = prob_gametes_haploid[1,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    prob_gametes_haploid_1 = prob_gametes_haploid[2,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    
    prob_gametes_diploid_0 = prob_gametes_diploid[1,] * ifelse(include_full_range==FALSE,0.1,1),
    prob_gametes_diploid_1 = prob_gametes_diploid[2,],
    prob_gametes_diploid_2 = prob_gametes_diploid[3,] * ifelse(include_full_range==FALSE,0.1,1),
    
    prob_gametes_triploid_0 = prob_gametes_triploid[1,] * ifelse(include_triploid_fertility==TRUE, 1, 0),
    prob_gametes_triploid_1 = prob_gametes_triploid[2,] * ifelse(include_full_range==FALSE,0.1,1) * ifelse(include_triploid_fertility==TRUE, 1, 0),
    prob_gametes_triploid_2 = prob_gametes_triploid[3,] * ifelse(include_full_range==FALSE,0.1,1) * ifelse(include_triploid_fertility==TRUE, 1, 0),
    prob_gametes_triploid_3 = prob_gametes_triploid[4,] * ifelse(include_full_range==FALSE,0.1,1) * ifelse(include_triploid_fertility==TRUE, 1, 0),
    
    prob_gametes_tetraploid_0 = prob_gametes_tetraploid[1,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    prob_gametes_tetraploid_1 = prob_gametes_tetraploid[2,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    prob_gametes_tetraploid_2 = prob_gametes_tetraploid[3,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    prob_gametes_tetraploid_3 = prob_gametes_tetraploid[4,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    prob_gametes_tetraploid_4 = prob_gametes_tetraploid[5,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    
    prob_apomixis_unreduced_gametes_haploid = runif(n=num_samples, min=0, max=1) * ifelse(include_apomixis==TRUE,1,0) * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    prob_apomixis_unreduced_gametes_diploid = runif(n=num_samples, min=0, max=1) * ifelse(include_apomixis==TRUE,1,0),
    prob_apomixis_unreduced_gametes_triploid = runif(n=num_samples, min=0, max=1) * ifelse(include_apomixis==TRUE,1,0),
    prob_apomixis_unreduced_gametes_tetraploid = runif(n=num_samples, min=0, max=1) * ifelse(include_apomixis==TRUE,1,0) * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_full_range==FALSE,0.1,1),
    
    n_indivs = round(runif(n=num_samples, min=10, max=1000),digits = 0),#sample(c(10,1000),size=num_samples,replace=TRUE),
    
    num_gametes_per_parent=10,
    
    n_iterations = 100,
    n_points_to_average = 20,
    save_plot=TRUE,
    save_time_series=TRUE,
    progress=FALSE
    
    )
  
  if(include_survival_variation==FALSE)
  {
    df_out$prob_survival_offspring_haploid=1
    df_out$prob_survival_offspring_diploid=1
    df_out$prob_survival_offspring_triploid=1
    df_out$prob_survival_offspring_tetraploid=1
  }
  
  return(df_out)
}







run_simulations <- function(params_simulations)
{
  # run all combinations
  outputs_all = mclapply(1:nrow(params_simulations), function(i)
  {
    system(paste("echo 'now processing:",i,"'"))
    #print(sprintf('%.3f',i/nrow(params_simulations)))
    
    args_this = params_simulations %>% 
      slice(i) %>%
      mutate(id=i)
    
    set.seed(1)
    
    result_this = do.call("do_simulation",
                          args=args_this)
    
  
    return(result_this)
    
  }, mc.cores=NUM_CORES)
  outputs_all_df = rbindlist(outputs_all)
  outputs_all_df_combined = cbind(outputs_all_df %>% 
                                    select("id", starts_with("mean_n"),starts_with("sd_n")), 
                                  params_simulations)
  
  if (!file.exists(params_simulations$dir_out[1]))
  {
    dir.create(params_simulations$dir_out[1])
  }
  write.csv(outputs_all_df_combined, row.names = FALSE, file=sprintf('%s/outputs.csv',params_simulations$dir_out[1]))
}

# generate parameters and run
set.seed(1)
params_simulations_all = sample_params(num_samples = NUM_SAMPLES,
                                       dir_out='outputs_all')
run_simulations(params_simulations_all)

set.seed(2)
params_simulations_no_apomixis = sample_params(num_samples = NUM_SAMPLES, include_apomixis = FALSE,
                                       dir_out='outputs_no_apomixis')
run_simulations(params_simulations_no_apomixis)

set.seed(3)
params_simulations_no_haploid_tetraploid = sample_params(num_samples = NUM_SAMPLES, include_haploid_tetraploid = FALSE,
                                               dir_out='outputs_no_haploid_tetraploid')
run_simulations(params_simulations_no_haploid_tetraploid)

set.seed(4)
params_simulations_no_offspring_survival_variation = sample_params(num_samples = NUM_SAMPLES, include_survival_variation = FALSE,
                                                    dir_out='outputs_no_offspring_survival_variation')
run_simulations(params_simulations_no_offspring_survival_variation)

set.seed(5)
params_simulations_no_triploid_fertility = sample_params(num_samples = NUM_SAMPLES, include_triploid_fertility = FALSE,
                                                         dir_out='outputs_no_triploid_fertility')
run_simulations(params_simulations_no_triploid_fertility)

set.seed(6)
params_simulations_no_survival_variation_no_apomixis = sample_params(num_samples = NUM_SAMPLES, include_survival_variation = FALSE, include_apomixis = FALSE, 
                                                         dir_out='outputs_no_survival_variation_no_apomixis')
run_simulations(params_simulations_no_survival_variation_no_apomixis)

# set.seed(7)
# params_simulations_no_full_range = sample_params(num_samples = NUM_SAMPLES, include_full_range = FALSE,
#                                                          dir_out='outputs_no_full_range')
# run_simulations(params_simulations_no_full_range)

