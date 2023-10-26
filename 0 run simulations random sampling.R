source('simulation.R')
library(parallel)

# set up simulations
DEBUG = FALSE

if (DEBUG==TRUE)
{
  NUM_ITERATIONS = 20
  NUM_SAMPLES = 1
  NUM_CORES = 1  
} else
{
  NUM_ITERATIONS = 100
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
                          include_offspring_survival_variation=TRUE,
                          include_haploid_tetraploid=TRUE,
                          include_triploid_fertility=TRUE,
                          include_parent_survival=TRUE)
{
  prob_gametes_1 = replicate(n=num_samples,pick_n_less_than_sum_one(n=2),simplify='matrix')
  prob_gametes_2 = replicate(n=num_samples,pick_n_less_than_sum_one(n=3),simplify='matrix')
  prob_gametes_3 = replicate(n=num_samples,pick_n_less_than_sum_one(n=4),simplify='matrix')
  prob_gametes_4 = replicate(n=num_samples,pick_n_less_than_sum_one(n=5),simplify='matrix')
  
  df_out = data.frame(
    dir_out=dir_out,
    
    prob_survival_offspring_1 = runif(n=num_samples, min=0, max=1) * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    prob_survival_offspring_2 = runif(n=num_samples, min=0, max=1),
    prob_survival_offspring_3 = runif(n=num_samples, min=0, max=1),
    prob_survival_offspring_4 = runif(n=num_samples, min=0, max=1) * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    
    prob_survival_parent_1 = runif(n=num_samples, min=0, max=1) * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_parent_survival==TRUE, 1, 0),
    prob_survival_parent_2 = runif(n=num_samples, min=0, max=1) * ifelse(include_parent_survival==TRUE, 1, 0),
    prob_survival_parent_3 = runif(n=num_samples, min=0, max=1) * ifelse(include_parent_survival==TRUE, 1, 0),
    prob_survival_parent_4 = runif(n=num_samples, min=0, max=1) * ifelse(include_haploid_tetraploid==TRUE, 1, 0) * ifelse(include_parent_survival==TRUE, 1, 0),
    
    prob_gametes_1_0 = prob_gametes_1[1,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    prob_gametes_1_1 = prob_gametes_1[2,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    
    prob_gametes_2_0 = prob_gametes_2[1,],
    prob_gametes_2_1 = prob_gametes_2[2,],
    prob_gametes_2_2 = prob_gametes_2[3,],
    
    prob_gametes_3_0 = prob_gametes_3[1,] * ifelse(include_triploid_fertility==TRUE, 1, 0),
    prob_gametes_3_1 = prob_gametes_3[2,] * ifelse(include_triploid_fertility==TRUE, 1, 0),
    prob_gametes_3_2 = prob_gametes_3[3,] * ifelse(include_triploid_fertility==TRUE, 1, 0),
    prob_gametes_3_3 = prob_gametes_3[4,] * ifelse(include_triploid_fertility==TRUE, 1, 0),
    
    prob_gametes_4_0 = prob_gametes_4[1,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    prob_gametes_4_1 = prob_gametes_4[2,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    prob_gametes_4_2 = prob_gametes_4[3,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    prob_gametes_4_3 = prob_gametes_4[4,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    prob_gametes_4_4 = prob_gametes_4[5,] * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    
    prob_apomixis_1 = runif(n=num_samples, min=0, max=1) * ifelse(include_apomixis==TRUE,1,0) * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    prob_apomixis_2 = runif(n=num_samples, min=0, max=1) * ifelse(include_apomixis==TRUE,1,0),
    prob_apomixis_3 = runif(n=num_samples, min=0, max=1) * ifelse(include_apomixis==TRUE,1,0),
    prob_apomixis_4 = runif(n=num_samples, min=0, max=1) * ifelse(include_apomixis==TRUE,1,0) * ifelse(include_haploid_tetraploid==TRUE, 1, 0),
    
    n_indivs = round(runif(n=num_samples, min=10, max=1000),digits = 0),
    
    num_gametes_per_parent=10,
    num_offspring_per_apomixis=10,
    
    n_iterations = NUM_ITERATIONS,
    n_points_to_average = 20,
    save_plot=TRUE,
    save_time_series=TRUE,
    progress=FALSE
    
    )
  
  if(include_offspring_survival_variation==FALSE)
  {
    df_out$prob_survival_offspring_1=1
    df_out$prob_survival_offspring_2=1
    df_out$prob_survival_offspring_3=1
    df_out$prob_survival_offspring_4=1
  }
  
  return(df_out)
}







run_simulations <- function(params_simulations)
{
  # run all combinations
  outputs_all = mclapply(1:nrow(params_simulations), function(i)
  {
    #system(paste("echo 'now processing:",i,"'"))
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

  if (!file.exists(params_simulations$dir_out[1]))
  {
    dir.create(params_simulations$dir_out[1])
  }
  write.csv(outputs_all_df, row.names = FALSE, file=sprintf('%s/outputs.csv',params_simulations$dir_out[1]))
}

# generate parameter combinations
df_parameters_raw = expand.grid(num_samples = NUM_SAMPLES, 
            include_triploid_fertility=c(0,1),
            include_apomixis=c(0,1),
            include_haploid_tetraploid=c(0,1),
            include_offspring_survival_variation=c(0,1),
            include_parent_survival=c(0,1))

case_name = df_parameters_raw %>%
  select(-num_samples) %>%
  mutate(dir_out = gsub("include_","",apply(., 1, function(x) { paste(names(.), x, sep="=",collapse=",") }))) %>%
  mutate(dir_out = paste("output", dir_out, sep=",")) %>%
  pull(dir_out)
         
df_parameters = df_parameters_raw %>% mutate(dir_out = case_name)   

# run each parameter combination
lapply(1:nrow(df_parameters), function(i) { 
  print(i)
  params_this = do.call("sample_params", df_parameters[i,])
  
  run_simulations(params_this)
})
