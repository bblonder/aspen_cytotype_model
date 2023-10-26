library(ggplot2)
library(progress)
library(dplyr)
library(data.table)
library(stringr)
library(tidyr)

do_simulation <- function(dir_out = sprintf('outputs_%s',Sys.Date()),
                          id=paste(letters[runif(20)*26],collapse=""),
                          n_indivs = 1000,
                          n_iterations = 100,
                          
                          prob_survival_parent_1=0,
                          prob_survival_parent_2=0,
                          prob_survival_parent_3=0, 
                          prob_survival_parent_4=0, 
    
                          prob_survival_offspring_1=0,
                          prob_survival_offspring_2=1,
                          prob_survival_offspring_3=0, 
                          prob_survival_offspring_4=0, 

                          prob_gametes_1_0=0.0,
                          prob_gametes_1_1=0.0,
                          
                          prob_gametes_2_0=0.0,
                          prob_gametes_2_1=0.97,
                          prob_gametes_2_2=0.03,
                          
                          prob_gametes_3_0=0.0,
                          prob_gametes_3_1=0.0,
                          prob_gametes_3_2=0.0,
                          prob_gametes_3_3=0.0,
                          
                          prob_gametes_4_0=0.0,
                          prob_gametes_4_1=0.03,
                          prob_gametes_4_2=0.97,
                          prob_gametes_4_3=0.00,
                          prob_gametes_4_4=0.00,
                          
                          prob_apomixis_1=0.0,
                          prob_apomixis_2=0.0,
                          prob_apomixis_3=0.0,
                          prob_apomixis_4=0.0,
                          
                          n_points_to_average=10,
                          
                          num_gametes_per_parent=10,
                          num_offspring_per_apomixis=10,
                          
                          progress=TRUE,
                          
                          set_seed=TRUE,
                          
                          save_plot=FALSE,
                          save_time_series=FALSE)
{
  if (set_seed==TRUE)
  {
    set.seed(1)
  }
  # get arguments
  argg <- as.list(environment())
  argg_string = paste(paste(names(argg),argg,sep="="),collapse=" ")
  print(argg_string)
  
  # make a transition matrix for each state, assuming no higher cytotypes created as gametes
  param_gamete_transitions = matrix(nrow=4,ncol=5,data=c(
    prob_gametes_1_0, prob_gametes_1_1, 0, 0, 0, # n=1
    prob_gametes_2_0, prob_gametes_2_1, prob_gametes_2_2, 0, 0, # n=2
    prob_gametes_3_0, prob_gametes_3_1, prob_gametes_3_2, prob_gametes_3_3, 0, # n=3
    prob_gametes_4_0, prob_gametes_4_1, prob_gametes_4_2, prob_gametes_4_3, prob_gametes_4_4  # n=4
    ),
    byrow = TRUE,
    dimnames = list(1:4,0:4))
  
  prob_gametes_summed = rowSums(param_gamete_transitions,  na.rm=T)
  if(!all(prob_gametes_summed <= 1)) # make sure we don't have values that are too large
  {
    warning('not all gamete row sums <=1')
  }
  
  # trim small values
  # param_gamete_transitions[abs(param_gamete_transitions) < 1e-6] = 0

  param_survival_offspring = c(prob_survival_offspring_1, 
                               prob_survival_offspring_2, 
                               prob_survival_offspring_3, 
                               prob_survival_offspring_4) # n=1,2,3,4
  
  param_survival_parent = c(prob_survival_parent_1, 
                               prob_survival_parent_2, 
                               prob_survival_parent_3, 
                               prob_survival_parent_4) # n=1,2,3,4
  
  param_apomixis = c(prob_apomixis_1, 
                     prob_apomixis_2,
                     prob_apomixis_3,
                     prob_apomixis_4)
  
  stopifnot(all(param_gamete_transitions>=0)) 
  stopifnot(all(param_gamete_transitions<=1)) 
  stopifnot(all(param_survival_offspring>=0)) 
  stopifnot(all(param_survival_parent>=0)) 
  stopifnot(all(param_apomixis>=0))
  stopifnot(n_iterations >= n_points_to_average)
  
  # start with an all diploid population
  population = rep(2, n_indivs)

  # initialize the counts table
  output_counts_all = matrix(nrow=n_iterations,ncol=4,data=NA)
  output_counts_gametes = matrix(nrow=n_iterations,ncol=5,data=NA)
  output_counts_offspring_sexual = matrix(nrow=n_iterations,ncol=4,data=NA)
  output_counts_offspring_asexual = matrix(nrow=n_iterations,ncol=4,data=NA)
  output_counts_offspring_survived = matrix(nrow=n_iterations,ncol=4,data=NA)
  output_counts_parent_survived = matrix(nrow=n_iterations,ncol=4,data=NA)
  
  if (progress==TRUE)
  {
    pb = progress_bar$new(total=n_iterations)
  }
  for (i in 1:n_iterations)
  {
    if (progress==TRUE)
    {
      pb$tick()
    }
    # randomly assign gametes based on transition matrix
    if (!all(is.na(population)))
    {
      gametes = unlist(apply(param_gamete_transitions[population,], 1, function(x) { 
        if (all(x==0 | is.na(x)))
        {
          return(NA)
        }
        else
        {
          # pick a ploidy level based on the probabilities
          return(sample(0:4, size=num_gametes_per_parent, prob=x, replace=TRUE))
        }
      }))
    }
    else
    {
      gametes = rep(NA, n_indivs)
    }

    # remove non-gametes from the pool
    gametes = na.omit(gametes)
    
    if (length(gametes>0))
    {
      # randomly mate and calculate ploidy of offspring
      offspring_sexual = sample(gametes,size=length(gametes), replace=FALSE) + sample(gametes,size=length(gametes), replace=FALSE)
      # cut off ploidy level at 1-4 (assume mortality of anything higher or lower)
      offspring_sexual[(offspring_sexual > 4) | (offspring_sexual==0)] = NA
    }
    else
    {
      offspring_sexual = NULL
    }

    # do apomixis
    offspring_asexual = rep(population, num_offspring_per_apomixis)
    offspring_asexual_keep = (runif(n=length(offspring_asexual),min=0,max=1) < rep(param_apomixis[population], num_offspring_per_apomixis))
    offspring_asexual = offspring_asexual[offspring_asexual_keep==TRUE]

    # remove NA offspring
    offspring_all = na.omit(c(offspring_sexual, offspring_asexual))
    if(length(offspring_all)==0)
    {
      offspring_all = rep(NA, size=n_indivs)
    }
    else
    {
      # resample values within existing offspring to bring population to correct size
      offspring_all = sample(x=offspring_all, size=n_indivs, replace=TRUE)
    }
      
    if (all(is.na(offspring_all)) | all(param_survival_offspring[offspring_all]==0))
    {
      #population_new = rep(NA, size=n_indivs)
      offspring_survived = NULL
    }
    else
    {
      # make a new adult population sampling based on survival of offspring
      offspring_survived = sample(offspring_all, size=n_indivs, replace=TRUE, 
                              prob=param_survival_offspring[offspring_all])
      
      #population_new = offspring_survived
    }
    
    # do parent survival
    parent_survived_keep = (runif(n=length(population),min=0,max=1) < param_survival_parent[population])
    parent_survived = na.omit(population[parent_survived_keep==TRUE])

    # make the final new population
    num_offspring_to_keep = n_indivs - length(parent_survived)
    # if there are offspring
    if (num_offspring_to_keep > 0 & length(offspring_survived > 0))
    {
      offspring_survived_keep = sample(offspring_survived, size=num_offspring_to_keep, replace=TRUE)
    }
    else
    {
      offspring_survived_keep = NULL
    }
    # join the surviving parents and surviving offspring
    population_new = c(parent_survived, offspring_survived_keep)
    
    # get population to right size and check for empty
    if (length(population_new)==0 | all(is.na(population_new)))
    {
      population_new = rep(NA, n_indivs)
    }
    else
    {
      population_new = sample(population_new, size=n_indivs, replace=TRUE)
    }

    # create the next time step
    population = population_new
    
    # summarize outputs
    if (!all(is.na(population)))
    {
      output_counts_all[i,] = tabulate(population,nbins=4)
    }
    else
    {
      output_counts_all[i,] = rep(0,ncol(output_counts_all))
    }
    
    fraction_offspring = num_offspring_to_keep / n_indivs
    
    output_counts_gametes[i,] = as.vector(table(factor(gametes,levels=c(0:4))))
    output_counts_offspring_sexual[i,] = as.vector(table(factor(offspring_sexual,levels=c(1:4))))
    output_counts_offspring_asexual[i,] = as.vector(table(factor(offspring_asexual,levels=c(1:4)))) 
    output_counts_offspring_survived[i,] = as.vector(table(factor(offspring_survived,levels=c(1:4)))) * fraction_offspring # correction factor due to parent survival
    output_counts_parent_survived[i,] = as.vector(table(factor(parent_survived,levels=c(1:4))))
  }
  
  # check output dir
  if (!file.exists(dir_out))
  {
    dir.create(dir_out)
  }

  
  # make summary data frame
  output_counts_save = data.frame(output_counts_all, output_counts_gametes, output_counts_offspring_sexual, output_counts_offspring_asexual, output_counts_offspring_survived, output_counts_parent_survived)
  names(output_counts_save) = c(paste("counts_all",c(1:4),sep="_"), paste("counts_gametes",c(0:4),sep="_"), paste("counts_offspring_sexual",c(1:4),sep="_"), paste("counts_offspring_asexual",c(1:4),sep="_"), paste("counts_offspring_survived",c(1:4),sep="_"), paste("counts_parent_survived",c(1:4),sep="_"))
  
  # make summary variables
  counts_all_mean = apply(tail(output_counts_all, n_points_to_average), 2, mean)
  names(counts_all_mean) = paste("counts_all_mean",1:4,sep="_")
  
  counts_all_sd = apply(tail(output_counts_all, n_points_to_average), 2, sd)
  names(counts_all_sd) = paste("counts_all_sd",1:4,sep="_") 

  counts_gametes_mean = apply(tail(output_counts_gametes, n_points_to_average), 2, sd)
  names(counts_gametes_mean) = paste("counts_gametes_mean",0:4,sep="_")   
  
  counts_offspring_sexual_mean = apply(tail(output_counts_offspring_sexual, n_points_to_average), 2, sd)
  names(counts_offspring_sexual_mean) = paste("counts_offspring_sexual_mean",1:4,sep="_")   
  
  counts_offspring_asexual_mean = apply(tail(output_counts_offspring_asexual, n_points_to_average), 2, sd)
  names(counts_offspring_asexual_mean) = paste("counts_offspring_asexual_mean",1:4,sep="_")  
  
  counts_offspring_survived_mean = apply(tail(output_counts_parent_survived, n_points_to_average), 2, sd)
  names(counts_offspring_survived_mean) = paste("counts_offspring_survived_mean",1:4,sep="_")  
  
  counts_parent_survived_mean = apply(tail(output_counts_parent_survived, n_points_to_average), 2, sd)
  names(counts_parent_survived_mean) = paste("counts_parent_survived_mean",1:4,sep="_")  
  
  # make output data frame
  df_output = cbind(argg,data.frame(t(c(counts_all_mean, 
                                        counts_all_sd, 
                                        counts_gametes_mean, 
                                        counts_offspring_sexual_mean, 
                                        counts_offspring_asexual_mean,
                                        counts_parent_survived_mean))))
  
  # make a plot
  if (save_plot==TRUE)
  {
    df_plot = output_counts_save %>% 
      mutate(time=1:n()) %>%
      pivot_longer(!time) %>%
      mutate(type=gsub("_[0-4]","",gsub("counts_","",name))) %>%
      mutate(cytotype=gsub("[^0-9]*","",name))
    
    g_plot = ggplot(df_plot, 
                aes(x=time,y=value,color=cytotype,group=name)) +
      geom_line() +
      theme_bw() +
      xlab('Time') +
      ylab('Number') +
      facet_wrap(~type,scales='free_y') +
      scale_color_brewer(palette='Set2')
    
    ggsave(g_plot,file=sprintf("%s/ts_%s.pdf",dir_out,id),width=12,height=8)
  }
  
  if (save_time_series==TRUE)
  {
    write.csv(output_counts_save, file=sprintf("%s/ts_%s.csv",dir_out,id), row.names=FALSE)
  }
  
  # return final data frame summary
  return(df_output)
}

# example simulation
# do_simulation(
#               prob_survival_offspring_1 = 0.1,
#               prob_survival_offspring_3 = 1,
#               prob_survival_offspring_4 = 0.0,
# 
#               prob_gametes_2_0 = 0.1,
#               prob_gametes_2_1 = 0.75,
#               prob_gametes_2_2 = 0.05,
#               prob_gametes_3_3 = 0.0,
# 
#               prob_apomixis_3 = 0.1,
# 
#               prob_survival_parent_3 = 0.1,
# 
#               progress = TRUE,
#               save_plot = TRUE,
#               save_time_series = TRUE,
# 
#               n_indivs = 1000,
#               n_iterations = 100,
#               n_points_to_average = 5,
# )


