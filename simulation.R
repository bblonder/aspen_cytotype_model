library(ggplot2)
library(reshape2)
library(progress)
library(dplyr)
library(data.table)
library(stringr)

do_simulation <- function(dir_out = sprintf('outputs_%s',Sys.Date()),
                          id=paste(letters[runif(20)*26],collapse=""),
                          n_indivs = 1000,
                          n_iterations = 100,
    
                          prob_survival_offspring_haploid=0,
                          prob_survival_offspring_diploid=1,
                          prob_survival_offspring_triploid=0, 
                          prob_survival_offspring_tetraploid=0, 

                          prob_gametes_haploid_0=0.0,
                          prob_gametes_haploid_1=0.0,
                          
                          prob_gametes_diploid_0=0.0,
                          prob_gametes_diploid_1=0.97,
                          prob_gametes_diploid_2=0.03,
                          
                          prob_gametes_triploid_0=0.0,
                          prob_gametes_triploid_1=0.0,
                          prob_gametes_triploid_2=0.0,
                          prob_gametes_triploid_3=0.0,
                          
                          prob_gametes_tetraploid_0=0.0,
                          prob_gametes_tetraploid_1=0.03,
                          prob_gametes_tetraploid_2=0.97,
                          prob_gametes_tetraploid_3=0.00,
                          prob_gametes_tetraploid_4=0.00,
                          
                          prob_apomixis_unreduced_gametes_haploid=0.0,
                          prob_apomixis_unreduced_gametes_diploid=0.0,
                          prob_apomixis_unreduced_gametes_triploid=0.0,
                          prob_apomixis_unreduced_gametes_tetraploid=0.0,
                          
                          n_points_to_average=10,
                          
                          num_gametes_per_parent=10,
                          
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
    prob_gametes_haploid_0, prob_gametes_haploid_1, 0, 0, 0, # n=1
    prob_gametes_diploid_0, prob_gametes_diploid_1, prob_gametes_diploid_2, 0, 0, # n=2
    prob_gametes_triploid_0, prob_gametes_triploid_1, prob_gametes_triploid_2, prob_gametes_triploid_3, 0, # n=3
    prob_gametes_tetraploid_0, prob_gametes_tetraploid_1, prob_gametes_tetraploid_2, prob_gametes_tetraploid_3, prob_gametes_tetraploid_4  # n=4
    ),
    byrow = TRUE,
    dimnames = list(1:4,0:4))
  
  prob_gametes_summed = rowSums(param_gamete_transitions,  na.rm=T)
  if(!all(prob_gametes_summed <= 1)) # make sure we don't have values that are too large
  {
    warning('not all gamete row sums <=1')
  }
  #print(param_gamete_transitions)
  
  # trim small values
  # param_gamete_transitions[abs(param_gamete_transitions) < 1e-6] = 0

  param_survival_offspring = c(prob_survival_offspring_haploid, 
                               prob_survival_offspring_diploid, 
                               prob_survival_offspring_triploid, 
                               prob_survival_offspring_tetraploid) # n=1,2,3,4
  
  param_apomixis = c(prob_apomixis_unreduced_gametes_haploid, 
                     prob_apomixis_unreduced_gametes_diploid,
                     prob_apomixis_unreduced_gametes_triploid,
                     prob_apomixis_unreduced_gametes_tetraploid)
  
  stopifnot(all(param_gamete_transitions>=0)) 
  stopifnot(all(param_gamete_transitions<=1)) 
  stopifnot(all(param_survival_offspring>=0)) 
  stopifnot(all(param_apomixis>=0))
  stopifnot(n_iterations >= n_points_to_average)
  
  # start with an all diploid population
  population = rep(2, n_indivs)

  # initialize the counts table
  output_counts_all = matrix(nrow=n_iterations,ncol=4,data=NA)
  
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
      offspring = sample(gametes,size=length(gametes), replace=FALSE) + sample(gametes,size=length(gametes), replace=FALSE)
      # cut off ploidy level at 1-4 (assume mortality of anything higher or lower)
      offspring[(offspring > 4) | (offspring==0)] = NA
      offspring = na.omit(offspring)
      # in case there were no offspring produced (e.g. no non-NA gametes)
      if(length(offspring)==0)
      {
        offspring = rep(NA, times=n_indivs)
      }
      else
      {
        # resample values within existing offspring to bring population to correct size
        offspring = sample(x=offspring, size=n_indivs, replace=TRUE)
      }
      # do apomixis/parent survival
      # i.e. overwrite offspring with parents' ploidy level
      apomixis_probs = param_apomixis[population]
      apomixis_occurred = (runif(n=length(population),min=0,max=1) < apomixis_probs)
      offspring[which(apomixis_occurred==TRUE)] = population[which(apomixis_occurred==TRUE)]
      
      # remove NA offspring
      offspring = na.omit(offspring)
      if(length(offspring)==0)
      {
        offspring = rep(NA, size=n_indivs)
      }
      else
      {
        # resample values within existing offspring to bring population to correct size
        offspring = sample(x=offspring, size=n_indivs, replace=TRUE)
      }
        
      if (all(is.na(offspring)) | all(param_survival_offspring[offspring]==0))
      {
        population_new = rep(NA, size=n_indivs)
      }
      else
      {
        # make a new adult population sampling based on survival of offspring
        population_new = sample(offspring, size=n_indivs, replace=TRUE, 
                                prob=param_survival_offspring[offspring])
      }
    }
    else
    {
      population_new = rep(NA, n_indivs)
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
  }
  
  # check output dir
  
  if (!file.exists(dir_out))
  {
    dir.create(dir_out)
  }
  
  
  # make a plot
  if (save_plot==TRUE)
  {
    df_for_plotting = reshape2::melt(output_counts_all,id.vars=1:4)
    names(df_for_plotting) = c('time','ploidy_level','count')
    g = ggplot(df_for_plotting, 
           aes(x=time,y=count/n_indivs,color=factor(ploidy_level),group=ploidy_level)) +
      geom_line() +
      theme_bw() +
      scale_color_manual(values=c('gray','blue','red','green'),name='Ploidy level') +
      xlab('Time') +
      ylab('Fraction of individuals') +
      ggtitle(str_wrap(argg_string,60)) + 
      theme(plot.title = element_text(size = 6))

    ggsave(g,file=sprintf("%s/ts_%s.pdf",dir_out,id),width=10,height=10)
  }
  
  if (save_time_series==TRUE)
  {
    output_counts_all = data.frame(output_counts_all)
    names(output_counts_all) = paste("ploidy",c(1:4),sep="_")
    write.csv(output_counts_all, file=sprintf("%s/ts_%s.csv",dir_out,id), row.names=FALSE)
  }
  
  fracs_final_mean = apply(tail(output_counts_all / n_indivs, n_points_to_average), 2, mean)
  names(fracs_final_mean) = paste("mean_n_",1:4,sep="")
 
  fracs_final_sd = apply(tail(output_counts_all / n_indivs, n_points_to_average), 2, sd)
  names(fracs_final_sd) = paste("sd_n_",1:4,sep="") 
  
  return(cbind(argg,data.frame(t(c(fracs_final_mean, fracs_final_sd)))))
}

# example simulation
# do_simulation(
#               prob_survival_offspring_haploid = 0.1,
#               prob_survival_offspring_triploid = 1,
#               prob_survival_offspring_tetraploid = 0.0,
# 
#               prob_gametes_diploid_0 = 0.1,
#               prob_gametes_diploid_1 = 0.75,
#               prob_gametes_diploid_2 = 0.05,
#               prob_gametes_triploid_3 = 0.0,
# 
#               prob_apomixis_unreduced_gametes_triploid = 0.01,
#   
#               progress = TRUE,
#               save_plot = TRUE,
#               save_time_series = TRUE,
# 
#               n_indivs = 1000,
#               n_iterations = 50,
#               n_points_to_average = 5,
# )


