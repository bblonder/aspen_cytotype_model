source('simulation.R')
library(dplyr)

simulation_to_optimize <- function(params, ...)
{
	result = do_simulation(prob_survival_offspring_triploid=params[1], 
	                       prob_gametes_diploid_2=params[2],
	                       prob_apomixis_unreduced_gametes_triploid=params[3],
	                       ...)
	result_mean = result %>% select(contains("mean_n")) %>% as.numeric
	print(result_mean)
	target = c(0,0.7,0.3,0)
	
	distance = sum((result_mean - target)^2)
	
	print(params)
	print(distance)
	return(distance)
}


w = optim(fn= simulation_to_optimize, par=c(0.1,0.1,0.1),lower=c(0,0,0),upper=c(1,1,1),method="L-BFGS-B")
