library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(ggpubr)
library(ranger)
library(pdp)
library(stringr)

# add distance targets
TARGET = c(0.0, 0.5, 0.5, 0)
DELTA_TARGET = 0.05



nice_names = c(
  dir_out='Scenario',
  in_target='In-target',
  distance='Target distance (d)',
  n_indivs='M',

  prob_apomixis_1='B1',
  prob_apomixis_2='B2',
  prob_apomixis_3='B3',
  prob_apomixis_4='B4',

  prob_gametes_1_0='A1,0',
  prob_gametes_1_1='A1,1',
  
  prob_gametes_2_0='A2,0',
  prob_gametes_2_1='A2,1',
  prob_gametes_2_2='A2,2',
  
  prob_gametes_3_0='A3,0',
  prob_gametes_3_1='A3,1',
  prob_gametes_3_2='A3,2',
  prob_gametes_3_3='A3,3',

  prob_gametes_4_0='A4,0',
  prob_gametes_4_1='A4,1',
  prob_gametes_4_2='A4,2',
  prob_gametes_4_3='A4,3',
  prob_gametes_4_4='A4,4',

  prob_survival_offspring_1='C1',
  prob_survival_offspring_2='C2',
  prob_survival_offspring_3='C3',
  prob_survival_offspring_4='C4',
  
  prob_survival_parent_1='D1',
  prob_survival_parent_2='D2',
  prob_survival_parent_3='D3',
  prob_survival_parent_4='D4'

  )


# make output directories
dir_results = 'stats'

if (!file.exists(dir_results))
{
  dir.create(dir_results)
}

# load in files
fn_all = file.path(dir(pattern='output'),'outputs.csv')
df_all = rbindlist(lapply(fn_all, read.csv))

names_predictors = df_all %>% 
  select(prob_survival_parent_1:prob_apomixis_4, n_indivs) %>% 
  names

# make nicer name
rewrite_names <- function(name)
{
  names_shortened = sapply(strsplit(name,',')[[1]][-1], function(x) { 
    ifelse(length(grep("=0",x)>0),NA,x)
  })
  names_final = paste(na.omit(names_shortened), collapse=", ")
  names_final = gsub("_", " ", names_final)
  names_final = gsub("=1", "", names_final)
  #names_final = gsub("^$","gametes only", names_final)
  
  return(names_final)
}

rewritten_names = data.frame(dir_out=unique(df_all$dir_out),
                             scenario=sapply(unique(df_all$dir_out), rewrite_names))
rewritten_names$nvar = sapply(rewritten_names$dir_out, str_count, "=1")
rewritten_names$scenario = factor(rewritten_names$scenario, levels=rewritten_names$scenario[order(rewritten_names$nvar)],ordered=TRUE)

df_all = df_all %>% left_join(rewritten_names,by='dir_out')



# add fractional data
data_fractions_all = df_all %>% 
  select(contains("counts_all_mean"))
data_fractions_all = apply(data_fractions_all,1, function(x) x/sum(x)) %>% 
  t %>% 
  as.data.frame
names(data_fractions_all) = paste("frac_all_mean",1:4,sep="_")
# put in 0s for NaN cases
data_fractions_all[is.na(data_fractions_all)] <- 0
df_all = cbind(df_all, data_fractions_all)


# add distance to target
data_fractions_all_target = matrix(nrow=nrow(data_fractions_all),ncol=ncol(data_fractions_all),data=TARGET,byrow=TRUE)

df_all$distance = apply(data_fractions_all - data_fractions_all_target, 1, function(x) { sum(x^2)  })
df_all$in_target = df_all$distance < DELTA_TARGET
# anything else that didn't work out doesn't count
df_all$in_target[is.na(df_all$in_target)] = 0
df_all$in_target = factor(df_all$in_target)



# make target graph
# figure out number of parameters per scenario


df_variances = df_all %>% 
  select(scenario, all_of(names_predictors)) %>% 
  group_by(scenario) %>% 
  summarize_if(is.numeric, var)

num_scenarios = df_all %>% 
  group_by(scenario) %>% 
  tally() %>%
  ungroup %>%
  pull(n) %>%
  max

num_vars = ((df_variances %>% 
  select_if(is.numeric)) > 0) %>%
  rowSums



df_counts = df_all %>% group_by(scenario, in_target, .drop=FALSE) %>%
  summarize(count=n()) %>%
  left_join(data.frame(scenario=df_variances$scenario, num_vars=num_vars),by='scenario') %>%
  filter(in_target==1) %>%
  mutate(frac=count/num_scenarios) %>%
  mutate(frac.per.var = frac^(1/num_vars))

ordering_counts = df_counts %>% 
  arrange(frac.per.var) %>% 
  pull(scenario)

df_counts = df_counts %>%
  mutate(scenario=factor(scenario,levels=ordering_counts,ordered=TRUE))


g_in_target = df_counts %>%
  select(scenario, frac, frac.per.var) %>%
  pivot_longer(!scenario) %>%
  mutate(name=factor(name,levels=c('frac','frac.per.var'),labels=c('No','Yes'))) %>%
  ggplot(aes(x=scenario,y=value,fill=name)) +
  geom_bar(stat='identity',position='dodge') +
  coord_flip() +
  scale_fill_manual(values=c('darkgoldenrod4','cornsilk3'),name='Normalized') +
  ylab('Fraction of cases') +
  xlab('Scenario') +
  theme_bw() +
  theme(legend.position='bottom')
ggsave(g_in_target, file=sprintf('%s/g_in_target.pdf',dir_results),width=10,height=6)
ggsave(g_in_target, file=sprintf('%s/g_in_target.png',dir_results),width=10,height=6)


# summarize results
df_counts %>% summary


#show mean_n
df_outcomes = df_all %>% 
  select(contains("frac_all_mean"), distance, scenario) %>% 
  group_by(scenario) %>%
  sample_n(500) %>%
  mutate(row=1:n()) %>%
  ungroup %>%
  pivot_longer(!c(distance,scenario,row)) %>%
  mutate(name=gsub("frac_all_mean_","",name))

g_mean_n = ggplot(df_outcomes %>% mutate(name=gsub("mean_n_","",name)), aes(x=name,y=value,color=distance,group=row)) +
  geom_line(alpha=0.1) + 
  geom_jitter(width=0.1) +
  geom_point(alpha=0.1,size=0.25) +
  theme_bw() +
  facet_wrap(~scenario,ncol=4,labeller = labeller(scenario = label_wrap_gen(width = 60))) +
  scale_color_viridis_c(direction=-1,name='Distance (d)') +
  xlab("Cytotype") +
  ylab("Frequency")
ggsave(g_mean_n, file=sprintf('%s/g_mean_n.pdf',dir_results),width=15,height=12)
ggsave(g_mean_n, file=sprintf('%s/g_mean_n.png',dir_results),width=15,height=12)





 # make histogram
num_cases = df_all %>% group_by(scenario) %>% tally() %>% summarize(max.n=max(n)) %>% pull(max.n)

g_hist = ggplot(df_all, aes(x=distance)) +
  facet_wrap(~scenario,scales='free_y',ncol=4,labeller = labeller(scenario = label_wrap_gen(width = 60))) +
  theme_bw() +
  geom_histogram(bins=100) +
  geom_vline(xintercept = DELTA_TARGET, color='red') +
  theme(legend.position='none') +
  xlab('Target distance (d)') + 
  ylab('Number of cases') +
  geom_text(data=df_counts %>% filter(in_target==FALSE) %>% mutate(frac=sprintf('%.2f%% in-target',100*(1-count/num_cases))),
            aes(x=1,y=200,label=frac),color='blue')
ggsave(g_hist, file=sprintf('%s/g_hist.pdf',dir_results),width=15,height=12)
ggsave(g_hist, file=sprintf('%s/g_hist.png',dir_results),width=15,height=12)






# look at 1-d histograms
# look for independent inferences on predictors
df_longer = df_all %>% 
  select(distance, in_target, scenario, all_of(names_predictors)) %>% 
  pivot_longer(!c(distance,scenario,in_target)) %>% 
  mutate(nice_name=nice_names[name]) %>% 
  mutate(in_target=factor(in_target))

g_density_new = ggplot(df_longer, aes(x=scenario, y=value,color=in_target,fill=in_target)) +
  facet_wrap(~name,scales='free',ncol=4, labeller=as_labeller(nice_names)) +
  geom_violin(draw_quantiles=c(0.5),alpha=0.5) + 
  coord_flip() +
  theme_bw() +
  xlab('Scenario') + ylab('Value') +
  scale_color_manual(values=c('gray','blue'),name='In-target') + 
  scale_fill_manual(values=c('gray','blue'),name='In-target')

ggsave(g_density_new, file=sprintf('%s/g_density_new.pdf',dir_results),width=30,height=40)
ggsave(g_density_new, file=sprintf('%s/g_density_new.png',dir_results),width=30,height=40)


g_density_all = ggplot(df_longer %>% 
                         filter(scenario=="triploid fertility, apomixis, haploid tetraploid, offspring survival variation, parent survival") %>%
                         mutate(nice_name=nice_names[name]), 
                       aes(x=nice_name, y=value,color=in_target,fill=in_target)) +
  geom_violin(draw_quantiles=c(0.5),alpha=0.5) + 
  facet_wrap(~nice_name,scales='free') +
  coord_flip() +
  theme_bw() +
  xlab('Scenario') + ylab('Value') +
  scale_color_manual(values=c('gray','blue'),name='In-target') + 
  scale_fill_manual(values=c('gray','blue'),name='In-target') +
  theme(axis.title.y=element_blank(),axis.text.y=element_blank()) +
  ggtitle("triploid fertility, apomixis, haploid tetraploid, offspring survival variation, parent survival")
ggsave(g_density_all, file=sprintf('%s/g_density_all.pdf',dir_results),width=9,height=8)
ggsave(g_density_all, file=sprintf('%s/g_density_all.png',dir_results),width=9,height=8)


# find cases to plot
id_example = df_all %>% 
  mutate(row_id=1:n()) %>%
  filter(dir_out=="output,triploid_fertility=1,apomixis=1,haploid_tetraploid=1,offspring_survival_variation=1,parent_survival=1") %>%
  filter(n_indivs > 800 & distance<0.05 & prob_gametes_2_1 > 0.5)
print(id_example[1,] %>% select(id, row_id))

plot_ts <- function(fn_ts)
{
  nice_names_time_series = c(all='Parents, P(t)',
                             gametes='Gametes, G(t)',
                             offspring_asexual='Offspring (asexual), Oasexual(t)',
                             offspring_sexual='Offspring (sexual), Osexual(t)',
                             offspring_survived='Offspring (survived), Os(t)',
                             parent_survived='Parents (survived), Ps(t)')
  
  df_ts = read.csv(fn_ts)
  
  df_ts = df_ts %>% 
    as.data.frame %>%
    mutate(time=1:n()) %>%
    pivot_longer(!time) %>%
    mutate(type=gsub("_[0-4]","",gsub("counts_","",name))) %>%
    mutate(cytotype=gsub("[^0-9]*","",name))

  g = ggplot(df_ts,
                  aes(x=time,y=value,color=cytotype,group=name)) +
    geom_line() +
    theme_bw() +
    xlab('Time') +
    ylab('Number') +
    facet_wrap(~type,scales='free_y',labeller=as_labeller(nice_names_time_series)) +
    scale_color_manual(values=c('black','blue','red','orange','gray'),drop=FALSE,name='Cytotype')
  
  return(g)
}
g_ts_in_target = plot_ts(sprintf('output,triploid_fertility=1,apomixis=1,haploid_tetraploid=1,offspring_survival_variation=1,parent_survival=1/ts_%d.csv',id_example$id[1]))

# plot parameters
plot_matrix <- function(m, subtract_one=FALSE)
{
  df = m %>%
    as.data.frame.table %>%
    mutate_if(is.factor, as.integer) %>%
    mutate(i=Var2,j=Var1)
  
  if (subtract_one==TRUE)
  {
    df$i = df$i - 1
  }
  
  g = ggplot(df, aes(x=i,y=j,fill=Freq,label=sprintf("%.2f",Freq))) +
    geom_tile() +
    scale_y_reverse() +
    geom_label(fill=gray(1,0.5),size=2) +
    theme_bw() +
    scale_fill_viridis_c(option='viridis',name='Value',limits=c(0,1)) +
    coord_equal() + 
    theme(axis.title.x=element_blank(),axis.title.y=element_blank()) +
    theme(legend.position='none')
  
  if (nrow(m)==1)
  {
    g = g + theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),axis.title.y=element_blank())
  }
  
  return(g)
}

plot_params <- function(row_id_this, df_this)
{
  A_this = matrix(nrow=4,ncol=5,data=0)
  A_this[1,1] = df_this$prob_gametes_1_0[row_id_this]
  A_this[1,2] = df_this$prob_gametes_1_1[row_id_this]
  A_this[2,1] = df_this$prob_gametes_2_0[row_id_this]
  A_this[2,2] = df_this$prob_gametes_2_1[row_id_this]
  A_this[2,3] = df_this$prob_gametes_2_2[row_id_this]
  A_this[3,1] = df_this$prob_gametes_3_0[row_id_this]
  A_this[3,2] = df_this$prob_gametes_3_1[row_id_this]
  A_this[3,3] = df_this$prob_gametes_3_2[row_id_this]
  A_this[3,4] = df_this$prob_gametes_3_3[row_id_this]
  A_this[4,1] = df_this$prob_gametes_4_0[row_id_this]
  A_this[4,2] = df_this$prob_gametes_4_1[row_id_this]
  A_this[4,3] = df_this$prob_gametes_4_2[row_id_this]
  A_this[4,4] = df_this$prob_gametes_4_3[row_id_this]
  A_this[4,5] = df_this$prob_gametes_4_4[row_id_this]

  B_this = matrix(nrow=1,ncol=4,data=0)
  B_this[1,1] = df_this$prob_apomixis_1[row_id_this]
  B_this[1,2] = df_this$prob_apomixis_2[row_id_this]
  B_this[1,3] = df_this$prob_apomixis_3[row_id_this]
  B_this[1,4] = df_this$prob_apomixis_4[row_id_this]
  
  C_this = matrix(nrow=1,ncol=4,data=0)
  C_this[1,1] = df_this$prob_survival_offspring_1[row_id_this]
  C_this[1,2] = df_this$prob_survival_offspring_2[row_id_this]
  C_this[1,3] = df_this$prob_survival_offspring_3[row_id_this]
  C_this[1,4] = df_this$prob_survival_offspring_4[row_id_this]
  
  D_this = matrix(nrow=1,ncol=4,data=0)
  D_this[1,1] = df_this$prob_survival_parent_1[row_id_this]
  D_this[1,2] = df_this$prob_survival_parent_2[row_id_this]
  D_this[1,3] = df_this$prob_survival_parent_3[row_id_this]
  D_this[1,4] = df_this$prob_survival_parent_4[row_id_this]
  
  g_A_this = plot_matrix(A_this,subtract_one = TRUE) + ggtitle('Gamete probabilities (A)')
  g_B_this = plot_matrix(B_this) + ggtitle('Apomixis probabilities (B)')
  g_C_this = plot_matrix(C_this) + ggtitle('Offspring survival probabilities (C)')
  g_D_this = plot_matrix(D_this) + ggtitle('Parent survival probabilities (D)')
  
  g_final = ggarrange(g_A_this, g_B_this, g_C_this, g_D_this,
            nrow=4,
            heights=c(2,1,1,1))
  
  return(g_final)
}

params_in_target = plot_params(id_example$row_id[1],df_all)

g_ts_example = ggarrange(params_in_target, g_ts_in_target,
                         labels=c("(a)","(b)"),
                         widths=c(2,6))
ggsave(g_ts_example, file=sprintf('%s/g_ts_example.pdf',dir_results),width=13,height=10)
ggsave(g_ts_example, file=sprintf('%s/g_ts_example.png',dir_results),width=13,height=10)






# rename and split data
df_all_split = df_all %>% 
  group_by(scenario) %>%
  select(scenario, in_target, distance, all_of(names_predictors))
# split into pieces
df_all_split = df_all_split %>% 
  group_split()

# fit random forests
rf_list = lapply(X=df_all_split, FUN=function(df_this) {
  cat('.')
  
  freqs = table(df_this$in_target)
  print(freqs)
  
  m_rf = ranger(formula(sprintf("distance ~ %s",paste(names_predictors,collapse=" + "))), 
                data=df_this,
                mtry=5,
                num.trees=1000,
                verbose=TRUE,
                #classification=TRUE,
                importance='impurity',
                #case.weights=1/df_this$distance
                case.weights=ifelse(df_this$in_target=="1",freqs["1"]/freqs["0"],1)
                )
  return(m_rf)
})
names(rf_list) = sapply(df_all_split, function(x) { x$scenario[1] })
#saveRDS(rf_list, file='rf_list.Rdata')

# post-hoc explanation of relationship among variables
df_importance = cbind(scenario=unique(df_all$scenario), do.call("rbind",lapply(rf_list, function(x) { data.frame(t(importance(x))) }))) %>%
  pivot_longer(!scenario) %>%
  group_by(scenario) %>%
  mutate(value_relative = value/max(value,na.rm=T)) %>%
  mutate(nice_name = nice_names[name])

g_rf_importance = ggplot(df_importance, aes(x=nice_name,y=value,fill=value_relative)) +
  geom_bar(stat='identity') +
  facet_wrap(~scenario,scales='free',ncol=4,labeller = labeller(scenario = label_wrap_gen(width = 60))) +
  theme_bw() +
  coord_flip() +
  scale_fill_viridis_c(direction=-1) + 
  scale_x_discrete(limits = rev(sort(unique(df_importance$nice_name)))) +
  ylab('Importance (Gini impurity)') +
  xlab('Parameter') +
  theme(legend.position='none')
ggsave(g_rf_importance, file=sprintf('%s/g_rf_importance.pdf',dir_results),width=15,height=25)
ggsave(g_rf_importance, file=sprintf('%s/g_rf_importance.png',dir_results),width=15,height=25)


g_rf_importance_all = ggplot(df_importance %>%
                               filter(scenario=="triploid fertility, apomixis, haploid tetraploid, offspring survival variation, parent survival"), 
                             aes(x=nice_name,y=value,fill=value_relative)) +
  geom_bar(stat='identity') +
  theme_bw() +
  coord_flip() +
  scale_fill_viridis_c(direction=-1) + 
  scale_x_discrete(limits = rev(sort(unique(df_importance$nice_name)))) +
  ylab('Importance (Gini impurity)') +
  xlab('Parameter') +
  theme(legend.position='none') + 
  ggtitle("triploid fertility, apomixis, haploid tetraploid, offspring survival variation, parent survival")
ggsave(g_rf_importance_all, file=sprintf('%s/g_rf_importance_all.pdf',dir_results),width=9,height=7)
ggsave(g_rf_importance_all, file=sprintf('%s/g_rf_importance_all.png',dir_results),width=9,height=7)



# get r2 values
rf_r2 = sapply(rf_list, function(x) {x$r.squared}) %>%
  as.data.frame
rf_r2$scenario = factor(row.names(rf_r2),levels=levels(ordering_counts),ordered=TRUE)
names(rf_r2)[1] = "r.squared"

g_rf_r2 = ggplot(rf_r2, aes(x=scenario,y=r.squared)) +
  geom_bar(stat='identity') +
  coord_flip() +
  theme_bw() +
  xlab('Scenario') +
  ylab('Random forest R2')
ggsave(g_rf_r2, file=sprintf('%s/g_rf_r2.pdf',dir_results),width=9,height=7)
ggsave(g_rf_r2, file=sprintf('%s/g_rf_r2.png',dir_results),width=9,height=7)

rf_r2$r.squared %>% summary

# plot params

make_dot_plot <- function(scenario_this="triploid fertility, apomixis, haploid tetraploid, offspring survival variation, parent survival",
                          xvar,
                          yvar,
                          include_log_limits=TRUE,
                          include_lines=TRUE,
                          include_title=TRUE,
                          xlab_this,
                          ylab_this)
{
  df_this = df_all %>%
    mutate(fecundity_4 = prob_gametes_4_1 + prob_gametes_4_2 + prob_gametes_4_3 + prob_gametes_4_4) %>%
    mutate(fecundity_3 = prob_gametes_3_1 + prob_gametes_3_2 + prob_gametes_3_3) %>%
    mutate(fecundity_2 = prob_gametes_2_1 + prob_gametes_2_2) %>%
    filter(scenario==scenario_this)
  
  g = ggplot(df_this, aes_string(x=xvar,y=yvar,color="distance")) +
    geom_point(alpha=0.5,aes(size=ifelse(distance<0.05,1.5,0.5))) +
    theme_bw() +
    scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.025,0.1,0.5,1),name='Distance (d)',limits=c(0,1.75)) +
    scale_size(name='In-target',guide='none',range=c(0.5,3)) +
    xlab(xlab_this) + 
    ylab(ylab_this)
  
  if (include_title==TRUE)
  {
    g = g + ggtitle(str_wrap(scenario_this, 70))
  }
  
  if (include_log_limits==TRUE)
  {
    g = g + 
      scale_x_log10(limits=c(1e-2,1e2)) +
      scale_y_log10(limits=c(1e-2,1e2))
  }
  if (include_lines==TRUE)
  {
    g = g + 
      geom_hline(yintercept = 1) +
      geom_vline(xintercept = 1)
  }
  
  return(g)
}



# all
g_dots_all_1 = make_dot_plot(scenario_this="triploid fertility, apomixis, haploid tetraploid, offspring survival variation, parent survival",
                                include_title = TRUE,
                                xvar="prob_survival_parent_3/prob_survival_parent_2",
                                yvar="prob_survival_offspring_3/prob_survival_offspring_2",
                                xlab_this="D3/D2",
                                ylab_this="C3/C2")

g_dots_all_2 = make_dot_plot(scenario_this="triploid fertility, apomixis, haploid tetraploid, offspring survival variation, parent survival",
                                include_title = FALSE,
                                 xvar="fecundity_3/fecundity_2",
                                 yvar="prob_survival_offspring_3/prob_survival_offspring_2",
                                 xlab_this="Relative triploid fecundity, (A31 + A32 + A33)/(A21+A22)",
                                 ylab_this="C3/C2")

g_dots_realistic_1 = make_dot_plot(scenario_this="triploid fertility, offspring survival variation, parent survival",
                                   include_title = TRUE,
                             xvar="prob_survival_parent_3/prob_survival_parent_2",
                             yvar="prob_survival_offspring_3/prob_survival_offspring_2",
                             xlab_this="D3/D2",
                             ylab_this="C3/C2")

g_dots_realistic_2 = make_dot_plot(scenario_this="triploid fertility, offspring survival variation, parent survival",
                                   include_title = FALSE,
                             xvar="fecundity_3/fecundity_2",
                             yvar="prob_survival_offspring_3/prob_survival_offspring_2",
                             xlab_this="Relative triploid fecundity, (A31 + A32 + A33)/(A21+A22)",
                             ylab_this="C3/C2")

ggsave(ggarrange(g_dots_all_1, g_dots_all_2, 
                 g_dots_realistic_1, g_dots_realistic_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels='auto',
                 legend='bottom'),
       file=sprintf('%s/g_dots_all.pdf',dir_results), 
       width=12,height=13) 
ggsave(ggarrange(g_dots_all_1, g_dots_all_2, 
                 g_dots_realistic_1, g_dots_realistic_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels='auto',
                 legend='bottom'),
       file=sprintf('%s/g_dots_all.png',dir_results), 
       width=12,height=13) 



g_dots_base_1 = make_dot_plot(scenario_this="",
                              include_title = TRUE,
                                   xvar="prob_gametes_2_2",
                                   yvar="prob_gametes_2_1",
                                   include_log_limits = FALSE,
                                   xlab_this="A22",
                                   ylab_this="A21")

g_dots_base_2 = make_dot_plot(scenario_this="triploid fertility",
                              include_title = TRUE,
                              xvar="fecundity_3",
                              yvar="fecundity_2",
                              include_log_limits = FALSE,
                              xlab_this="Triploid fecundity, (A31 + A32 + A33)",
                              ylab_this="Diploid fecundity, (A21 + A22)")

g_dots_base_3 = make_dot_plot(scenario_this="haploid tetraploid",
                              include_title = TRUE,
                              xvar="prob_gametes_4_2",
                              yvar="prob_gametes_2_2",
                              include_log_limits = FALSE,
                              xlab_this="A42",
                              ylab_this="A22")

g_dots_base_4 = make_dot_plot(scenario_this="apomixis",
                              include_title = TRUE,
                              xvar="prob_apomixis_3/prob_apomixis_2",
                              yvar="prob_gametes_2_2/prob_gametes_2_1",
                              include_log_limits = TRUE,
                              xlab_this="B3/B2",
                              ylab_this="A22/A21")

ggsave(ggarrange(g_dots_base_1, g_dots_base_2, g_dots_base_3, g_dots_base_4,
          align='hv',
          common.legend = TRUE,
          labels='auto',
          legend='bottom'),      
      file=sprintf('%s/g_dots_base.pdf',dir_results), 
      width=12,height=13) 
ggsave(ggarrange(g_dots_base_1, g_dots_base_2, g_dots_base_3, g_dots_base_4,
                 align='hv',
                 common.legend = TRUE,
                 labels='auto',
                 legend='bottom'),      
       file=sprintf('%s/g_dots_base.png',dir_results), 
       width=12,height=13) 




# do PDPs
# do_pdp_2d <- function(m_rf, pred_x, pred_y, grid_resolution=5)
# {
#   grid = expand.grid(x=seq(0,1,length.out=grid_resolution),y=seq(0,1,length.out=grid_resolution))
#   names(grid) = c(pred_x, pred_y)
#   
#   pdp_this = partial(m_rf,
#               pred.var=c(pred_x,pred_y),
#               pred.grid=grid,
#               train=df_all %>% filter(dir_out=='all') %>% sample_n(1000),
#               progress = TRUE,
#               plot = FALSE)
#   
#   g_this = ggplot(pdp_this, aes_string(x=pred_x,y=pred_y,fill="yhat")) + 
#     geom_raster() +
#     scale_fill_viridis_c(name='Distance (d)') +
#     theme_bw() +
#     xlab(nice_names[pred_x]) +
#     ylab(nice_names[pred_y])
#   
#   return(g_this)
# }

do_pdp_1d <- function(scenario_this, pred_x, grid_resolution=5)
{
  grid = expand.grid(x=seq(0,1,length.out=grid_resolution))
  names(grid) = c(pred_x)
  
  pdp_this = partial(rf_list[[scenario_this]],
                     pred.var=c(pred_x),
                     pred.grid=grid,
                     train=df_all %>% 
                       filter(scenario==scenario_this) %>%
                       sample_n(1000), # for speed
                     progress = TRUE,
                     plot = FALSE)
  
  return(pdp_this)
}

# make PDPs
make_pdps_1d_for_scenario <- function(scenario_this)
{
  pdps_1d = lapply(setdiff(names_predictors,"n_indivs"), do_pdp_1d,
         scenario_this=scenario_this, 
         grid_resolution = 20)
  
  # plot PDPs
  pdps_1d_joined = rbindlist(lapply(pdps_1d, function(x) { 
    df = data.frame(x,var=nice_names[names(x)[1]])
    names(df) = c("x","yhat","var")
    return(df) }))
  pdps_1d_joined = pdps_1d_joined %>%
    group_by(var) %>%
    mutate(is.min=factor(yhat==min(yhat))) %>%
    mutate(group=substr(var,1,1))
  
  pdps_1d_list = by(pdps_1d_joined, pdps_1d_joined$group, function(x) {
    g_pdps_1d_this = ggplot(x, aes(x=x,y=yhat,color=is.min,group=var)) +
      theme_bw() +
      geom_point() +
      geom_line(color='black') +
      facet_wrap(~var,ncol=4) +
      scale_color_manual(values=c('black','red')) +
      theme(legend.position='none') +
      ylab("Distance (d)") +
      xlab("Value")
  })
  
  ggsave(ggarrange(plotlist=pdps_1d_list,ncol=1,labels='auto',heights=c(3,1,1,1)),
         file=sprintf('%s/g_pdps_1d_%s.png',dir_results,scenario_this), 
         width=8,height=15) 
}

make_pdps_1d_for_scenario("triploid fertility, apomixis, haploid tetraploid, offspring survival variation, parent survival")

