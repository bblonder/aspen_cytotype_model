library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(rgl)
library(ggpubr)
library(ranger)
library(pdp)

# add distance targets
TARGET = c(0.0, 0.5, 0.5, 0)
DELTA_TARGET = 0.05

nice_names = c(
  dir_out='Scenario',
  in_target='In-target',
  distance='Distance',
  n_indivs='M',
  prob_apomixis_unreduced_gametes_diploid='B2',
  prob_apomixis_unreduced_gametes_haploid='B1',
  prob_apomixis_unreduced_gametes_tetraploid='B4',
  prob_apomixis_unreduced_gametes_triploid='B3',
  prob_gametes_diploid_0='A2,0',
  prob_gametes_diploid_1='A2,1',
  prob_gametes_diploid_2='A2,2',
  prob_gametes_haploid_0='A1,0',
  prob_gametes_haploid_1='A1,1',
  prob_gametes_tetraploid_0='A4,0',
  prob_gametes_tetraploid_1='A4,1',
  prob_gametes_tetraploid_2='A4,2',
  prob_gametes_tetraploid_3='A4,3',
  prob_gametes_tetraploid_4='A4,4',
  prob_gametes_triploid_0='A3,0',
  prob_gametes_triploid_1='A3,1',
  prob_gametes_triploid_2='A3,2',
  prob_gametes_triploid_3='A3,3',
  prob_survival_offspring_diploid='C2',
  prob_survival_offspring_haploid='C1',
  prob_survival_offspring_tetraploid='C4',
  prob_survival_offspring_triploid='C3'
  )


# make output directories
dir_results = 'stats'

if (!file.exists(dir_results))
{
  dir.create(dir_results)
}

# load in files
fn_all = file.path(dir(pattern='outputs*'),'outputs.csv')
df_all = rbindlist(lapply(fn_all, read.csv)) %>%
  mutate(dir_out=gsub("outputs_","",dir_out))



data_mean_n = df_all %>% 
  select(contains("mean_n"))
stopifnot(ncol(data_mean_n)==4)

data_mean_n_target = matrix(nrow=nrow(data_mean_n),ncol=ncol(data_mean_n),data=TARGET,byrow=TRUE)

df_all$distance = apply(data_mean_n - data_mean_n_target, 1, function(x) { sum(x^2)  })
df_all$in_target = df_all$distance < DELTA_TARGET

df_counts = df_all %>% group_by(dir_out, in_target, .groups='keep') %>%
  summarize(count=n())

g_in_target = ggplot(df_counts,aes(x=dir_out,y=count,fill=in_target)) +
  geom_bar(stat='identity') +
  #facet_wrap(~dir_out) + 
  theme_bw() +
  coord_flip() +
  theme(legend.position='none')
ggsave(g_in_target, file=sprintf('%s/g_in_target.pdf',dir_results),width=8,height=6)



#show mean_n
df_outcomes = df_all %>% 
  select(contains("mean_n"), distance, dir_out) %>% 
  group_by(dir_out) %>%
  sample_n(500) %>%
  mutate(row=1:n()) %>%
  ungroup %>%
  pivot_longer(!c(distance,dir_out,row))

g_mean_n = ggplot(df_outcomes %>% mutate(name=gsub("mean_n_","",name)), aes(x=name,y=value,color=distance,group=row)) +
  geom_line(alpha=0.1) + 
  geom_jitter(width=0.1) +
  geom_point(alpha=0.1,size=0.25) +
  theme_bw() +
  facet_wrap(~dir_out) +
  scale_color_viridis_c(direction=-1,name='Distance (d)') +
  xlab("Cytotype") +
  ylab("Frequency")
ggsave(g_mean_n, file=sprintf('%s/g_mean_n.pdf',dir_results),width=10,height=6)
ggsave(g_mean_n, file=sprintf('%s/g_mean_n.png',dir_results),width=10,height=6)





 # make histogram
g_hist = ggplot(df_all, aes(x=distance)) +
  facet_wrap(~dir_out,scales='free_y') +
  theme_bw() +
  geom_histogram(bins=100) +
  geom_vline(xintercept = DELTA_TARGET, color='red') +
  theme(legend.position='none') +
  xlab('Distance from target (d)') + 
  ylab('Number of cases') +
  geom_text(data=df_counts %>% filter(in_target==FALSE) %>% mutate(frac=sprintf('%.2f%% in-target',100*(1-count/20000))),
            aes(x=1,y=200,label=frac),color='blue')
ggsave(g_hist, file=sprintf('%s/g_hist.pdf',dir_results),width=10,height=5)
ggsave(g_hist, file=sprintf('%s/g_hist.png',dir_results),width=10,height=5)

# ggplot(df_all, 
#        aes(x=mean_n_2,y=mean_n_3,color=distance)) + 
#   geom_point() + 
#   facet_wrap(~dir_out) +
#   scale_color_viridis_c() +
#   theme_bw()





# look at 1-d histograms
# look for independent inferences on predictors
df_quantitative = df_all %>% 
  select(mean_n_1:mean_n_4, 
         distance,
         in_target,
         dir_out,
         prob_survival_offspring_haploid:n_indivs)

names_predictors = df_quantitative %>% 
  select(prob_survival_offspring_haploid:prob_apomixis_unreduced_gametes_tetraploid, n_indivs) %>% 
  names

df_longer = df_quantitative %>% 
  select(distance, in_target, dir_out, all_of(names_predictors)) %>% 
  pivot_longer(!c(distance,dir_out,in_target))
# 
# g_density = ggplot(df_longer, aes(x=in_target,y=value,color=in_target,fill=in_target)) +
#   facet_wrap(~dir_out+name,scales='free',nrow=length(unique(df_all$dir_out))) +
#   geom_violin(draw_quantiles = c(0.5),alpha=0.5) +
#   theme_bw() + 
#   scale_color_manual(values=c('gray','blue')) + 
#   scale_fill_manual(values=c('gray','blue')) 
#   
# ggsave(g_density, file=sprintf('%s/g_density.pdf',dir_results),width=70,height=25, limitsize=FALSE)


g_density_new = ggplot(df_longer %>% mutate(nice_name=nice_names[name]), aes(color=in_target,fill=in_target,x=dir_out,y=value)) +
  facet_wrap(~nice_name,scales='free',ncol=4) +
  geom_violin(alpha=0.5,draw_quantiles = c(0.5)) +
  coord_flip() +
  theme_bw() +
  xlab('Scenario') + ylab('Value') +
  scale_x_discrete(drop=FALSE,name='In-target') +
  scale_color_manual(values=c('gray','blue')) + 
  scale_fill_manual(values=c('gray','blue'))
ggsave(g_density_new, file=sprintf('%s/g_density_new.pdf',dir_results),width=15,height=15)
ggsave(g_density_new, file=sprintf('%s/g_density_new.png',dir_results),width=15,height=15)



# find cases to plot
id_example_23 = df_all %>% 
  filter(dir_out=='all') %>%
  filter(distance<0.01) %>%
  slice_head(n=1)
print(id_example_23)

id_example_4 = df_all %>% 
  filter(dir_out=='all') %>%
  filter(mean_n_4>0.7) %>%
  slice_head(n=1)
print(id_example_4)

plot_ts <- function(fn_ts)
{
  df_ts = read.csv(fn_ts)
  
  df_ts = apply(df_ts, 1, function(x) {x/sum(x)})
  
  df_ts = df_ts %>% 
    t %>%
    as.data.frame %>%
    mutate(time=1:n()) %>%
    pivot_longer(!time) %>%
    mutate(name = gsub("ploidy_","",name))
  
  g = ggplot(df_ts, 
             aes(x=time,y=value,color=name,group=name)) +
    geom_line(linewidth=2,alpha=0.75) +
    theme_bw() +
    scale_color_manual(values=c('black','blue','red','orange'),name='Cytotype') +
    xlab('Time') +
    ylab('Frequency') +
    ylim(0,1)
  return(g)
}
g_ts_in_target_0 = plot_ts('outputs_all/ts_36.csv')
g_ts_in_target_1 = plot_ts('outputs_all/ts_623.csv')

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
  A_this[1,1] = df_this$prob_gametes_haploid_0[row_id_this]
  A_this[1,2] = df_this$prob_gametes_haploid_1[row_id_this]
  A_this[2,1] = df_this$prob_gametes_diploid_0[row_id_this]
  A_this[2,2] = df_this$prob_gametes_diploid_1[row_id_this]
  A_this[2,3] = df_this$prob_gametes_diploid_2[row_id_this]
  A_this[3,1] = df_this$prob_gametes_triploid_0[row_id_this]
  A_this[3,2] = df_this$prob_gametes_triploid_1[row_id_this]
  A_this[3,3] = df_this$prob_gametes_triploid_2[row_id_this]
  A_this[3,4] = df_this$prob_gametes_triploid_3[row_id_this]
  A_this[4,1] = df_this$prob_gametes_tetraploid_0[row_id_this]
  A_this[4,2] = df_this$prob_gametes_tetraploid_1[row_id_this]
  A_this[4,3] = df_this$prob_gametes_tetraploid_2[row_id_this]
  A_this[4,4] = df_this$prob_gametes_tetraploid_3[row_id_this]
  A_this[4,5] = df_this$prob_gametes_tetraploid_4[row_id_this]

  B_this = matrix(nrow=1,ncol=4,data=0)
  B_this[1,1] = df_this$prob_apomixis_unreduced_gametes_haploid[row_id_this]
  B_this[1,2] = df_this$prob_apomixis_unreduced_gametes_diploid[row_id_this]
  B_this[1,3] = df_this$prob_apomixis_unreduced_gametes_triploid[row_id_this]
  B_this[1,4] = df_this$prob_apomixis_unreduced_gametes_tetraploid[row_id_this]
  
  C_this = matrix(nrow=1,ncol=4,data=0)
  C_this[1,1] = df_this$prob_survival_offspring_haploid[row_id_this]
  C_this[1,2] = df_this$prob_survival_offspring_diploid[row_id_this]
  C_this[1,3] = df_this$prob_survival_offspring_triploid[row_id_this]
  C_this[1,4] = df_this$prob_survival_offspring_tetraploid[row_id_this]
  
  g_A_this = plot_matrix(A_this,subtract_one = TRUE) + ggtitle('Gamete probabilities (A)')
  g_B_this = plot_matrix(B_this) + ggtitle('Apomixis probabilities (B)')
  g_C_this = plot_matrix(C_this) + ggtitle('Survival probabilities (C)')
  
  g_final = ggarrange(g_A_this, g_B_this, g_C_this,
            nrow=3,
            heights=c(2,1,1))
  
  return(g_final)
}

params_in_target_0 = plot_params(36,df_all)
params_in_target_1 = plot_params(623,df_all)

g_ts_example = ggarrange(params_in_target_1, g_ts_in_target_1,params_in_target_0, g_ts_in_target_0, 
                         labels=c("(a)","","(b)",""),
                         widths=c(2,3))
ggsave(g_ts_example, file=sprintf('%s/g_ts_example.pdf',dir_results),width=10,height=10)
ggsave(g_ts_example, file=sprintf('%s/g_ts_example.png',dir_results),width=10,height=10)

# df_example_params = df_all %>% 
#   mutate(row=1:n()) %>%
#   group_by(dir_out) %>% 
#   slice_head(n=1)
# 
# for (i in 1:nrow(df_example_params))
# {
#   ggsave(plot_params(i,df_example_params),file=sprintf('%s/g_params_row=%d_%s.pdf',dir_results, df_example_params$row[i], df_example_params$dir_out[i]),width=5,height=10)
# }


# 
# 
# # pca doesn't work very well
# pca = prcomp(df_all %>% select(all_of(names_predictors)),center = TRUE, scale=TRUE)
# df_pca = data.frame(in_target=df_quantitative$in_target, pca$x[,1:3], dir_out=df_all$dir_out)
# 
# g_pca_12 = ggplot(df_pca, aes(x=PC1, y=PC2,color=in_target)) +
#   geom_point(size=0.2,alpha=0.1) +
#   theme_bw() +
#   scale_color_manual(values=c('gray','blue')) +
#   facet_wrap(~dir_out,nrow=1)
# 
# g_pca_23 = ggplot(df_pca, aes(x=PC2, y=PC3,color=in_target)) +
#   geom_point(size=0.2,alpha=0.1) +
#   theme_bw() +
#   scale_color_manual(values=c('gray','blue')) +
#   facet_wrap(~dir_out,nrow=1)
# 
# ggsave(ggarrange(g_pca_12, g_pca_23,align='hv',nrow=2, common.legend = TRUE,legend='bottom'),file=sprintf('%s/g_pca.pdf',dir_results),width=14,height=8)
# 
# #plot3d(data_pca$PC1, data_pca$PC2, data_pca$PC3,col=rainbow(15)[cut(data_pca$distance,breaks=10)])
# 
# #biplot(pca)






# rename and split data
df_all_split = df_all %>% 
  group_by(dir_out) %>%
  select(dir_out, in_target, distance, all_of(names_predictors))
# split into pieces
df_all_split = df_all_split %>% 
  group_split()

# fit random forests
rf_list = lapply(X=df_all_split, FUN=function(df_this) {
  cat('.')
  
  freqs = table(df_this$in_target)
  
  m_rf = ranger(formula(sprintf("distance ~ %s",paste(names_predictors,collapse=" + "))), 
                data=df_this,
                mtry=5,
                num.trees=1000,
                verbose=TRUE,
                #classification=TRUE,
                importance='impurity',
                #case.weights=1/df_this$distance
                case.weights=ifelse(df_this$in_target==1,freqs["FALSE"]/freqs["TRUE"],1)
                )
  return(m_rf)
})
names(rf_list) = sapply(df_all_split, function(x) { x$dir_out[1] })

# post-hoc explanation of relationship among variables
df_importance = cbind(dir_out=unique(df_all$dir_out), do.call("rbind",lapply(rf_list, function(x) { data.frame(t(importance(x))) }))) %>%
  pivot_longer(!dir_out) %>%
  group_by(dir_out) %>%
  mutate(value_relative = value/max(value,na.rm=T)) %>%
  mutate(nice_name = nice_names[name])

g_rf_importance = ggplot(df_importance, aes(x=nice_name,y=value,fill=value_relative)) +
  geom_bar(stat='identity') +
  facet_wrap(~dir_out,scales='free',ncol=2) +
  theme_bw() +
  coord_flip() +
  scale_fill_viridis_c(direction=-1) + 
  scale_x_discrete(limits = rev(sort(unique(df_importance$nice_name)))) +
  xlab('Importance (Gini impurity)') +
  ylab('Parameter') +
  theme(legend.position='none')
ggsave(g_rf_importance, file=sprintf('%s/g_rf_importance.pdf',dir_results),width=8,height=12)
ggsave(g_rf_importance, file=sprintf('%s/g_rf_importance.png',dir_results),width=8,height=12)




# plot params

# all
g_scenario_all_1 = df_all %>%
  filter(dir_out=='all') %>%
  ggplot(aes(x=prob_survival_offspring_triploid/prob_survival_offspring_diploid,y=prob_apomixis_unreduced_gametes_triploid/prob_apomixis_unreduced_gametes_diploid,color=distance)) +
  geom_point(size=0.3,alpha=0.5) +
  theme_bw() +
  scale_x_log10(limits=c(1e-2,1e2)) +
  scale_y_log10(limits=c(1e-2,1e2)) +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1),name='Distance (d)') +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("C3/C2") + ylab("B3/B2")

g_scenario_all_2 = df_all %>%
  filter(dir_out=='all') %>%
  ggplot(aes(x=prob_gametes_diploid_2,y=prob_gametes_diploid_1,color=distance)) +
  geom_point(size=0.3,alpha=0.5) +
  theme_bw() +
  xlim(0,1) + ylim(0,1) +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1),name='Distance (d)') +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("A2,2") + ylab("A2,1")

ggsave(ggarrange(g_scenario_all_1, g_scenario_all_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels='auto',
                 legend='bottom'),
       file=sprintf('%s/g_scenario_all.pdf',dir_results), 
       width=12,height=7) 
ggsave(ggarrange(g_scenario_all_1, g_scenario_all_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels='auto',
                 legend='bottom'),
       file=sprintf('%s/g_scenario_all.png',dir_results), 
       width=12,height=7) 


# no apomixis
g_scenario_no_apomixis_1 = df_all %>%
  filter(dir_out=='no_apomixis') %>%
  ggplot(aes(x=prob_survival_offspring_triploid/prob_survival_offspring_diploid,y=prob_survival_offspring_tetraploid/prob_survival_offspring_diploid,color=distance)) +
  geom_point(size=0.3) +
  theme_bw() +
  scale_x_log10(limits=c(1e-2,1e2)) +
  scale_y_log10(limits=c(1e-2,1e2)) +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1),name='Distance (d)') +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("C3/C2") + ylab("C4/C1")

g_scenario_no_apomixis_2 = df_all %>%
  filter(dir_out=='no_apomixis') %>%
  ggplot(aes(x=prob_gametes_diploid_1,y=prob_gametes_diploid_2,color=distance)) +
  geom_point(size=0.3,alpha=0.5) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1),name='Distance (d)') +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("A2,2") + ylab("A2,1")

ggsave(ggarrange(g_scenario_no_apomixis_1, g_scenario_no_apomixis_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels='auto',
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_apomixis.pdf',dir_results), 
       width=12,height=7) 
ggsave(ggarrange(g_scenario_no_apomixis_1, g_scenario_no_apomixis_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels='auto',
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_apomixis.png',dir_results), 
       width=12,height=7) 

# no haploid or tetraploid
g_scenario_no_haploid_tetraploid_1 = df_all %>%
  filter(dir_out=='no_haploid_tetraploid') %>%
  ggplot(aes(x=prob_survival_offspring_triploid/prob_survival_offspring_diploid,y=prob_apomixis_unreduced_gametes_triploid/prob_apomixis_unreduced_gametes_diploid,color=distance)) +
  geom_point(size=0.3) +
  theme_bw() +
  scale_x_log10(limits=c(1e-2,1e2)) +
  scale_y_log10(limits=c(1e-2,1e2)) +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("C3/C2") + ylab("B3/B2")

g_scenario_no_haploid_tetraploid_2 = df_all %>%
  filter(dir_out=='no_haploid_tetraploid') %>%
  ggplot(aes(x=prob_gametes_diploid_2,y=prob_gametes_diploid_1,color=distance)) +
  geom_point(size=0.3) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("A2,2") + ylab("A2,1")

ggsave(ggarrange(g_scenario_no_haploid_tetraploid_1, g_scenario_no_haploid_tetraploid_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels = 'auto',
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_haploid_tetraploid.pdf',dir_results), 
       width=12,height=7) 
ggsave(ggarrange(g_scenario_no_haploid_tetraploid_1, g_scenario_no_haploid_tetraploid_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels = 'auto',
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_haploid_tetraploid.png',dir_results), 
       width=12,height=7) 

# no survival variation
g_scenario_no_survival_variation_1 = df_all %>%
  filter(dir_out=='no_survival_variation') %>%
  ggplot(aes(x=prob_apomixis_unreduced_gametes_triploid/prob_apomixis_unreduced_gametes_diploid,y=prob_apomixis_unreduced_gametes_tetraploid/prob_apomixis_unreduced_gametes_haploid,color=distance)) +
  geom_point(size=0.3) +
  theme_bw() +
  scale_x_log10(limits=c(1e-2,1e2)) +
  scale_y_log10(limits=c(1e-2,1e2)) +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("B3/B2") + ylab("B4/B1")

g_scenario_no_survival_variation_2 = df_all %>%
  filter(dir_out=='no_survival_variation') %>%
  ggplot(aes(x=prob_gametes_diploid_2,y=prob_gametes_diploid_1,color=distance)) +
  geom_point(size=0.3) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("A2,2") + ylab("A2,1")

ggsave(ggarrange(g_scenario_no_survival_variation_1, g_scenario_no_survival_variation_2, 
                 align='hv',
                 common.legend = TRUE,
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_survival_variation.pdf',dir_results), 
       width=12,height=7) 
ggsave(ggarrange(g_scenario_no_survival_variation_1, g_scenario_no_survival_variation_2, 
                 align='hv',
                 common.legend = TRUE,
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_survival_variation.png',dir_results), 
       width=12,height=7) 

# no triploid fertility
g_scenario_no_triploid_fertility_1 = df_all %>%
  filter(dir_out=='no_triploid_fertility') %>%
  ggplot(aes(x=prob_survival_offspring_triploid/prob_survival_offspring_diploid,y=prob_apomixis_unreduced_gametes_triploid/prob_apomixis_unreduced_gametes_diploid,color=distance)) +
  geom_point(size=0.3,alpha=0.5) +
  theme_bw() +
  scale_x_log10(limits=c(1e-2,1e2)) +
  scale_y_log10(limits=c(1e-2,1e2)) +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1),name='Distance (d)') +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("C3/C2") + ylab("B3/B2")

g_scenario_no_triploid_fertility_2 = df_all %>%
  filter(dir_out=='no_triploid_fertility') %>%
  ggplot(aes(x=prob_gametes_diploid_2,y=prob_gametes_diploid_1,color=distance)) +
  geom_point(size=0.3) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("A2,2") + ylab("A2,1")

ggsave(ggarrange(g_scenario_no_triploid_fertility_1, g_scenario_no_triploid_fertility_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels='auto',
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_triploid_fertility.pdf',dir_results), 
       width=12,height=7) 
ggsave(ggarrange(g_scenario_no_triploid_fertility_1, g_scenario_no_triploid_fertility_2, 
                 align='hv',
                 common.legend = TRUE,
                 labels='auto',
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_triploid_fertility.png',dir_results), 
       width=12,height=7) 



# no survival variation no apomixis (just gametes
g_scenario_no_survival_variation_no_apomixis_1 = df_all %>%
  filter(dir_out=='no_survival_variation_no_apomixis') %>%
  ggplot(aes(x=prob_gametes_diploid_2/prob_gametes_diploid_1,y=prob_gametes_tetraploid_4/prob_gametes_tetraploid_2,color=distance)) +
  geom_point(size=0.3) +
  theme_bw() +
  scale_x_log10(limits=c(1e-2,1e2)) +
  scale_y_log10(limits=c(1e-2,1e2)) +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("A2,2/A2,1") + ylab("A4,4/A4,2")

g_scenario_no_survival_variation_no_apomixis_2 = df_all %>%
  filter(dir_out=='no_survival_variation_no_apomixis') %>%
  ggplot(aes(x=prob_gametes_triploid_3/prob_gametes_diploid_1,y=prob_gametes_triploid_1/prob_gametes_diploid_1,color=distance)) +
  geom_point(size=0.3) +
  theme_bw() +
  scale_x_log10(limits=c(1e-2,1e2)) +
  scale_y_log10(limits=c(1e-2,1e2)) +
  scale_color_gradientn(colors=c('purple','red','orange','yellow','gray'),values=c(0,0.05,0.1,0.5,1)) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  xlab("A3,3/A2,1") + ylab("A3,1/A2,1")

ggsave(ggarrange(g_scenario_no_survival_variation_no_apomixis_1, g_scenario_no_survival_variation_no_apomixis_2, 
                 align='hv',
                 common.legend = TRUE,
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_survival_variation_no_apomixis.pdf',dir_results), 
       width=12,height=7) 
ggsave(ggarrange(g_scenario_no_survival_variation_no_apomixis_1, g_scenario_no_survival_variation_no_apomixis_2, 
                 align='hv',
                 common.legend = TRUE,
                 legend='bottom'),
       file=sprintf('%s/g_scenario_no_survival_variation_no_apomixis.png',dir_results), 
       width=12,height=7) 





# do PDPs
do_pdp_2d <- function(m_rf, pred_x, pred_y, grid_resolution=5)
{
  grid = expand.grid(x=seq(0,1,length.out=grid_resolution),y=seq(0,1,length.out=grid_resolution))
  names(grid) = c(pred_x, pred_y)
  
  pdp_this = partial(m_rf,
              pred.var=c(pred_x,pred_y),
              pred.grid=grid,
              train=df_all %>% filter(dir_out=='all') %>% sample_n(1000),
              progress = TRUE,
              plot = FALSE)
  
  g_this = ggplot(pdp_this, aes_string(x=pred_x,y=pred_y,fill="yhat")) + 
    geom_raster() +
    scale_fill_viridis_c(name='Distance (d)') +
    theme_bw() +
    xlab(nice_names[pred_x]) +
    ylab(nice_names[pred_y])
  
  return(g_this)
}

do_pdp_1d <- function(m_rf, pred_x, grid_resolution=5)
{
  grid = expand.grid(x=seq(0,1,length.out=grid_resolution))
  names(grid) = c(pred_x)
  
  pdp_this = partial(m_rf,
                     pred.var=c(pred_x),
                     pred.grid=grid,
                     train=df_all %>% filter(dir_out=='all') %>% sample_n(1000),
                     progress = TRUE,
                     plot = FALSE)
  
  # g_this = ggplot(pdp_this, aes_string(x=pred_x,y="yhat")) + 
  #   geom_point() +
  #   geom_line() +
  #   theme_bw() +
  #   xlab(nice_names[pred_x]) +
  #   ylab("Distance (d)")
  
  return(pdp_this)
}

pdp_gametes_diploid_12 = do_pdp(m_rf = rf_list[["all"]], 
       pred_x = "prob_gametes_diploid_1", 
       pred_y = "prob_gametes_diploid_2",
       grid_resolution=7)
pdp_gametes_triploid_13 = do_pdp(m_rf = rf_list[["all"]], 
                             pred_x = "prob_gametes_triploid_1", 
                             pred_y = "prob_gametes_triploid_3",
                             grid_resolution=7)
pdp_gametes_tetraploid_42 = do_pdp(m_rf = rf_list[["all"]], 
                                 pred_x = "prob_gametes_tetraploid_2", 
                                 pred_y = "prob_gametes_tetraploid_4",
                                 grid_resolution=7)
pdp_gametes_tetraploid_diploid = do_pdp(m_rf = rf_list[["all"]], 
                                   pred_x = "prob_gametes_tetraploid_2", 
                                   pred_y = "prob_gametes_diploid_1",
                                   grid_resolution=7)

ggarrange(pdp_gametes_diploid_12, pdp_gametes_triploid_13, pdp_gametes_tetraploid_42, pdp_gametes_tetraploid_diploid,
          common.legend = TRUE)


do_pdp_1d(m_rf = rf_list[["all"]], 
          pred_x = "prob_gametes_triploid_1",
          grid_resolution = 10)


names_gametes = df_all %>% 
  select(contains("prob_gametes")) %>% 
  names

pdps_1d = lapply(names_gametes, do_pdp_1d,
       m_rf = rf_list[["all"]], 
       grid_resolution = 20)

pdps_1d_joined = rbindlist(lapply(pdps_1d, function(x) { 
  df = data.frame(x,var=nice_names[names(x)[1]])
  names(df) = c("x","yhat","var")
  return(df) }))
pdps_1d_joined = pdps_1d_joined %>%
  group_by(var) %>%
  mutate(is.min=yhat==min(yhat))

g_pdps_1d = ggplot(pdps_1d_joined, aes(x=x,y=yhat,color=is.min,group=var)) +
  theme_bw() +
  geom_point() +
  geom_line(color='black') +
  facet_wrap(~var) +
  scale_color_manual(values=c('black','red')) +
  theme(legend.position='none') +
  ylab("Distance (d)") +
  xlab("Value")
ggsave(g_pdps_1d,file=sprintf('%s/g_pdps_1d.png',dir_results), 
       width=8,height=8) 

# full range
# haven't figured out clear parameters



# try filtering simulation outputs to only the realistic cases
# maybe more efficient to pre-set the parameters to this restricted range via sampling
# data_realistic = df_all %>%
#           filter(dir_out=='outputs_all_2023-10-09') %>%
#           filter(in_target==TRUE)
# 					filter(prob_gametes_haploid_1 < 0.2)
# 					
# nrow(data_realistic)/nrow(data)
# data_realistic$distance




ggplot(df_all %>% filter(dir_out=='outputs_restricted_range_2023-10-09'), aes(x=prob_gametes_diploid_2,y=prob_gametes_tetraploid_3,color=distance)) + geom_point(size=0.5) +
  scale_color_gradientn(colors=c('purple','red','orange','lightgray'),values=c(0,0.01,0.1,1))

