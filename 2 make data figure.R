library(dplyr)
library(terra)
library(ggplot2)
library(basemaps)
library(sf)
library(terra)
library(RStoolbox)
library(ggspatial)
library(ggpubr)

data = read.csv('/Users/benjaminblonder/Documents/ASU/aspen remote sensing/2019/data analysis 2020/aspen data site-level processed 30 Mar 2020.csv')

data_grid = data %>%
  filter(Point_Type=='Grid')#%>%
  #filter(Watershed=='Middle East River')
    
cytotype_counts_by_genotype = data_grid %>%
  group_by(Genotype) %>%
  slice_head(n=1) %>%
  ungroup %>%
  group_by(Watershed,Ploidy_level) %>%
  tally() %>%
  mutate(Ploidy_level=ifelse(is.na(Ploidy_level),'Unknown',Ploidy_level)) %>%
  mutate(Cytotype=factor(Ploidy_level,levels=c('Haploid','Diploid','Triploid','Tetraploid','Unknown'),labels=c('1','2','3','4','Unknown'),ordered=TRUE))

ggplot(cytotype_counts_by_genotype, aes(x=Watershed,y=n,fill=Cytotype)) +
  theme_bw() +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('black','blue','red','orange','gray'),drop=FALSE) +
  ylab('Count') +
  xlab('Grid name') +
  coord_flip()


plot_watershed <- function(ws)
{
  points_this = data_grid %>%
    filter(Watershed==ws) %>%
    select(x=X.UTM,y=Y.UTM,Ploidy_level) %>%
    st_as_sf(coords=c("x","y"),crs="+proj=utm +zone=13 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") %>%
    mutate(Ploidy_level=ifelse(is.na(Ploidy_level),'Unknown',Ploidy_level)) %>%
    mutate(Cytotype=factor(Ploidy_level,levels=c('Haploid','Diploid','Triploid','Tetraploid','Unknown'),labels=c('1','2','3','4','Unknown'),ordered=TRUE))
  
  map_this = basemap_raster(ext = st_bbox(points_this), map_service = "mapbox", map_type = "satellite", map_token='pk.eyJ1IjoiYmVuamFtaW5ibG9uZGVyYmVya2VsZXkiLCJhIjoiY2xudGd3MWc3MDNvczJrcXh4dmQzdjVyMCJ9.4LCZn7T11Fj4EkOINtPMvg')
  
  points_this_projected = st_transform(points_this, crs=crs(map_this))
  
  g = ggRGB(map_this,r=1,g=2,b=3,scale=FALSE) +
    geom_sf(data=points_this_projected,aes(color=Cytotype),size=3) +
    scale_color_manual(values=c('black','blue','red','orange','gray'),drop=FALSE) +
    theme_void() +
    annotation_scale(location = "tr", width_hint =0.5, pad_x = unit(0.7, "cm"))
  
  return(g)
}

plots_all_list = lapply(unique(data_grid$Watershed), plot_watershed)

if(!file.exists('stats'))
{
  dir.create('stats')
}
ggsave(ggarrange(plotlist=plots_all_list,labels='auto',common.legend = TRUE,legend='bottom'),
       file='stats/g_maps.pdf',width=10,height=10)
ggsave(ggarrange(plotlist=plots_all_list,labels='auto',common.legend = TRUE,legend='bottom'),
       file='stats/g_maps.png',width=10,height=10)

unique(data_grid$Watershed)
