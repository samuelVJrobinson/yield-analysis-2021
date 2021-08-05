#Compare filtered + unfiltered models

library(tidyverse)
library(ggpubr)
theme_set(theme_classic())

# Canola models ------------------------

# croptype <- 'wheat'
croptype <- 'canola'

path1 <- paste0('./Data/postSamples_',croptype,'.Rdata')
path2 <- paste0('/media/rsamuel/Storage/other/posteriorSamples_older/postSamples_',croptype,'.Rdata')

#Filtered
load(path1)
samp_filt <- samp
#Unfiltered
load(path2)
samp_unfilt <- samp
rm(samp); gc()

ylabMean <- 'Mean Yield (T/ha)'
ylabSD <- 'SD Yield'
alphaVal <- 0.1 #Transparancy 
qs <- c(0.1,0.5,0.9) #Quantiles
scalecols <- c('blue','red')

pAreaSamp <- lapply(list(samp_filt,samp_unfilt),function(samp){
  lapply(samp,function(x) x$pArea) %>% bind_rows(.id = 'sample') %>% 
    group_by(pArea) %>% 
    summarise(predMean = quantile(predMean, qs), predLogSD = exp(quantile(predLogSD, qs)), q = qs) %>% 
    ungroup() %>% 
    mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
    pivot_wider(names_from=q,values_from=c(predMean,predLogSD)) %>% 
    data.frame}) %>% 
  set_names('Filtered','Unfiltered') %>% bind_rows(.id='Data')

p1 <- ggplot(pAreaSamp,aes(x=pArea))+
  geom_ribbon(aes(ymax=predMean_upr,ymin=predMean_lwr,fill=Data),alpha=0.3)+
  geom_line(aes(y=predMean_med,col=Data),size=1)+
  labs(x='Polygon Area',y=ylabMean)+
  scale_color_manual(values=scalecols)+scale_fill_manual(values=scalecols)

p2 <- ggplot(pAreaSamp,aes(x=pArea))+
  geom_ribbon(aes(ymax=predLogSD_upr,ymin=predLogSD_lwr,fill=Data),alpha=0.3)+
  geom_line(aes(y=predLogSD_med,col=Data),size=1)+
  labs(x='Polygon Area',y=ylabSD)+
  scale_color_manual(values=scalecols)+
  scale_fill_manual(values=scalecols)+
  coord_cartesian(ylim=c(0,1))

coverSamp <- lapply(list(samp_filt,samp_unfilt),function(samp){
  lapply(samp,function(x) x$coverDist %>% bind_rows(.id = 'dist_type')) %>% 
    bind_rows(.id = 'sample') %>% group_by(dist_type,dist) %>% 
    summarise(predMean = quantile(predMean, qs), predLogSD = exp(quantile(predLogSD, qs)), q = qs) %>% 
    ungroup() %>% 
    mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
    pivot_wider(names_from=q,values_from=c(predMean,predLogSD)) %>% 
    data.frame}) %>% 
  set_names('Filtered','Unfiltered') %>% bind_rows(.id='Data')

p3 <- coverSamp %>% filter(dist<200) %>% 
  ggplot(aes(x=dist)) + 
  geom_ribbon(aes(ymax=predMean_upr,ymin=predMean_lwr,fill=Data),alpha=0.3) + 
  geom_line(aes(y=predMean_med,col=Data)) + #Meta-model
  geom_hline(yintercept = 0,col='black',linetype='dashed')+
  facet_wrap(~dist_type,scales='free') + labs(x='Distance',y=ylabMean) +
  scale_color_manual(values=scalecols)+
  scale_fill_manual(values=scalecols)

p4 <- coverSamp %>% filter(dist<200) %>% 
  ggplot(aes(x=dist)) + 
  geom_ribbon(aes(ymax=predLogSD_upr,ymin=predLogSD_lwr,fill=Data),alpha=0.3) + 
  geom_line(aes(y=predLogSD_med,col=Data)) + #Meta-model
  facet_wrap(~dist_type,scales='free') + 
  labs(x='Distance',y=ylabSD) +
  scale_color_manual(values=scalecols)+
  scale_fill_manual(values=scalecols)+
  coord_cartesian(ylim=c(0,5))


rSamp <- lapply(list(samp_filt,samp_unfilt),function(samp){
  lapply(samp,function(x) x$r) %>% bind_rows(.id = 'sample') %>% 
    group_by(r) %>% 
    summarise(predMean = quantile(predMean, qs), predLogSD = exp(quantile(predLogSD, qs)), q = qs) %>% 
    ungroup() %>% 
    mutate(q=factor(q,labels=c('lwr','med','upr'))) %>% 
    pivot_wider(names_from=q,values_from=c(predMean,predLogSD)) %>% 
    data.frame}) %>% 
  set_names('Filtered','Unfiltered') %>% bind_rows(.id='Data')

p5 <- ggplot(rSamp,aes(x=r))+
  geom_ribbon(aes(ymax=predMean_upr,ymin=predMean_lwr,fill=Data),alpha=0.3)+
  geom_line(aes(y=predMean_med,col=Data),size=1)+
  geom_hline(yintercept = 0,col='black',linetype='dashed')+
  labs(x='Harvest sequence',y=ylabMean)+
  scale_color_manual(values=scalecols)+scale_fill_manual(values=scalecols)

p6 <- ggplot(rSamp,aes(x=r))+
  geom_ribbon(aes(ymax=predLogSD_upr,ymin=predLogSD_lwr,fill=Data),alpha=0.3)+
  geom_line(aes(y=predLogSD_med,col=Data),size=1)+
  labs(x='Harvest sequence',y=ylabSD)+
  scale_color_manual(values=scalecols)+scale_fill_manual(values=scalecols) +
  coord_cartesian(ylim=c(0,5))

(p <- ggarrange(p1,p3,p5,p2,p4,p6,ncol=3,nrow=2,common.legend = TRUE, legend = 'bottom'))
ggsave(paste0('./Figures/compareFiltered_',croptype,'.png'),p,height=6,width=16,dpi=350)

(p <- ggarrange(p3+facet_wrap(~dist_type,nrow=1),p4+facet_wrap(~dist_type,nrow=1),ncol=1,nrow=2,common.legend = TRUE, legend = 'bottom'))
ggsave(paste0('./Figures/compareFiltered_',croptype,'_b.png'),p,height=6,width=16,dpi=350)
