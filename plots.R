sp_list <- read_csv('species_list.csv') 
# get shorter names to label the species
newsplabs <- data.frame(species=rownames(nmds_spp$species)) %>% left_join(sp_list) %>% pull(plot)
#NMDS
ordiplot(nmds_spp, type = 'none', las = 1)
ordiellipse(nmds_spp, groups = community_matrix_sites$class2, draw = 'polygon', col = c('#3b674b','#8ab07e','#d8cb39'), alpha = 0.5)
# ordiellipse(nmds_spp, groups = site_type_2, draw = 'polygon', border = FALSE, col = c('forestgreen','gold','firebrick'), alpha = 0.4)
points(nmds_spp, display = 'sites', select = community_matrix_sites$class2=='Forest', col = '#3b674b', pch = 16)
points(nmds_spp, display = 'sites', select = community_matrix_sites$class2=='Reforestation', col = '#8ab07e', pch = 16)
points(nmds_spp, display = 'sites', select = community_matrix_sites$class2=='Natural regeneration', col = '#d8cb39', pch = 16)

text(nmds_spp, display = 'species', labels = newsplabs)
pnmds <- recordPlot()

# Simpson diversity
tibble(sim=simpson_index, class = factor(site_type_2, levels = c('Forest','Reforestation','Natural regeneration'))) %>% 
  ggplot()+geom_boxplot(aes(class, sim, fill = class), show.legend = FALSE)+
  labs(x='Site type',y='Diversity')+
  # scale_fill_manual(values = c('#134e13','#228b22','#31c831'))+
  scale_fill_manual(values = c('#3b674b','#8ab07e','#d8cb39'))+
  theme_classic(base_size = 18) -> psim

# Rarefaction index
db_rar %>% mutate(class = factor(class, levels = c('Forest','Reforestation','Natural regeneration'))) %>% 
  ggplot()+geom_boxplot(aes(class, ri, fill=class), show.legend = FALSE)+
  labs(x='Site type',y='Rarefaction index')+
  # scale_fill_manual(values = c('#134e13','#31c831','#228b22'))+
  scale_fill_manual(values = c('#3b674b','#8ab07e','#d8cb39'))+
  theme_classic(base_size = 18) -> prar

cowplot::plot_grid(prar,psim)
# CCA
cca_plot <- cca(cca_spp~PC1+PC2, data=prin_comps)
plot(cca_plot, type = 'none', display = c('sp'), cex.axis=1.5, cex.lab=1.5, las=1)
text(cca_plot, display = 'bp', col = 'firebrick', labels = c('Vegetation', 'Tree cover'))
text(cca_plot, display = 'species', cex = 1.2)

# forest plot glm of detections
ggplot(glm_freq_pred)+
  geom_pointrange(aes(x = fct_reorder(Species,.x = fit), y = fit, ymin = fit-se.fit, ymax = fit+se.fit, color = class2), position = position_dodge(width = 0.5))+
  scale_color_manual(values = c('#3b674b','#8ab07e','#d8cb39'))+
  theme_classic(base_size = 16)+
  labs(x = 'Species', y = 'Detections', color = 'Site class')+
  coord_flip()

# forest plot, raw detection frequency data
sp3cat <- db_freq %>% distinct(class2, Species) %>% count(Species) %>% filter(n==3)
db_freq %>% slice(-148) %>% 
  mutate(Species = factor(Species, levels = sp_list$species, labels = sp_list$Common)) %>% 
  # filter(Species %in% sp3cat$Species) %>% 
  # group_by(class2, Species) %>% 
  group_by(cavity_type, Species) %>% 
  summarise(meanfreq=mean(freq), se = sd(freq)/sqrt(n())) %>% 
  ggplot(aes(fct_reorder(Species,meanfreq, mean), 7*meanfreq))+geom_pointrange(aes(ymin = 7*(meanfreq-se), ymax = 7*(meanfreq+se), color=class2), position = position_dodge(width = 0.5))+
  scale_color_manual(values = c('#3b674b','#8ab07e','#d8cb39'))+
  theme_classic(base_size = 16)+
  labs(x = '', y = 'Detection frequency (per week)', color = 'Site class')+
  coord_flip()

# forest plot, overall frequency + frequency of interaction
db_behav  %>% mutate(interaxn = entrance | interaction | mark | eat | exit) %>% 
  left_join(select(station_data, Station,work_time)) %>% 
  filter(!(Station == "N_4" & Species =="coati")) %>% 
  group_by(Station, Species) %>% mutate(inter_freq=n()/work_time) %>% 
  distinct(Species,cavity_type,inter_freq) %>% 
  group_by(Species, cavity_type) %>% summarise(mean_freq = mean(inter_freq), se = sd(inter_freq/sqrt(n()))) %>% 
  ggplot()+
  geom_pointrange(aes(y=fct_reorder(Species,mean_freq,max), x=7*mean_freq, xmin = 7*(mean_freq-se), xmax=7*(mean_freq+se), color = cavity_type), show.legend = FALSE, position = position_dodge(width = 0.5))+
  theme_minimal()+
  scale_x_reverse()+
  geom_vline(xintercept = 0)+
  scale_color_manual(values = c("#3b674b","#d8cb39"))+
  labs(x = "Detection frequency (per week)",
       y = "",
       color = "Cavity type")

# combine the two
# get order of species according to max detection frequency
sp_order_forestplot <- db_behav %>% count(Station, Species) %>% 
  left_join(select(station_data, Station, work_time)) %>% 
  mutate(freq = n/work_time) %>% 
  group_by(Species) %>% summarise(m = mean(freq)) %>% 
  arrange(m) %>% select(Species) %>% left_join(sp_list, by = c('Species'='species'))

freq_fplot <- db_behav %>% 
  filter(!(Station == "N_4" & Species =="coati")) %>% 
  count(Station, Species, cavity_type) %>% 
  left_join(select(station_data, Station, work_time)) %>% 
  mutate(freq = n/work_time) %>% 
  group_by(Species, cavity_type) %>% summarise(m = mean(freq), se = sd(freq/sqrt(n())))
%>% 
  mutate(Species = factor(Species, levels = sp_order_forestplot$Species))

ifreq_fplot <- db_behav  %>% mutate(interaxn = entrance | interaction | mark | eat | exit) %>% 
  filter(!(Station == "N_4" & Species =="coati"), interaxn) %>% 
  count(Station, Species, cavity_type) %>% 
  left_join(select(station_data, Station, work_time)) %>% 
  mutate(freq = n/work_time) %>% 
  group_by(Species, cavity_type) %>% summarise(m = mean(freq), se = sd(freq/sqrt(n()))) 
%>% 
  mutate(Species = factor(Species, levels = sp_order_forestplot$Species))
#### base r forest plot ####
baseR_freq_data <- left_join(sp_order_forestplot, freq_fplot) %>% left_join(ifreq_fplot, by=c('Species', 'cavity_type')) %>% 
  pivot_wider(names_from = cavity_type, values_from = m.x:se.y)

plot.new()
par(mar = c(5,10,2,2)+0.5, cex = 1.1)
plot.window(xlim = 7*c(-0.1,0.1), ylim = c(0,20))
axis(side = 2, at = 0:20, labels = baseR_freq_data$Common, las=1, tick = FALSE)
axis(1, at = seq(-0.6,0.6,0.2,), labels = c(seq(0.6,0,-0.2), seq(0.2,0.6,0.2)))
title(xlab="Frequency (per week)")
abline(h = 0:20, lty=2, col='gray80')
abline(v = 0)
mtext(text = c("Natural", "Artificial"), side = 3, at = c(-0.4, 0.4), cex = 1.3)
rect(xleft = 0, ybottom = 0:20-0.4, xright = 7*baseR_freq_data$m.x_Artificial, ytop = (0:20)+0.4, col = '#d8cb39')
rect(xleft = 0, ybottom = 0:20-0.4, xright = 7*baseR_freq_data$m.y_Artificial, ytop = (0:20)+0.4, col = '#e8e18b')
rect(xleft = 0, ybottom = 0:20-0.4, xright = -7*baseR_freq_data$m.x_Natural, ytop = (0:20)+0.4, col = '#3b674b')
rect(xleft = 0, ybottom = 0:20-0.4, xright = -7*baseR_freq_data$m.y_Natural, ytop = (0:20)+0.4, col = '#60a479')

# order species according to guild
sp_order_forestplot2 <- data.frame(Species=c('paca','agouti','rodent','squirrel_uid',# small herbivores
                          'tapir','peccary_whitelipped','peccary_collared',#large herb
                          'coati','raccoon','opossum_common','opossum_foureyed',# omnivores
                          'tamandua','armadillo_ninebanded','skunk_striped',#insectivores
                          'puma','ocelot','margay','tayra',# predators
                          'capuchin_whitefaced','squirrel_monkey',#primates
                          'curassow_great'#bird
                          ))
                          
baseR_freq_data2 <- read.csv('freq_interact_data.csv') %>% 
  left_join(sp_list, by = c('Species'='species')) %>% 
  mutate(ifreq = freq*interact/100) %>% 
  select(cavity_type,Species,Latin, Common, plot, freq,ifreq) %>% 
  pivot_wider(names_from = cavity_type, values_from = c(freq,ifreq))

baseR_freq_data2 <- left_join(sp_order_forestplot2,baseR_freq_data2)
                          
plot.new()
par(mar = c(5,10.5,2,2)+0.5, cex = 1.1)
plot.window(xlim = 5*c(-0.1,0.1), ylim = c(0,20))
axis(side = 2, at = 0:20, labels = baseR_freq_data2$Common, las=1, tick = FALSE)
axis(1,at = seq(-0.4,0.4,0.2), labels = c(seq(0.4,0,-0.2), seq(0.2,0.4,0.2)))
title(xlab="Frequency (per week)")
abline(h = 0:20, lty=2, col='gray80')
abline(v = 0)
mtext(text = c("Natural", "Artificial"), side = 3, at = c(-0.3, 0.3), cex = 1.3)
rect(xleft = 0, ybottom = 0:20-0.4, xright = 7*baseR_freq_data2$freq_Artificial, ytop = (0:20)+0.4, col = '#d8cb39')
rect(xleft = 0, ybottom = 0:20-0.4, xright = 7*baseR_freq_data2$ifreq_Artificial, ytop = (0:20)+0.4, col = '#e8e18b')
rect(xleft = 0, ybottom = 0:20-0.4, xright = -7*baseR_freq_data2$freq_Natural, ytop = (0:20)+0.4, col = '#3b674b')
rect(xleft = 0, ybottom = 0:20-0.4, xright = -7*baseR_freq_data2$ifreq_Natural, ytop = (0:20)+0.4, col = '#60a479')

#### time of interaction ####
db_behav %>% filter(entrance | interaction | eat | mark, Species %in% c("coati", "peccary_collared", "agouti", "rodent", "ocelot", "opossum_common")) %>% # filter to keep interactions
  mutate(Species = factor(Species, levels = baseR_freq_data2$Species, labels = baseR_freq_data2$Common)) %>% 
  group_by(Species, cavity_type) %>% summarise(meandur = mean(event_dur), se = sd(event_dur)/sqrt(n())) %>% # get mean and st. err.
  ggplot(aes(y = Species))+ # forest plot
  geom_pointrange(aes(x = meandur, xmax = meandur+se, xmin=meandur-se, color = cavity_type), position = position_dodge(width = 0.5))+
  theme_classic(base_size = 16)+
  scale_color_manual(values = c('#d8cb39','#3b674b'))+
  labs(x = 'Interaction duration (s)', 
       color = 'Cavity type')
#plot with panels, for six species that interact with both
db_behav %>% filter(entrance | interaction | eat | mark, Species %in% c("coati", "peccary_collared", "agouti", "rodent", "ocelot", "opossum_common")) %>% # filter to keep interactions
  mutate(Species = factor(Species, levels = baseR_freq_data2$Species, labels = baseR_freq_data2$Common)) %>% 
  ggplot(aes(cavity_type, event_dur, fill = cavity_type))+ # forest plot
  geom_boxplot(show.legend = F)+
  facet_wrap(~Species, scales = 'free_y')+
  theme_classic(base_size = 16)+
  "theme(strip.background = element_blank())+"
  scale_fill_manual(values = c('#d8cb39','#3b674b'))+
  labs(y = 'Interaction duration (s)', x='Cavity type')
