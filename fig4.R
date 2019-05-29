### Setup ####
#load packages and set working directory
data_wd <- 'C:\\Users\\John\\Google Drive\\Mouse_Microbiome_Shared_Files\\Data\\'
wd <- 'C:\\Users\\John\\Documents\\msu\\gut alternative states model'
setwd(data_wd)
library(tidyverse)
library(scales)
library(cowplot)

#load data 
#note we are using ZOTUs, not OTUs
otus <- readRDS(file = paste0(wd, 'zotus.RDS')) #load OTU table for plotting abundance-weighted coverage
tax_mouse <- readRDS(file = paste0(wd, 'tax_mouse_SILVA_zotus.RDS')) #tax mouse data only
meta <- readRDS(file = paste0(wd, 'meta.RDS'))

#change hour of sampling so that they are identical across runs (sometimes they are close but not identical)
meta$hour <- ifelse(meta$hour == 24, 25, meta$hour)
meta$hour <- ifelse(meta$hour == 48, 49, meta$hour)
meta$hour <- ifelse(meta$hour == 72, 73, meta$hour)
meta$hour <- ifelse(meta$hour == 96, 97, meta$hour)

#remove poorly sampled mice
meta <- meta %>% group_by(mouse) %>% filter(min(hour) < 4 & max(hour) >= 240)
otus <- filter(otus, sample %in% meta$sample)

# initial abundances of E and anaerobic mutualists
min_abun <- 260267551769

#a number that is well below the threshold of detection
min_abun <- 1000000000

# join metadata, remove faulty qpcr samples, calculate densities, summarize by mouse and hour, set "zero abundances" to minimum thresholds, calculate median abundances for each group, plot
p <- otus %>%
  left_join(meta, by = "sample") %>%
  filter(!is.na(cells_per_gram)) %>%
  mutate(abun = abun * cells_per_gram) %>%
  left_join(tax_mouse, by = "otu") %>%
  group_by(mouse, Abx, hour) %>%
  summarise(Clostriales = sum(abun[Phylum %in% c('Firmicutes','Bacteroidetes')], na.rm = T),
            Enterobacteriaceae = sum(abun[Family == 'Enterobacteriaceae'], na.rm = T)) %>%
  mutate(Clostriales = ifelse(Clostriales == 0, min_abun, Clostriales),
         Enterobacteriaceae = ifelse(Enterobacteriaceae == 0, min_abun, Enterobacteriaceae)) %>% 
  group_by(hour, Abx) %>% 
  summarise(C = median(Clostriales), E = median(Enterobacteriaceae)) %>%
  mutate(Abx = ifelse(Abx, 'Antibiotics-treated mice', 'Untreated control mice')) %>% 
  ggplot(aes(x = C, y = E, color = Abx)) + 
    geom_point() + 
    geom_path() + 
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_color_manual(values = c('#FF5555','#00AA00'), name = '') + 
    labs(x = 'Firmicutes and Bacteriodetes', y = 'Enterobacteriaceae') + 
    theme_classic() + 
    theme(
      legend.position = 'bottom',
      strip.background = element_blank(),
      plot.margin = unit(c(0, 0, 0, 0), "pt"),
      plot.title = element_text(hjust = 0.5))

setwd(wd)
#Print
p
ggsave('Fig4.pdf', width = 5.15, height = 4.75)

#print with lables
p + geom_text(aes(label = hour), color = 'black')
ggsave('Fig4_labeled.pdf', width = 5.15, height = 4.75)

