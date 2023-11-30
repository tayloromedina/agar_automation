
### split by acid and base

acids <- plate_means %>% 
  filter(sample %in% c('control', 'pH 2', 'pH 3', 'pH 4', 'pH 5'))
acid_stats <- titer_stats %>% 
  filter(sample %in% c('control', 'pH 2', 'pH 3', 'pH 4', 'pH 5'))

acidscols <- c('darkgrey', 'firebrick', 'coral', 'goldenrod', 'forestgreen')

ggplot() +
  
  geom_point(data = acids, aes(x = time, y = titer_mean, color = sample), alpha = 0.6, size = 5,
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.3)) +
  
  geom_point(data = acid_stats, aes(x = time, y = avg_titer, color = sample), shape = '\u2015', size = 8,
             position = position_dodge(0.3), show.legend = FALSE) +
  geom_line(data = acid_stats, aes(x = time, y = avg_titer, color = sample), alpha = 0.3, show.legend = FALSE) +
  
  geom_errorbar(data = acid_stats, aes(x = time, ymin = avg_titer - sd_titer, ymax = avg_titer + sd_titer, color = sample),
                position = position_dodge(0.3), width = 0) +
  
  geom_hline(data = experiment, yintercept = detectionlimit, linetype = '3313', size = 1.2) +
  
  annotate(geom = 'text', label = 'detection limit', x = 0, y = detectionlimit, vjust = -1) +
  
  scale_y_log10() +
  scale_color_manual(values = acidscols) +
  labs(x = experiment$x,
       y = experiment$y) +
  theme(panel.border = element_rect(colour = "black", fill= NA, size = 0.75), # black plot border
        panel.background = element_blank(), # white panel background
        panel.grid.major = element_line(color = 'lightgrey', linetype = 'dashed'),
        text = element_text(size = 15)) # light grey dashed plot lines




bases <- plate_means %>% 
  filter(sample %in% c('control', 'pH 8', 'pH 9', 'pH 10', 'pH 11'))
base_stats <- titer_stats %>% 
  filter(sample %in% c('control', 'pH 8', 'pH 9', 'pH 10', 'pH 11'))

basecols <- c('darkgrey', 'black','darkorchid3', 'turquoise', 'dodgerblue')

ggplot() +
  
  geom_point(data = bases, aes(x = time, y = titer_mean, color = sample), alpha = 0.6, size = 5,
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.3)) +
  
  geom_point(data = base_stats, aes(x = time, y = avg_titer, color = sample), shape = '\u2015', size = 8,
             position = position_dodge(0.3), show.legend = FALSE) +
  geom_line(data = base_stats, aes(x = time, y = avg_titer, color = sample), alpha = 0.3, show.legend = FALSE) +
  
  geom_errorbar(data = base_stats, aes(x = time, ymin = avg_titer - sd_titer, ymax = avg_titer + sd_titer, color = sample),
                position = position_dodge(0.3), width = 0) +
  geom_hline(data = experiment, yintercept = detectionlimit, linetype = '3313', size = 1.2) +
  
  annotate(geom = 'text', label = 'detection limit', x = 0, y = detectionlimit, vjust = -1) +
  
  scale_y_log10() +
  scale_color_manual(values = basecols) +
  labs(x = experiment$x,
       y = experiment$y) +
  theme(panel.border = element_rect(colour = "black", fill= NA, size = 0.75), # black plot border
        panel.background = element_blank(), # white panel background
        panel.grid.major = element_line(color = 'lightgrey', linetype = 'dashed'),
        text = element_text(size = 15)) # light grey dashed plot lines +



### rate constant plots
ratecols <- c('firebrick', 'coral', 'goldenrod', 'forestgreen', 'darkgrey', 'turquoise', 'dodgerblue', 'black','darkorchid3')
ggplot() +
  
  geom_bar(data = constants_edit, aes(x = sample, y = k, fill = sample), stat = 'identity', alpha = 0.6) +
  
  geom_errorbar(data = combined, aes(x = sample, ymin = x2_5_percent, ymax = x97_5_percent), width = 0.2) +
  
  geom_segment(data = kmax.calc, aes(x = sample, y = k, xend = sample, yend = k + 0.08 * k), arrow = arrow(),
               lineend = 'round', linejoin = 'round', size = 1) +
  
  ggbreak::scale_y_break(c(12, 30)) +
  labs(x = 'pH',
       y = 'k (hr-1)') +
  scale_fill_manual(values = ratecols) +
  theme(panel.border = element_rect(colour = "black", fill= NA, size = 0.75), # black plot border
        panel.background = element_blank(), # white panel background
        panel.grid.major = element_line(color = 'lightgrey', linetype = 'dashed'),
        text = element_text(size = 15))  # light grey dashed plot lines

ggplotstring2 <- paste(experiment$EXPERIMENT, experiment$DATE, 'rateconstant', sep = '_')
ggsave(here('output', paste(ggplotstring2, '.png')), width = 30, height = 25, units = 'cm')

