

# SCRIPT for automating double-layer plaque assay enumerations and plate-spotting enumerations

# libraries
library(tidyverse)
library(readxl)
library(here)
library(writexl)


# load data
enumeration <- read_excel(here('input', 'template.xlsx'), 1)
experiment <- read_excel(here('input', 'template.xlsx'), 2) 


# clean up the data
enumeration_clean <- enumeration %>% 
  janitor::clean_names() %>% 
  filter(!is.na(count))

platevolume <- experiment$VOLUME %>% 
  as.numeric()


# titer calculations
plate_titer <- enumeration_clean %>% 
  mutate(titer = count * 10 ^ dilution / platevolume)

plate_means <- plate_titer %>% 
  group_by(sample, time, replicate) %>% 
  summarize(titer_mean = psych::geometric.mean(titer))


# titer statistics
titer_stats <- plate_means %>% 
  group_by(sample, time) %>% 
  summarize(avg_titer = mean(titer_mean, na.rm = TRUE),
            sd_titer = sd(titer_mean, na.rm = TRUE)) 


# generates the excel sheet that has the titer information by sample
string <- paste(experiment$DATE, experiment$EXPERIMENT, 'titer', sep = '_')

detectionlimit <- experiment$`detection limit` # extracts the detection limit 

plate_means_dl <- plate_means %>% 
  filter(titer_mean != 0)

# replaces all '0' concentrations with the BDL 
plate_means <- plate_means %>% 
  mutate(titer_mean = ifelse(titer_mean == 0, detectionlimit / sqrt(2), titer_mean))

titer_stats <- titer_stats %>% 
  mutate(avg_titer = ifelse(avg_titer == 0, detectionlimit / sqrt(2), avg_titer))


# rate constant calculations 
c0 <- plate_means_dl %>% 
  filter(time == 0)

logs <- inner_join(plate_means_dl, c0, by = c('sample', 'replicate')) %>% 
  select(!time.y ) %>% 
  mutate(log_removal = - log(titer_mean.x / titer_mean.y)) 


# creates linear models to calculate a 1st order rate constant
models <- plyr::dlply(logs, 'sample', function (d)
  lm(log_removal ~ time.x, data = d))


# extracts all of the rate constants from the list
constants <- data.frame()

for (i in 1:length(models)) {
  
  fits <- broom::tidy(models[[i]]) 
  value <- fits$estimate[2]
  samples <- names(models)[i]
  combined <- c(samples, value)
  constants <- rbind(constants, combined)
  
}

colnames(constants) <- c('sample', 'k')
print(constants)


# get the 95% confidence interval of all of the models
fitted_models <- logs %>%
  group_by(sample) %>%
  do(model = lm(log_removal ~ time.x, data = .)) %>%
  ungroup()

t <- purrr::map_dfr(
  fitted_models$model,
  ~ as_tibble(confint(., "time.x", level = 0.95)))

combined <- bind_cols(fitted_models, t) %>% 
  mutate(sample = ifelse(sample == 'control', 7.4, sample),
         sample = ifelse(str_detect('pH', sample), str_remove(sample, 'pH '), sample),
         sample = as.numeric(sample)) %>% 
  janitor::clean_names()

# combined <- combined[order(combined$sample), ] %>% 
#   mutate(sample = as.factor(sample))


### generates the scatter plot to show data vs. sample time
ggplot() +
  
  geom_point(data = plate_means, aes(x = time, y = titer_mean, color = sample), alpha = 0.6, size = 5,
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.3)) +
  
  geom_point(data = titer_stats, aes(x = time, y = avg_titer, color = sample), shape = '\u2015', size = 8,
             position = position_dodge(0.3), show.legend = FALSE) +
  
  geom_line(data = titer_stats, aes(x = time, y = avg_titer, color = sample), alpha = 0.3, show.legend = FALSE) +
  
  geom_errorbar(data = titer_stats, aes(x = time, ymin = avg_titer - sd_titer, ymax = avg_titer + sd_titer, color = sample),
                position = position_dodge(0.3), width = 0) +
  
  geom_hline(data = experiment, yintercept = detectionlimit, linetype = '3313', size = 1.2) +
  
  annotate(geom = 'text', label = 'detection limit', x = 0, y = detectionlimit, vjust = -1) +
  
  scale_y_log10() +
  labs(x = experiment$x,
       y = experiment$y) +
  theme(panel.border = element_rect(colour = "black", fill= NA, size = 0.75), # black plot border
        panel.background = element_blank(), # white panel background
        panel.grid.major = element_line(color = 'lightgrey', linetype = 'dashed'),
        text = element_text(size = 15)) # light grey dashed plot lines

ggplotstring <- paste(experiment$EXPERIMENT, experiment$DATE, 'plot', sep = '_')

# saves the scatter plot - you can change the width and height and units if you want different plot sizes
ggsave(here('output', paste(ggplotstring, '.png')), width = 35, height = 25, units = 'cm')



# edits the rate constant data for the plot
constants_edit <- constants %>% 
  mutate(k = as.numeric(k),
         sample = ifelse(sample == 'control', 7.4, sample),
         sample = str_remove(sample, 'pH '),
         sample = as.character(sample)) 

# calculates the kmax of BDL samples
kmax <- constants_edit %>% 
  filter(is.na(k)) %>% 
  mutate(sample = sub('^', 'pH ', sample)) %>% 
  inner_join(titer_stats, by = 'sample') %>% 
  group_by(sample) %>% 
  filter(row_number() %in% c(1,2))

kmax.calc <- kmax %>% 
  select(sample, time, avg_titer) %>% 
  group_by(sample) %>% 
  mutate(k = -log(lead(avg_titer) / avg_titer) / lead(time)) %>% 
  filter(!is.na(k)) %>% 
  select(sample, k) %>% 
  mutate(sample = str_remove(sample, 'pH '))

# combines the fitted constants with the kmax constants  
constants_edit <- left_join(constants_edit, kmax.calc, by = 'sample') %>% 
  mutate(k = coalesce(k.x, k.y)) %>% 
  select(sample, k) 

# the original code was written with pH in mind, so this section fixes the outputs to be more robust
if (str_detect('pH', experiment$EXPERIMENT)) {
  
  next 
  
} else {
  
  constants_edit <- constants_edit %>%
    mutate(sample = ifelse(sample == 7.4, 'control', sample))
  
  combined <- combined %>% 
    mutate(sample = ifelse(sample == '7.4', 'control', sample))
  
}


# writes the excel sheet that has the titer information and the rate constant
write_xlsx(list('titer' = titer_stats, 'k' = constants_edit), here('output', paste0(string, '.xlsx')))

# plots the rate constant data as a function of pH 
ggplot() +
  
  geom_bar(data = constants_edit, aes(x = sample, y = k, fill = sample), stat = 'identity', alpha = 0.6) +
  
  geom_errorbar(data = combined, aes(x = sample, ymin = x2_5_percent, ymax = x97_5_percent), width = 0.2) +
  
  geom_segment(data = kmax.calc, aes(x = sample, y = k, xend = sample, yend = k + 0.08 * k), arrow = arrow(),
               lineend = 'round', linejoin = 'round', size = 1) +
  
  labs(x = 'pH',
       y = 'k (hr-1)') +
  theme(panel.border = element_rect(colour = "black", fill= NA, size = 0.75), # black plot border
        panel.background = element_blank(), # white panel background
        panel.grid.major = element_line(color = 'lightgrey', linetype = 'dashed'),
        text = element_text(size = 15))  # light grey dashed plot lines

ggplotstring2 <- paste(experiment$EXPERIMENT, experiment$DATE, 'rateconstant', sep = '_')
ggsave(here('output', paste(ggplotstring2, '.png')), width = 30, height = 25, units = 'cm')


# plots the c/c0 to see if there is any significant difference for each replicate
cc0 <- logs %>% 
  mutate(logremoval = log_removal) %>% 
  select(sample, time.x, logremoval, replicate) %>% 
  mutate(day = case_when(
    replicate < 3 ~ 1,
    replicate < 5 ~ 2,
    replicate < 7 ~ 3,
    replicate < 9 ~ 4,
    replicate < 11 ~ 5,
    replicate < 13 ~ 6
  ))

ggplot(data = cc0, aes(x = time.x, y = logremoval, color = as.factor(day))) +
  geom_point(alpha = 0.6, size = 5) +
  geom_smooth(se = F) +
  facet_wrap(~ sample) +
  labs(x = experiment$x,
       y = 'ln(C/C0)',
       color = 'day') +
  theme(panel.border = element_rect(colour = "black", fill= NA, size = 0.75), # black plot border
        panel.background = element_blank(), # white panel background
        panel.grid.major = element_line(color = 'lightgrey', linetype = 'dashed'),
        text = element_text(size = 15)) # light grey dashed plot lines
  
