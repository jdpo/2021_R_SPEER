#swimming time in days

load("./2021_SPEER/stats.RData")

library(tidyverse)

swimtime <- stats %>% 
  group_by(ID) %>% 
  summarize(start = min(datetime_min),
            end = max(datetime_max),
            avg = mean(speed_mean)) %>% 
  mutate(time_d = as.numeric(end-start),
         time_h = as.numeric((end-start)*24),
         distance = (avg*(time_h*60*60))/1000)

# treatment times 
stats_treat <- stats %>% 
  mutate(treat = paste(ID, ntemp, npress, sep = "_"))

stats_treat <- stats_treat[order(stats_treat$run, stats_treat$tunnel, stats_treat$datetime_min),]

p <- c(nrow(stats_treat))  
p[1] <- 1
for (i in 2:nrow(stats_treat)){
  ifelse(stats_treat$treat[i] != stats_treat$treat[i-1], p[i] <- p[i-1]+1  , p[i] <- p[i-1])}
stats_treat$diff <- p 

stats_treat <- stats_treat %>% 
  group_by(diff) %>% 
  mutate(diffbeg = as.numeric(datetime_min - min(datetime_min)),
         diffmax_1 = as.numeric(max(diffbeg)),
         dk = 0) %>% 
  ungroup() %>% 
  group_by(treat) %>% 
  mutate(diffmax_2 = max(diffmax_1),
         dk = ifelse(diffmax_1 != diffmax_2, "d", "k")) %>% 
  filter(dk == "k")


treattime <- stats_treat %>%
  mutate(treat = paste(ID, ntemp, npress, sep = "_")) %>% 
  group_by(treat) %>% 
  summarize(start = min(datetime_min),
            end = max(datetime_max)) %>% 
  mutate(time_d = as.numeric(end-start),
         time_h = as.numeric((end-start)*24))

treattime <- treattime %>%  
  left_join((stats_treat %>% select(ID, ntemp, npress, treat)), by = "treat") %>%
  distinct() %>% 
  select(ID, ntemp, npress, start, end, time_d, time_h, -treat)

