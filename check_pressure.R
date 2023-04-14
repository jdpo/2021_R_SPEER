test <- summary %>% 
  filter(press_min > 1.5 & press_min < 7.5)


data_wrong <-  data[data$no %in% test$no ,]

graph_all_press <- ggplot(data_wrong, aes(x=m_time, y = press)) +
  geom_point()+
  facet_wrap(~no)

ggsave("check_pressure.pdf", graph_all_press, device = "pdf", height = 7, width = 8, dpi = 300)
