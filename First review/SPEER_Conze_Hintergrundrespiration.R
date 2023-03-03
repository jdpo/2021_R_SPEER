{
{
library(nlme)
library(dplyr)
library(ggplot2)
library(broom)
library(tidyverse)
library(data.table)
library(lubridate)
library(reshape2)
options(scipen = 999)
}
  
######### create a single df with all measurements and order ################
{
  source_directory <- "C:/Users/pohlmann/Desktop/SPEER/Hintergrundrespiration/Conze/edited"          #HERE  DEFINE DIRECTORY WHERE YOUR SOURCE DATA IS (NO FILENAME!). DIRECTORY SHOULD ONLY CONTAIN FILES FROM ONE MEASUREMENT, BUT THEY NEED TO BE IN A SINGLE DIRECTORY WITHOUT SUBFOLDERS
  
  #bind all csvs  
  setwd(source_directory)    
  mycsv <- dir(pattern=".csv")
  files <- list.files(pattern="*.csv")
  alldata <- do.call(rbind,lapply(files, fread))                    
  
  #rename columns
  alldata <- alldata %>%
    rename(date = Datum,
           time = Zeit,
           conc = O2.Gehalt..mg.l.,
           sat = O2.Sättigung....,
           temp = Temperatur...C.,
           press = Druck..bar.,
           speed = Strömung..m.s.,
           flow = Durchfluss..l.h.,
           mode = Modus,
           press_s = Druck.Soll..bar.,
           speed_s = Strömung.Soll..m.s.,
           temp_s = Temperatur.Soll...C.,
           front = Lichtschranke.Vorne,
           back = Lichtschranke.Hinten,
           measuring = Messung
    )
  
  #order columns
  alldata <- as.data.frame(alldata)
  col_order <- c("measurement", "tunnel", "date", "time", "conc", "sat", "temp", "temp_s", "press", "press_s", "speed", "speed_s", "front", "back", "UTC", "measuring", "flow", "mode")
  alldata <- alldata[, col_order]
  
  #clean environment
  rm(list = ls()[!ls() %in% c("alldata")])
  alldata$Index <- seq.int(nrow(alldata))
}  

##############################################sort & clear #################################
{
{
backup1 <- alldata
alldata <- backup1

alldata$date <- as.Date(alldata$date, "%d.%m.%Y")
alldata$datetime <- as.POSIXct(paste(alldata$date, alldata$UTC), format="%Y-%m-%d %H:%M:%S", tz = "UTC")

alldata <- alldata[with(alldata, order(measurement, tunnel, datetime)),]

alldata <- alldata%>%
  mutate(diff = datetime - lag(datetime, default = first(datetime)))

{
  m <- numeric(nrow(alldata)) 
  m[1] <- alldata$measuring[1]
  
  for (i in 2:nrow(alldata)){ 
    if (alldata$measurement[i] != alldata$measurement[i-1] | alldata$tunnel[i] != alldata$tunnel[i-1]) { 
      m[i] <- 0                  
    } else { 
      if (alldata$flow[i] == 0 && alldata$flow[i-1] != 0) { 
        m[i] <- m[i-1]+1                   
      } else {
        if (alldata$flow[i] == 0 && alldata$flow[i-1] == 0 && alldata$diff[i]>1800){
          m[i] <- m[i-1]+1                                    
        } else {
          m[i] <- m[i-1]
        }
        
      }
      
    }
    
  }
  alldata$no <- m
  alldata$cycle_no <- paste(alldata$measurement, alldata$tunnel, alldata$no, sep = "_")
  
}



backup2 <- alldata
alldata <- backup2

#remove all values with erroneous measurements
alldata <- alldata[!(alldata$no > 3),]
alldata <- alldata[!(alldata$temp > 25 | alldata$temp < 2),]
alldata <- alldata[!(alldata$press < 0 | alldata$press > 10),]
alldata <- alldata[!(alldata$conc > 20),]
alldata <- alldata[!(alldata$measuring == 0),]
alldata <- alldata[!(alldata$flow > 0),]
}

#------------------------------plot raw data to check-----------------------------#

{
#create time scale per measurement and tunnel
alldata$mt <- paste(alldata$measurement, alldata$tunnel, sep = "_")

alldata <- alldata %>%
  group_by(mt) %>%
  mutate(starttime_mt = min(datetime))

alldata$mt_time <- as.numeric(alldata$datetime-alldata$starttime_mt)

options(scipen = 999) #turn off scientific notations
library(ggplot2)
library(scico)
library(egg)

#facet lable names
tunnel.names <- c("tunnel 1", "tunnel 2", "tunnel 3")
names(tunnel.names) <- c("1", "2", "3")

ggplot(alldata, aes(x=mt_time, y=conc))+ 
  geom_point(aes(col = temp, shape = as.factor(alldata$no)), size=0.9)+ 
  theme_bw()+
  scale_colour_scico(palette = "batlow", name = "temp (°C)")+
  labs(x = "time (in sec)", y = expression((O[2]-concentration ~(mg %*% l^-1))))+
  facet_grid(measurement ~ tunnel, labeller = labeller(tunnel = tunnel.names))+
  guides(shape=FALSE)

#remove redundant columns
alldata <- alldata[,! names(alldata) %in% c("measuring", "flow", "mt_time", "mt", "starttime_mt", "no"), drop = F]

rm(list = ls()[!ls() %in% c("alldata", "backup1", "backup2")])
}

#---------------------------------FORMAT REMAINING DATA-------------------------------#

{
  #create a column with time in cycle
  alldata <- alldata %>%
    group_by(cycle_no) %>%
    mutate(starttime_c = min(datetime))
  
  alldata$c_time <- ((alldata$datetime-alldata$starttime_c)/3600)
  alldata <- alldata[,! names(alldata) %in% "starttime_c", drop = F]
  
  
  
  #create a column with timedifference between two measurements
  {
    d <- numeric(nrow(alldata)) 
    d[1] <- NA
    
    for (i in 2:nrow(alldata)){ 
      if (alldata$cycle_no[i] == alldata$cycle_no[i-1]) { 
        d[i] = alldata$c_time[i] - alldata$c_time[i-1]                   
      } else {
        d[i] = NA                                    
      }
      
    } 
    alldata$timediff <- d
  }
  
  #backup alldata
  backup3 <- alldata
  
  #clean environment
  rm(list = ls()[ls() %in% c("d", "i")])
  
}
}

############## RUN REGRESSION, CREATE A SUMMARY FILE WHERE EACH MEASUREMENT CYCLE IS REDUCED TO ONE ROW, CONTAINING SLOPE AND ADDITIONAL DATA #########################

{
  library(broom)
  
  alldata <- group_by(alldata, cycle_no)
  
  #------------------CREATE TABLES WITH MAXIMA, MINIMA, SLOPE ESTIMATES etc. AND MERGE THEM---------------#
  
  #create tables with max, min, mean, SD and counts for each cycle  
  stats <- alldata %>% 
    group_by(cycle_no) %>% 
    summarise(datetime_min = min(datetime), 
              conc_min = min(conc),
              sat_min = min(sat),
              temp_s_min = min(temp_s),
              temp_min = min(temp),
              press_s_min = min(press_s),
              press_min = min(press),
              speed_s_min = min(speed_s),
              speed_min = min(speed),
              tunnel = min(tunnel),
              c_time_max = max(c_time),
              datetime_max = max(datetime), 
              conc_max = max(conc),
              sat_max = max(sat),
              temp_s_max = max(temp_s),
              temp_max = max(temp),
              press_s_max = max(press_s),
              press_max = max(press),
              speed_s_max = max(speed_s),
              speed_max = max(speed),
              conc_mean = mean(conc),
              timediff_max = max(timediff, na.rm = TRUE), #max(na.omit(timediff))also works, just be aware that omit ignores the full row whereas na.rm only ignores the NA for the given column. Makes no difference here but e.g. for colMeans na.omit will affect all columns!              
              sat_mean = mean(sat),
              temp_mean = mean(temp),
              press_mean = mean(press),
              speed_mean = mean(speed),
              conc_sd = sd(conc),
              sat_sd = sd(sat),
              temp_sd = sd(temp),
              press_sd = sd(press),
              speed_sd = sd(speed),
              n = sum(!is.na(conc)),
              temp_count = sum(!is.na(temp)),
              measurement = min(measurement)
    ) 
  
  #produce data frames with results of linear regression, tidy = parameter estimates, glance = regression statistics
  result_tidy <- do(alldata, tidy(lm(conc ~ c_time, data = .)))
  
  result_tidy_intercept <- result_tidy[which(result_tidy$term == "(Intercept)"),]
  colnames(result_tidy_intercept)[2:ncol(result_tidy_intercept)] <- (paste(colnames(result_tidy_intercept)[2:ncol(result_tidy_intercept)], "i", sep ="_"))
  
  result_tidy_slope <- result_tidy[which(result_tidy$term == "c_time"),]
  colnames(result_tidy_slope)[2:ncol(result_tidy_slope)] <- (paste(colnames(result_tidy_slope)[2:ncol(result_tidy_slope)], "s", sep ="_"))
  
  result_glance <- do(alldata, glance(lm(conc ~ c_time, data = .)))
  colnames(result_glance)[2:ncol(result_glance)] <- (paste(colnames(result_glance)[2:ncol(result_glance)], "lm", sep ="_"))
  
  #Merging all the different data.frames created above: maxs, mins, means, sds, counts tidy and glance by "no" (essentially VLOOKUP in excel, solution for different col_names see below)
  result <- merge(result_glance, result_tidy_intercept , by = "cycle_no", all = TRUE)
  result <- merge(result, result_tidy_slope, by = "cycle_no", all = TRUE)
  result <- merge(result, stats, by = "cycle_no", all = TRUE)
  
  #Check if values in temp_count match values in conc_count, if all match, remove column temp_count (conc_count is named "n" already), if not, give rows that do not match 
  if (nrow(as.matrix(which(result$temp_count-result$n != 0))) == 0){
    result <- result[,! names(result) %in% c("temp_count"), drop = F]
  }else{
    rownames(result[!result$temp_count %in% result$n, ])  
  } 
  
  
  #drop unnecessary columns and create the data.frame "summary"
  summary <- result[,! names(result) %in% c("r.squared_lm", "df_lm", "term_i", "term_s", "statistic_lm", "statistic_i", "statistic_s","logLik_lm","AIC_lm", "BIC_lm"), drop = F]
  
  #Sauerstoffverbrauch pro Liter berechnen 
  summary$resp <- (summary$estimate_s*-205)
  
  #order columns
  col_order <- c( "measurement", "datetime_min", "datetime_max", "c_time_max", "temp_min", "temp_mean", "temp_max", "temp_sd", "press_min", "press_mean", "press_max", "press_sd", "n", "resp", "std.error_s", "p.value_s", "estimate_i", "std.error_i", "p.value_i", "conc_min", "conc_mean", "conc_max", "conc_sd", "speed_min", "speed_mean", "speed_max", "speed_sd", "sat_min", "sat_mean", "sat_max", "sat_sd", "temp_s_min", "temp_s_max", "press_s_min", "press_s_max", "speed_s_min", "speed_s_max", "tunnel", "cycle_no", "adj.r.squared_lm", "sigma_lm", "p.value_lm", "deviance_lm", "df.residual_lm", "timediff_max")
  summary <- summary[, col_order]
  
  
  #clean environment
  rm(list = ls()[!ls() %in% c("summary")])
  
  #create summary backup
  summary_raw <- summary
  summary <- summary_raw
}

#####################RUN CHECKS AND OUTPUT FILES#######################

{
  
  output_removals <- "C:/Users/pohlmann/Desktop/SPEER/Hintergrundrespiration/Conze/removals.csv"
  output_summary <-   "C:/Users/pohlmann/Desktop/SPEER/Hintergrundrespiration/Conze/summary.csv"
  output_counts <-   "C:/Users/pohlmann/Desktop/SPEER/Hintergrundrespiration/Conze/counts.csv"
  
  NaNs <- summary %>% filter_all(any_vars(grepl("NaN", .)))
  if (nrow(NaNs) > 0){NaNs$condition <- "relevant columns contain NaNs"}
  summary <- summary[!summary$cycle_no %in% NaNs$cycle_no, ]
  
  #Run checks on given conditions and 
  r_squared <- summary[(summary$adj.r.squared_lm < 0.0),]
  if (nrow(r_squared) > 0){r_squared$condition <- "r-squared < 0.0"}
  
  resp <- summary[(summary$resp <= 0),]
  if (nrow(resp) > 0){resp$condition <- "Respiration <= 0"}
  
  n <- summary[(summary$n <= 30),]
  if (nrow(n) > 0){n$condition <- "no of datapoints < 30"} 
  
  distance <- summary[(summary$timediff_max >= 0.01666667),]
  if (nrow(distance) > 0){distance$condition <- "datapoints > 60sec apart"}
  
  duration <- summary[(summary$c_time_max < 0.5),]
  if (nrow(duration) > 0){duration$condition <- "measurement less than 30min"}
  
  temp <- summary[(summary$temp_max-summary$temp_min > 1.5),]
  if (nrow(temp) > 0){temp$condition <- "Temperature varies by more than 1.5°C"}
  
  press <- summary[(summary$press_mean < 0.8 | summary$press_mean > 1.2 & summary$press_mean < 7.8),]
  if (nrow(press) > 0){press$condition <- "Mean pressure more than 0.2 from 1 or 8 bar"}
  
  press_sd <- summary[(summary$press_sd > 0.1),]
  if (nrow(press_sd) > 0){press_sd$condition <- "Mean pressure sd > 0.1"}
  
  unreal <- summary[(summary$resp/mean(summary$resp) > 30),]
  if (nrow(unreal) > 0){unreal$condition <- "respiration 30x more then average of all resp"}
  
  sat <- summary[(summary$sat_max < 80),]
  if (nrow(sat) > 0){sat$condition <- "max o2-Saturation < 90%"}
  
  p_s <- summary[(summary$p.value_s > 0.05),]
  if (nrow(p_s) > 0){p_s$condition <- "slope not significant"}
  
  
  #create a list of the dataframes created above
  conditions <- list(r_squared, resp, n, distance, duration, temp, press, press_sd, unreal, sat, NaNs, p_s)
  
  #summarize all removals and write table
  all_removals_full <- do.call(rbind,conditions)
  removals_reason <- reshape2::dcast(all_removals_full, cycle_no ~ condition, value.var = "condition", fun.aggregate = length)
  
  all_removals <- all_removals_full[!duplicated(all_removals_full$cycle_no),]
  all_removals <- merge(all_removals, removals_reason, by = "cycle_no", all = TRUE)
  all_removals <- all_removals[,! names(all_removals) %in% c("condition"), drop = F]
  all_removals <- all_removals[order(all_removals$tunnel, all_removals$datetime_min),]
  all_removals <- all_removals[rowSums(is.na(all_removals)) != ncol(all_removals), ]
  write.table(all_removals, file = output_removals, row.names = FALSE, sep = ";")
  
  #create table with the sum of removals for respective reason
  count_sums <- all_removals %>%
    bind_rows(summarise(.,
                        across(where(is.numeric), sum),
                        across(where(is.character), ~"Total")))
  
  count_sums <- count_sums[nrow(count_sums),]
  count_sums <- t(count_sums)
  count_sums <- cbind(Condition = rownames(count_sums), count_sums)
  rownames(count_sums) <- 1:nrow(count_sums)
  count_sums <- as.data.frame(count_sums)
  names(count_sums)[2] <- "Total"
  count_sums <- count_sums[46:nrow(count_sums),]
  write.table(count_sums, file = output_counts, row.names = FALSE, sep = ";")
  
  #remove the data that did not meet conditions from summary
  summary <- summary[!summary$cycle_no %in% all_removals$cycle_no, ]
  summary <- summary[!(summary$tunnel == 3 & summary$measurement == 3),]
  write.table(summary, file = output_summary, row.names = FALSE, sep = ";")
  
  rm(list = ls()[!ls() %in% c("summary" , "summary_raw", "all_removals")])
}
}
#####################################VISUALIZE DATA##########################################

{
  #remove all but tunnel3 (only one with visible temperature dependency)
  
  
  options(scipen = 999) #turn off scientific notations
  library(ggplot2)
  setwd("C:/Users/pohlmann/Desktop/SPEER/Ergebnisse")
  devtools::install_github("thomasp85/scico")
  library(scico)
  library(egg)
  
  
  #facet lable names
  tunnel.names <- c("tunnel 1", "tunnel 2", "tunnel 3")
  names(tunnel.names) <- c("1", "2", "3")
  
  ggplot(summary, aes(x=datetime_min, y=resp))+ 
    geom_point(aes(col = temp_mean), size=0.9)+ 
    theme_bw()+
    scale_colour_scico(palette = "batlow", name = "temp (°C)")+
    labs(x = expression(Date), y = expression(O[2]-consumption ~(mg %*% h^-1)))+
    facet_grid(measurement ~ tunnel, labeller = labeller(tunnel = tunnel.names))+
    guides(shape=FALSE)
  
  ggplot(summary, aes(x=datetime_min, y=resp))+ 
    geom_point(aes(col = temp_mean), size=0.9)+ 
    theme_bw()+
    scale_colour_scico(palette = "batlow", name = "temp (°C)")+
    labs(x = expression(Date), y = expression(O[2]-consumption ~(mg %*% h^-1)))+
    facet_grid(. ~ tunnel, labeller = labeller(tunnel = tunnel.names))+
    guides(shape=FALSE)
  
}

#######################Test if HR changed during experiment##########################


lm <- lm(resp~as.factor(tunnel)-1, data=summary)
summary(lm)

summary <- summary %>%
  group_by(tunnel) %>%
  mutate(trainingstart = min(datetime_min))

summary$tunnel_time <- (summary$datetime_min-summary$trainingstart)/3600

lm.time <- lm(resp~tunnel_time*as.factor(tunnel), data=summary)
summary(lm.time)

lm.temp <- lm(temp_mean~as.factor(tunnel)-1, data=summary)
summary(lm.temp)

