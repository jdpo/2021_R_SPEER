{
{
setwd("../Sonstiges/Backup step 2")

  
#install/load packages
library(nlme)
library(dplyr)
library(ggplot2)
library(broom)
library(tidyverse)
library(data.table)
library(lubridate)
library(reshape2)
}



############################################## READ DATA AND BRING IN CORRECT FORMAT #####################################################################

{
#read all csv files and create a single data frame named "data", group by no and format datetime as POSIXct
files <- list.files(pattern="*.csv")
data <- do.call(rbind,lapply(files, fread))
data$datetime <- as.POSIXct(data$datetime, format = "%Y-%m-%d %H:%M:%S")
data <- group_by(data, no)
}
  
############################################# EXTRACT SWIM SPEEDS FOR EACH ID##########################################################################
  
{
  SetSpeeds <- data %>% 
    group_by(ID,speed_s) %>%
    summarise(setSpeeds = n(),
              Front = sum(front %in% c("TRUE")),
              Back = sum(back %in% c("TRUE"))
              )


write.table(SetSpeeds, file = "../../SetSpeeds.csv", row.names = FALSE, sep = ";") 
rm(list = ls()[ls() %in% c("SetSpeeds")]) 

}

############################################# RECALCULATE ABSOLUTE 02-CCONSUMPTION BASED ON SATURATION ##################################################

{
  data$nopress <- data$sat/(1+(0.032*data$press)/(1000))
  data$tempK <- log((298.15-data$temp)/(273.15+data$temp))
  data$solut <- exp(2.00856 + 3.224*data$tempK+3.99063*data$tempK^2+4.80299*data$tempK^3+0.978188*data$tempK^4+1.71069*data$tempK^5+34.5*((-0.00624097)+(-0.00693498)*data$tempK+(-0.00690358)*data$tempK^2+(-0.00429155)*data$tempK^3)+0.00000031168*34.5^2)
  data$concm <- (data$solut*44.659*data$nopress)/100 
  data$conc1 <- data$concm/31.25
  data$concpress <- data$conc1*(1+(0.032*data$press)/1000)
  
  data$conc <- data$concpress
  
  data$nopress <- NULL
  data$tempK <- NULL
  data$solut <- NULL
  data$concm <- NULL 
  data$conc1 <- NULL
  data$concpress <- NULL
}

####################################### PRODUCE DATA FRAMES WITH DIFFERENT NEEDED INFO AND MERGE TO ONE ###################################################

{
#produce data frames with results of linear regression, tidy = parameter estimates, glance = regression statistics
result_tidy <- do(data, tidy(lm(conc ~ m_time, data = .)))
result_tidy_intercept <- result_tidy[which(result_tidy$term == "(Intercept)"),]
colnames(result_tidy_intercept)[2:ncol(result_tidy_intercept)] <- (paste(colnames(result_tidy_intercept)[2:ncol(result_tidy_intercept)], "i", sep ="_"))
result_tidy_slope <- result_tidy[which(result_tidy$term == "m_time"),]
colnames(result_tidy_slope)[2:ncol(result_tidy_slope)] <- (paste(colnames(result_tidy_slope)[2:ncol(result_tidy_slope)], "s", sep ="_"))

result_glance <- do(data, glance(lm(conc ~ m_time, data = .)))
colnames(result_glance)[2:ncol(result_glance)] <- (paste(colnames(result_glance)[2:ncol(result_glance)], "lm", sep ="_"))

#create tables with max, min, mean, SD and counts for each cycle
stats <- data %>% 
  group_by(no) %>% 
  summarise(datetime_min = min(datetime), 
            conc_min = min(conc),
            sat_min = min(sat),
            ntemp = min(ntemp),
            temp_s_min = min(temp_s),
            temp_min = min(temp),
            npress = min(npress),
            press_s_min = min(press_s),
            press_min = min(press),
            speed_s_min = min(speed_s),
            speed_min = min(speed),
            run = min(run),
            tunnel = min(tunnel),
            ID = min(ID),
            m_time_max = max(m_time),
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
            timediff_max = max(timediff, na.rm = TRUE), #max(na.omit(timediff))also works, just be aware that omit ignores the full row whereas na.rm only ignores the NA for the given column. Makes npo difference here but e.g. for colMeans na.omit will affect all columns!              
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
            back_true = sum(back %in% c("TRUE")),
            back_false = sum(back %in% c("FALSE")),
            front_true = sum(front %in% c("TRUE")),
            front_false = sum(front %in% c("FALSE"))
            )   
            
            

#Merging all the different data.frames created above: maxs, mins, means, sds, counts tidy and glance by "no" (essentially VLOOKUP in excel, solution for different col_names see below)
result <- merge(result_glance, result_tidy_intercept , by = "no", all = TRUE)
result <- merge(result, result_tidy_slope, by = "no", all = TRUE)
result <- merge(result, stats, by = "no", all = TRUE)

#Check if values in temp_count match values in conc_count, if all match, remove column temp_count (conc_count is named "n" already), if not, give rows that do not match 
if (nrow(as.matrix(which(result$temp_count-result$conc != 0))) == 0){
    result <- result[,! names(result) %in% c("temp_count"), drop = F]
    }else{
    rownames(result[!result$temp_count %in% result$n, ])  
    } 

}




###############################################################ADD BIOMETRY COLUMNS and ORDER/TIDY NEW DATA FRAME #####################################################################

{
#drop unnecessary columns and create the data.frame "summary"
summary <- result[,! names(result) %in% c("r.squared_lm", "df_lm", "term_i", "term_s", "statistic_lm", "statistic_i", "statistic_s","logLik_lm","AIC_lm", "BIC_lm"), drop = F]


#ADD A COLUMN WITH TIME AFTER START OF "TREATMENT"
summary <- summary %>%
    group_by(ID) %>%
    mutate(trainingstart = min(datetime_min))

summary$ID_time <- (summary$datetime_min-summary$trainingstart)/3600


#add columns with biomasses and length
biomass <- read.csv("../..//Metadaten/biomasses.csv", header = T, sep =";")
summary <- merge(summary, biomass, by.x = c("ID"), by.y = c("ID"))
summary <- summary %>%
  group_by(ID) %>%
  mutate(max_ID_time = max(ID_time))

#calculate respiration per kg/h from estimate_s and remove background respiration  (at this point body mass after experiment, needs better solution)
summary$max_ID_time_num <- as.numeric(summary$max_ID_time)
summary$ID_time_num <- as.numeric(summary$ID_time)
summary$mass <- summary$startmass - ((summary$startmass - summary$endmass)*(summary$ID_time_num/summary$max_ID_time_num))
summary$resp_raw  <- (summary$estimate_s*-206)
hr <- read.csv("../../Hintergrundrespiration/Conze/HR.csv", header = T, sep =";")
summary <- merge(summary, hr, by.x = c("tunnel"), by.y = c("tunnel"))
summary$hr_calc <- summary$hr*1.83027^((summary$temp_mean-summary$hr_temp)/10)
summary$resp_absol <- summary$resp_raw-summary$hr_calc
summary$resp <- (summary$resp_raw - summary$hr_calc)/(summary$mass/1000)


#calculate speed in BL
summary$speedBL <- summary$speed_mean/(summary$length/100)

#order the "summary" data frame
summary <- summary[,! names(summary) %in% c("trainingstart", "max_ID_time_num", "ID_time_num", "max_ID_time","hr_calc", "resp_kg", "hr", "hr_temp"), drop = F]
col_order <- c("ID", "datetime_min", "datetime_max", "m_time_max", "resp", "resp_raw", "resp_absol", "ntemp", "temp_min", "temp_mean", "temp_max", "npress", "press_min", "press_mean", "press_max",  "n", "estimate_s", "std.error_s", "p.value_s", "conc_max", "conc_min", "conc_mean", "conc_sd", "temp_sd", "press_sd", "speed_min", "speed_max", "speed_mean", "speed_sd", "sat_min", "sat_max", "sat_mean", "sat_sd", "temp_s_min", "temp_s_max", "press_s_min", "press_s_max", "speed_s_min", "speed_s_max", "run", "tunnel", "no", "adj.r.squared_lm", "sigma_lm", "p.value_lm", "deviance_lm", "df.residual_lm", "estimate_i", "std.error_i", "p.value_i", "timediff_max", "ID_time","mass", "startmass", "endmass", "length", "stage", "speedBL", "back_true", "back_false", "front_true", "front_false")
summary <- summary[, col_order]

summary_raw1 <- summary
#summary <- summary_raw1
}


#################################################ADD COLUMN WITH TIME AFTER PRESSURE CHANGE#############################################################

{
  summary <- summary[order(summary$run, summary$tunnel, summary$datetime_min),]
  p <- POSIXct(nrow(summary))
  p[1] <- summary$datetime_min[1]
  for (i in 2:nrow(summary)){ 
    if (summary$npress[i] != summary$npress[i-1]){
      p[i] <- summary$datetime_min[i]
    } else {
      p[i] <- p[i-1]
    }
  }
  summary$startpress <- p    

summary$presstime <- (summary$datetime_min-summary$startpress)/3600
summary <- summary[,! names(summary) %in% "startpress", drop = F]
}

#####################################################ADD COLUMN WITH TIME AFTER TEMPERATURE CHANGE######################################################
{
  for (i in 1:nrow(summary)){
    if (summary$temp_s_min[i] > 17.5){
      summary$ntemp[i] <- 19
    } else if (summary$temp_s_min[i] < 11.5 && summary$temp_s_min[i] > 2){
      summary$ntemp[i] <- 12
    } else if (summary$temp_s_min[i] < 0.45){
      summary$ntemp[i] <- 19
    } else if ((summary$temp_s_min[i] < 0.6 && summary$temp_s_min[i] > 0.45)){
      summary$ntemp[i] <- 15.5
    } else {
      summary$ntemp[i] <- 15.5 
    }
    }
 

  summary <- summary %>%
    group_by(ID, ntemp) %>%
    mutate(trainingstart = min(datetime_min))
  
  summary$temptime <- (summary$datetime_min-summary$trainingstart)/3600
  summary <- summary[,! names(summary) %in% "trainingstart", drop = F]
  
  #add column with "fishnumber_temp"
  summary$ID_temp <- paste(summary$ID, summary$ntemp,sep = "_")
  summary$ID_temp_press <- paste(summary$ID_temp, summary$npress, sep = "_")  
}  
  



#####################################################CREATE A ROLLBACK OF SUMMARY#######################################################################

{
#backup summary
summary_raw2 <- summary
#summary <- summary_raw2
}

######################################SET NEW WD AND REMOVE ALL REDUNDANT ENTRIES IN THE ENVIRONMENT####################################################

{
rm(list = ls()[!ls() %in% c("summary" , "summary_raw1", "summary_raw2")])
setwd("../../Ergebnisse/Bereinigung/")
}

######################################TEST CONDITIONS, OUTPUT COUNTS AND CASES OF REMOVALS####################################################################

{
pacclimatization <- summary[(summary$presstime < 24),]
if (nrow(summary) > 0){pacclimatization$condition <- "pressure acclimatization"}

tacclimatization <- summary[(summary$temptime < 30),]
if (nrow(summary) > 0){tacclimatization$condition <- "temperature acclimatization"}

r_squared <- summary[(summary$adj.r.squared_lm < 0.9),]
if (nrow(r_squared) > 0){r_squared$condition <- "r-squared < 0.9"}
 
resp <- summary[(summary$resp <= 0),]
if (nrow(resp) > 0){resp$condition <- "Respiration <= 0"} 

n <- summary[(summary$n <= 60),]
if (nrow(n) > 0){n$condition <- "no of measurements < 60"} 

distance <- summary[(summary$timediff_max >= 0.01666666667),]
if (nrow(distance) > 0){distance$condition <- "datapoints > 60sec apart"}

duration <- summary[(summary$m_time_max < 0.5),]
if (nrow(duration) > 0){duration$condition <- "measurement less than 30min"}

temp <- summary[(summary$temp_max-summary$temp_min > 1),]
if (nrow(temp) > 0){temp$condition <- "Temperature varies by more than 1?C"}

press <- summary[(summary$press_mean-summary$npress > 0.1 | summary$press_mean-summary$npress < -0.1),]
if (nrow(press) > 0){press$condition <- "Mean pressure < 0.1 from nominal pressure"}

press_s <- summary[(summary$press_s_max-summary$press_s_min > 1),]
if (nrow(press_s) > 0){press_s$condition <- "min and max set pressure > 1bar apart"}

unreal <- summary[(summary$resp/mean(summary$resp) > 50),]
if (nrow(unreal) > 0){unreal$condition <- "respiration 50x more then average of all resp"}

sat <- summary[(summary$sat_max < 90),]
if (nrow(sat) > 0){sat$condition <- "max o2-Saturation < 90%"}

p_lm <- summary[(summary$p.value_lm > 0.05),]
if (nrow(p_lm) > 0){p_lm$condition <- "Model not significant"}

p_i <- summary[(summary$p.value_i > 0.05),]
if (nrow(p_i) > 0){p_i$condition <- "Intercept not significant"}

p_s <- summary[(summary$p.value_s > 0.05),]
if (nrow(p_s) > 0){p_s$condition <- "Slope not significant"}

#create a list of the dataframes created above
conditions <- list(tacclimatization, pacclimatization, r_squared, unreal, n, distance, duration, temp, press, press_s, resp, sat, p_lm, p_i, p_s)

#summarize all removals and write table
all_removals_full <- do.call(rbind,conditions)
removals_reason <- dcast(all_removals_full, no ~ condition, value.var = "condition", fun.aggregate = length)


removals_reason <- removals_reason %>%
  rename(removed_tacclimatization = "temperature acclimatization",
         removed_pacclimatization = "pressure acclimatization",
         r_squared = "r-squared < 0.9",
         removed_respiration = "Respiration <= 0",
         removed_n = "no of measurements < 60",
         removed_distance = "datapoints > 60sec apart",
         removed_duration = "measurement less than 30min",
         removed_temp = "Temperature varies by more than 1?C",
         removed_press = "Mean pressure < 0.1 from nominal pressure",
         removed_press_s = "min and max set pressure > 1bar apart",
         removed_resp = "respiration 50x more then average of all resp",
         removed_sat = "max o2-Saturation < 90%",
         removed_lm = "Model not significant",
         removed_i = "Intercept not significant",
         removed_s = "Slope not significant"
         
  )

all_removals <- all_removals_full[!duplicated(all_removals_full$no),]
all_removals <- merge(all_removals, removals_reason, by = "no", all = TRUE)
all_removals <- all_removals[,! names(all_removals) %in% c("condition"), drop = F]
all_removals <- all_removals[order(all_removals$run, all_removals$tunnel, all_removals$datetime_min),]
write.table(all_removals, file = "all_removals_R09.csv", row.names = FALSE, sep = ";") 

#remove redundant columns from summary
summary <- summary[,! names(summary) %in% c("p.value_i", "p.value_lm", "p.value_s"), drop = F]

rm(list = ls()[!ls() %in% c("summary" , "summary_raw1", "summary_raw2", "all_removals", "all_removals_full")])
}



######################################REMOVE RESPECTIVE DATA FROM "SUMMARY" AND CREATE SIMPLE TABLE WITH REMOVED NUMBERS########################################################

{
summary_counts_raw <- summary %>%
  group_by(ID_temp_press) %>%
  summarise(count=n())

summary <- summary[!summary$no %in% all_removals$no, ]

summary_counts_edited <- summary %>%
  group_by(ID_temp_press) %>%
  summarise(count=n())

summary_counts <- merge(summary_counts_raw, summary_counts_edited, by = "ID_temp_press")

summary_counts <- summary_counts %>%
  rename(count_raw = count.x,
         count_edited = count.y
  )

summary_counts$removals_total <- summary_counts$count_raw - summary_counts$count_edited

summary_counts_reasons <- dcast(all_removals_full, ID_temp_press ~ condition, value.var = "condition", length)
summary_counts <- merge(summary_counts, summary_counts_reasons, by = "ID_temp_press", all.x = TRUE)
summary_counts[is.na(summary_counts)] = 0

summary_counts <- rbind(summary_counts, summarise_all(summary_counts, ~if(is.integer(.) | is.numeric(.)) sum(.) else "Total"))

write.table(summary_counts, file = "summary_counts_R09.csv", row.names = FALSE, sep = ";")
rm(list = ls()[!ls() %in% c("summary" , "summary_raw1", "summary_raw2")])

}


############################################ORDER DATA FRAME, CALCULATE THE TOTAL OXYGEN AND  AND EXPORT#################################################

{
summary <- summary[order(summary$run, summary$tunnel, summary$datetime_min),]
setwd("../../Ergebnisse/")
write.table(summary_raw2, file = "summary_raw.csv", row.names = FALSE, sep = ";")
write.table(summary, file = "summary_R09.csv", row.names = FALSE, sep = ";")
}
}


######################################################
#####################UNUSED CODE######################
######################################################

#setwd("C:/Users/pohlmann/Desktop/SPEER/Metadaten")

#meta_summary <- summary %>% 
  #group_by(ID) %>% 
  #summarise(enddate = max(datetime_max), 
            #startdate = min(datetime_min)
  #)
#write.table(meta_summary, file = "meta_summary.csv", row.names = FALSE, sep = ";")

#meta_summary_temp <- summary %>% 
  #group_by(ID_temp) %>% 
  #summarise(enddate = max(datetime_max), 
   #         startdate = min(datetime_min)
  #)
#write.table(meta_summary_temp, file = "meta_summary_temp.csv", row.names = FALSE, sep = ";")

#meta_summary_temp_press <- summary %>% 
  #group_by(ID_temp_press) %>% 
  #summarise(enddate = max(datetime_max), 
  #          startdate = min(datetime_min)
  #)
#write.table(meta_summary_temp_press, file = "meta_summary_temp_press.csv", row.names = FALSE, sep = ";")

#data$ID_temp <- paste(data$ID, data$ntemp,sep = "_")
#data$ID_temp_press <- paste(data$ID_temp, data$npress, sep = "_")

#meta_data <- data %>% 
  #group_by(ID) %>% 
  #summarise(enddate = max(datetime), 
            #startdate = min(datetime)
  #)        
#write.table(meta_data, file = "meta_data.csv", row.names = FALSE, sep = ";")

#meta_data_temp <- data %>% 
 # group_by(ID_temp) %>% 
 # summarise(enddate = max(datetime), 
  #          startdate = min(datetime)
 # ) 
#write.table(meta_data_temp, file = "meta_data_temp.csv", row.names = FALSE, sep = ";")

#meta_data_temp_press <- data %>% 
 # group_by(ID_temp_press) %>% 
 # summarise(enddate = max(datetime), 
  #          startdate = min(datetime)
 # ) 
#write.table(meta_data_temp_press, file = "meta_data_temp_press.csv", row.names = FALSE, sep = ";")

#add columns with headers Intercept/Slope and values (Intercept)/t -according to the identifier in column term in result_tidy- to result_glance
#result_glance[,"intercept"] <- "(Intercept)"
#result_glance[,"slope"] <- "m_time"

#remove columns indicating slope or intercept, they are in a crossed table format now
#result <- result[,! names(result) %in% c("slope","intercept"), drop = F]

#cleanse data of outliers
#data <- data[!(data$press < 7.9 & data$press > 1.1),]
#data <- data[!(data$press > 8.1 | data$press < 0.9),]
#data <- data[!(data$sat > 110 | data$sat < 80 | data$temp > 22 | data$temp < 10),]
#data <- data[!(data$temp > 13 & data$temp < 14 | data$temp > 16 & data$temp < 18 | data$temp > 20),]

#rename the columns in the data frame "result" to make them clear
#result <- result %>%
#rename(estimate_i = estimate.x,
#se_i = std.error.x,
#fstat_i = statistic.y,
#p_i = p.value.y,
#estimate_s = estimate.y,
#se_s = std.error.y,
#fstat_s = statistic,
#p_s = p.value,
#fstat = statistic.x,
#p = p.value.x,
#)

#add a suffix to coloumn names indicating whether it is max, min, mean or SD
#colnames(mins) <- (paste(colnames(mins), "min", sep ="_"))
#colnames(maxs) <- (paste(colnames(maxs), "max", sep ="_"))
#colnames(means) <- (paste(colnames(means), "mean", sep ="_"))
#colnames(sds) <- (paste(colnames(sds), "sd", sep ="_"))
#colnames(count) <- (paste(colnames(count), "count", sep ="_"))

#add a suffix to coloumn names indicating whether it is max, min, mean or SD
#colnames(mins)[2:ncol(mins)] <- (paste(colnames(mins)[2:ncol(mins)], "min", sep ="_"))
#colnames(maxs)[2:ncol(maxs)] <- (paste(colnames(maxs)[2:ncol(maxs)], "max", sep ="_"))
#colnames(means)[2:ncol(means)] <- (paste(colnames(means)[2:ncol(means)], "mean", sep ="_"))
#colnames(sds)[2:ncol(sds)] <- (paste(colnames(sds)[2:ncol(sds)], "sd", sep ="_"))
#colnames(counts)[2:ncol(counts)] <- (paste(colnames(counts)[2:ncol(counts)], "count", sep ="_"))

#nrow(as.matrix(which(summary$ntemp-summary$temp_max < -1.5 | summary$ntemp-summary$temp_min > 1.5)))

#mins <- aggregate(cbind(conc, sat, ntemp, temp_s, temp, npress, press_s, press, speed_s, speed, run, tunnel, ID) ~ no, data, min)
#maxs <- aggregate(cbind(m_time, conc, sat, temp_s, temp, press_s, press, speed_s, speed) ~ no, data, max)
#means <- aggregate(cbind(conc, sat, temp, press, speed)~ no, data, mean)
#sds <- aggregate(cbind(conc, sat, temp, press, speed)~ no, data, sd)
#counts <- aggregate(cbind(conc, temp) ~ no, data, function (x) sum(!is.na(x)))

#add a suffix to coloumn names indicating whether it is max, min, mean or SD
#colnames(mins)[2:ncol(mins)] <- (paste(colnames(mins)[2:ncol(mins)], "min", sep ="_"))
#colnames(maxs)[2:ncol(maxs)] <- (paste(colnames(maxs)[2:ncol(maxs)], "max", sep ="_"))
#colnames(means)[2:ncol(means)] <- (paste(colnames(means)[2:ncol(means)], "mean", sep ="_"))
#colnames(sds)[2:ncol(sds)] <- (paste(colnames(sds)[2:ncol(sds)], "sd", sep ="_"))
#colnames(counts)[2:ncol(counts)] <- (paste(colnames(counts)[2:ncol(counts)], "count", sep ="_"))

#Merging all the different data.frames created above: maxs, mins, means, sds, counts tidy and glance by "no" (essentially VLOOKUP in excel, solution for different col_names see below)
#result <- merge(result, maxs, by = "no", all = TRUE)
#result <- merge(result, mins, by = "no", all = TRUE)
#result <- merge(result, means, by = "no", all = TRUE)
#result <- merge(result, sds, by = "no", all = TRUE)
#result <- merge(result, counts, by = "no", all = TRUE)

#Merging all the different data.frames created above: maxs, mins etc, tidy and glance - essentially VLOOKUP in excel
#1st step: = if "no" and "Intercept" (in relult_glance) match the values in "no" and "term" (result_tidy), these parts of the data frames will be merged. Same for Slope, just with "no"/"slope" and "no"/"term".  
#result <- merge(result_glance, result_tidy_intercept , by = "no", all = TRUE)
#result <- merge(result, result_tidy_slope, by.x = c("no"), by.y = c("no"))
#result <- merge(result, maxs, by.x = "no", by.y = "no_max")
#result <- merge(result, mins, by.x = "no", by.y = "no_min")
#result <- merge(result, means, by.x = "no", by.y = "no_mean")
#result <- merge(result, sds, by.x = "no", by.y = "no_sd")

#rename columns where necessary (e.g. ID was captured in the "mins" but doesn't need the extention)
#summary <- summary %>%
#rename(ID = ID_min,
# tunnel = tunnel_min,
# npress = npress_min,
# ntemp = ntemp_min,
# run = run_min,
# n = conc_count
#)

#check if values in two columns match
#rownames(result[!result$temp_count %in% result$n, ])

#compare two vectors for unmatching values (order: A,B displays what is in A but not B)
#test1 <- as.vector(colnames(summary))
#test2 <- as.vector(col_order)
#setdiff(test2,test1)

#write.table(as.data.frame(colnames(summary)), file = "header.csv", row.names = FALSE, sep = ";")

#summary <- summary[!(summary$estimate_s >= 0 | summary$n <= 60 | summary$adj.r.squared_lm < 0.5 | summary$temp_max-summary$temp_min > 1| summary$press_min < 0.8 | summary$press_max >8.2 | summary$press_min > 1.2 & summary$press_min < 7.8 | summary$press_max > 1.2 & summary$press_max < 7.8 | summary$timediff_max >= 0.01666666667 | summary$press_max-summary$press_min > 0.2),]

#count_pressure_maxmindev_above_point2 <- nrow(as.matrix(which(summary$press_min < 0.8 | summary$press_max > 8.2 | summary$press_min > 1.2 & summary$press_min < 7.8 | summary$press_max > 1.2 & summary$press_max < 7.8)))
#pressure_maxmindev_above_point2 <- summary[(summary$press_min < 0.8 | summary$press_max >8.2 | summary$press_min > 1.2 & summary$press_min < 7.8 | summary$press_max > 1.2 & summary$press_max < 7.8),]
#write.table(pressure_maxmindev_above_point2, file = "pressure_maxmindev_above_point2.csv", row.names = FALSE, sep = ";")

####################################unnessecary since 0 cases##############################
#count_time <- nrow(summary[(summary$datetime_max-summary$datetime_min < 0.5),])
#time <- summary[(summary$datetime_max-summary$datetime_min < 0.5),]
#time$condition <- "measurement less than 30min"
#write.table(time, file = "time.csv", row.names = FALSE, sep = ";")

#all_removals <- merge(all_removals, conditions[[i]][, c("no", colnames(conditions[[i]])[49])], by="no", all.x=TRUE)

#all_removals <- do.call(rbind,lapply(list.files(pattern="*.csv"), fread, sep = ";"))

#for (i in 1:length(conditions)){
#all_removals <- merge(all_removals, conditions[[i]][, c("no", colnames(conditions[[i]])[49])], by="no", all.x=TRUE)

#test <- all_removals_full %>%
  #group_by(ID_temp_press, condition) %>%
  #summarise(n = n()) %>% 
  #spread(condition, n, fill = 0)
#}



#summarize all counts from the above in one matrix and write table

#count_r_squared <- nrow(filter(summary, adj.r.squared_lm < 0.9))
#count_slope <- nrow(filter(summary, estimate_s >= 0))
#count_n <- nrow(filter(summary, n <= 60))
#count_distance <- nrow(filter(summary, timediff_max >= 0.01666666667))
#count_temp <- nrow(as.matrix(which(summary$temp_max-summary$temp_min > 1)))
#count_press <- nrow(as.matrix(which(summary$press_mean-summary$npress > 0.1 | summary$press_mean-summary$npress < -0.1)))
#count_press_s <- nrow(as.matrix(which(summary$press_s_max-summary$press_s_min > 1)))
#count_resp <- nrow(as.matrix(which(summary$estimate_s/mean(summary$estimate_s) > 50)))
#count_sat <- nrow(filter(summary, sat_max < 90))

#counts <- list(count_r_squared, count_slope, count_n, count_distance, count_temp, count_press, count_press_s, count_resp, count_sat)

#count_of_removed <- rbind(count_r_squared, count_slope, count_n, count_distance, count_temp, count_press, count_press_s, count_resp, count_sat)
#count_of_removed <- as.data.frame(cbind(Condition = rownames(count_of_removed), count_of_removed), row.names = FALSE) #add column with criterion (from index) and transform to dataframe
#count_of_removed <- count_of_removed %>% rename(Count = V2) # edit column headers
#write.table(count_of_removed, file = "C:/Users/pohlmann/Desktop/SPEER/Ergebnisse/Bereinigung/count_of_removed.csv", row.names = FALSE, sep = ";") #write table