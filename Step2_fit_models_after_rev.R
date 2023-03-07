#load data

#load("summary_after_rev.RData")
load("summary_first_rev.RData")

#use environment after running step 1

summary_backup <- summary
summary <- summary_backup

summary$ID_time <- as.numeric(summary$ID_time)
summary$run <-as.factor(summary$run)
summary$tunnel <-as.factor(summary$tunnel)
summary$ID <- as.factor(summary$ID)
summary$FR <- summary$front_true/(summary$front_false+summary$front_true)
summary$ID_press <- paste(summary$ID, summary$npress, sep ="_")


library(nlme)
library(lme4)
library(lmerTest)
library(cAIC4)
library(MuMIn)
library(dplyr)
library(lmtest)

#reduce dataset
summary_short <- summary[, names(summary) %in% c("length", "FR"," ntemp", "temp_mean", "press_mean", "ID_time", "speedBL", "run", "tunnel","ID", "resp", "npress"), drop = F]
summary_corr <- summary_backup[, names(summary_backup) %in% c("temp_mean", "press_mean", "speedBL", "run", "tunnel", "ID_time", "mass"), drop = F]
summary$mass_g <- summary$mass
summary$mass <- summary$mass_g/1000



########################################### SELECT RANDOM EFFECTS STRUCTURE ##########################################

lme.i <-  lme(log(resp) ~ temp_mean * press_mean + speedBL + run, 
              random = ~1|ID, 
              data = summary,
              method ="REML",
              #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
              correlation = corCAR1(form = ~ ID_time|ID))

lme.s <-  lme(log(resp) ~ temp_mean * press_mean + speedBL + run, 
            random = ~1+speedBL|ID, 
            data = summary,
            method ="REML",
            #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
            correlation = corCAR1(form = ~ ID_time|ID))

lme.t <-  lme(log(resp) ~ temp_mean * press_mean + speedBL + run, 
              random = ~1+temp_mean|ID, 
              data = summary,
              method ="REML",
              #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
              correlation = corCAR1(form = ~ ID_time|ID))

lme.p <-  lme(log(resp) ~ temp_mean * press_mean + speedBL + run, 
              random = ~1+press_mean|ID, 
              data = summary,
              method ="REML",
              #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
              correlation = corCAR1(form = ~ ID_time|ID))


lrtest(lme.i, lme.s, lme.t, lme.p)


summary(lme.s)
plot(lme.s)
qqnorm(resid(lme.s))
qqline(resid(lme.s))



######################################### SELECT FIXED EFFECTS STRUCTURE ##################################################

#fit best model by ML as quality check

lme.dredge <- lme(log(resp) ~ temp_mean * press_mean + speedBL + run, 
              random = ~1+speedBL|ID, #fill in the wanted RE structure
              data = summary,
              method ="ML",
              #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
              correlation = corCAR1(form = ~ ID_time|ID))

summary(lme.dredge)
plot(lme.dredge)
qqnorm(resid(lme.dredge))
qqline(resid(lme.dredge))

options(na.action = "na.fail")
select.lme <- dredge(lme.dredge, rank = "AIC", m.lim = c(3,NA), evaluate = T)
options(na.action = "na.omit")
setwd("../Ergebnisse")
write.table(select.lme, file = "select_lme.csv", row.names = FALSE, sep = ";")

#test if interaction is needed
lme.inter <- lme(log(resp) ~ temp_mean * press_mean + speedBL, 
                  random = ~1+speedBL|ID, #fill in the wanted RE structure
                  data = summary,
                  method ="ML",
                  #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
                  correlation = corCAR1(form = ~ ID_time|ID))

lme.nointer<- lme(log(resp) ~ temp_mean + press_mean + speedBL, 
                  random = ~1+speedBL|ID, #fill in the wanted RE structure
                  data = summary,
                  method ="ML",
                  #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
                  correlation = corCAR1(form = ~ ID_time|ID))


#compare whether interaction is significant
lrtest(lme.inter, lme.nointer)
anova(lme.inter, lme.nointer)

#test if run is needed
lme.run <- lme(log(resp) ~ temp_mean + press_mean + speedBL + run, 
                  random = ~1+speedBL|ID, #fill in the wanted RE structure
                  data = summary,
                  method ="ML",
                  #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
                  correlation = corCAR1(form = ~ ID_time|ID))

lme.norun <- lme(log(resp) ~ temp_mean + press_mean + speedBL , 
                 random = ~1+speedBL|ID, #fill in the wanted RE structure
                 data = summary,
                 method ="ML",
                 #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
                 correlation = corCAR1(form = ~ ID_time|ID))

summary(lme.run)
summary(lme.norun)
lrtest(lme.run, lme.norun)
AIC(lme.press, lme.nopress)



################################ GRAPH MAIN RESULTS ###########################################

lme.plot <- lme(log(resp) ~ temp_mean * press_mean + speedBL, 
                random = ~1+speedBL|ID, #fill in the wanted RE structure
                data = summary,
                method ="REML",
                #control = lmeControl(maxIter = 10000, msMaxIter = 100000, msTol = 0.005, tolerance = 0.005, pnlsTol = 0.005, pnlsMaxIter =1000, niterEM = 1000),
                correlation = corCAR1(form = ~ ID_time|ID))

summary(lme.plot)

ID_key <- data.frame(ID = unique(summary$ID),
                     Letter = letters[1:9])

summary <- summary %>% 
  left_join(ID_key, by = "ID") %>%
  rename(Letter_ID=Letter)

#----------------------Individual predictions plot------------------------#
{
  options(scipen = 999) #turn off scientific notations
  library(ggplot2)
  setwd("C:/Users/pohlmann/Desktop/SPEER/Ergebnisse")
  devtools::install_github("thomasp85/scico")
  library(scico)
  library(egg)
  library(nlme)
  
  #facet lable names
  tunnel.names <- c("tunnel 1", "tunnel 2", "tunnel 3")
  names(tunnel.names) <- c("1", "2", "3")
  
  run.names <- c("run 1", "run 2", "run 3")
  names(run.names) <- c("1", "2", "3")
  
  #define tag_facet with keeping strips
  tag_facet2 <- function(p, open = "(", close = ")", tag_pool = LETTERS, x = -Inf, y = Inf, 
                         hjust = -0.5, vjust = 1.5, fontface = 2, family = "", ...) {
    
    gb <- ggplot_build(p)
    lay <- gb$layout$layout
    tags <- cbind(lay, label = paste0(open, tag_pool[lay$PANEL], close), x = x, y = y)
    p + geom_text(data = tags, aes_string(x = "x", y = "y", label = "label"), ..., hjust = hjust, 
                  vjust = vjust, fontface = fontface, family = family, inherit.aes = FALSE)
  }

{
#-------------plot swimming speeds---------------#

  ggplot(summary, aes(speedBL))+
    geom_histogram(bins = 25)+
    facet_wrap(~ID_temp_press)

#----------------------Speed vs Resp-----------------------------#

  ggplot(summary, aes(x=speedBL, y=resp)) +
    geom_point(size=0.2)+ 
    theme_bw()+
    #geom_smooth(aes(linetype = as.factor(npress), colour = as.factor(ntemp)))+
    #geom_point(aes(x = (speed_s_min/(length/100)), colour = "red"))+
    #geom_point(aes(x = (speed_s_max/(length/100)), colour = "blue"))+
    #geom_vline(aes(xintercept = 0.6))+
    #geom_rect(aes(xmin=0.564, xmax = 0.667, ymin=-Inf, ymax=Inf),alpha = 0.004)+
    xlim(0.4, 0.8)+
    facet_grid(ntemp ~ ID_press)+
    guides(shape=FALSE)


#----------------------Speed vs Temp-----------------------------#

  ggplot(summary, aes(x=resp, y=speedBL)) +
    geom_point(size=0.2)+ 
    theme_bw()+
    #geom_smooth(aes(linetype = as.factor(npress), colour = as.factor(ntemp)))+
    #geom_point(aes(x = (speed_s_min/(length/100)), colour = "red"))+
    #geom_point(aes(x = (speed_s_max/(length/100)), colour = "blue"))+
    #geom_vline(aes(xintercept = 0.6))+
    #geom_rect(aes(xmin=0.564, xmax = 0.667, ymin=-Inf, ymax=Inf),alpha = 0.004)+
    #xlim(0.4, 0.8)+
    facet_wrap(~ID)+
    guides(shape=FALSE)
}
