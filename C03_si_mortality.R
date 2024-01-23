#-------------------------------------------------------------------------#
#--------------|             Mortality Analysis     |---------------------#
#--------------|             Suresh Parit           |---------------------#
#-------------------------------------------------------------------------#

#--------------|             1. Libraries     |---------------------#
#install.packages('demography')
#install.packages('tidyverse')
#install.packages('StMoMo')

library('demography')
library(ggplot2)
library(demography)
library(tidyverse)
library(StMoMo)
library(fds)

#--------------|             2. Data        |---------------------#

SPdata <- hmd.mx(country = "ESP", username = "suraj2611asi@gmail.com",
                  password = "Hmd@2024", label = "SPAIN")

SPdata$type
SPdata$label
SPdata$lambda
SPdata$year
SPdata$pop$female
SPdata$pop$male
SPdata$pop$total

SPdata$rate$female
#View(SPdata$rate$male)
SPdata$rate$total

SPdata$pop$female
#View(SPdata$rate$male)


#--------------|             3.  StMoMO : Stochastic Mortality Model        |---------------------#

male_data <- StMoMoData(SPdata,series = "male") 
female_data <- StMoMoData(SPdata,series = "female") 
persp(male_data$ages,male_data$years,male_data$Dxt/male_data$Ext,theta=-145
      ,phi=15,xlab='ages',ylab='years',zlab='rates',col="lightgreen",expand=0.5,shade=0.8,ticktype="detailed")
persp(female_data$ages,female_data$years,female_data$Dxt/female_data$Ext,theta=-145
      ,phi=15,xlab='ages',ylab='years',zlab='rates',col="lightgreen",expand=0.5,shade=0.8,ticktype="detailed")



actual<-extract.years(SPdata,1990:2021)
actual_male <- StMoMoData(actual,series = "male") # male series data 
actual_female <- StMoMoData(actual,series = "female") #  female series data


SPdata1<-extract.years(SPdata,1990:2009)

SP_male <- StMoMoData(SPdata1,series = "male") # male series data 
SP_female <- StMoMoData(SPdata1,series = "female") #  female series data

SP_male$type
summary(SP_male)

SP_female$type
summary(SP_female)


## central2initial : Transform StMoMoData from central to initial exposures. 
## initial exposures are approximated by transforming the available central exposures.
## Initial exposures are computed by adding one half of the deaths to the central exposures.


SP_male_1 = central2initial(SP_male)
summary(SP_male_1)

SP_female_1 = central2initial(SP_female)
summary(SP_female_1)


#================== Models

LC <- lc(link = "logit")                     #  Lee-Carter models
RH <- rh(link = "logit", cohortAgeFun = "1") # Renshaw and Haberman model: Lee-Carter with cohort effects
# APC <- apc(link = "logit")                 #  Age-Period-Cohort (APC)
CBD <- cbd()                                 #  The Cairns-Blake-Dowd model : no static age function and no cohort effect.
M7 <- m7()                                   # Quadratic CBD model with cohort effects

# PLAT  combines the CBD model with some features of the Lee-Carter model to produce
# a model that is suitable for full age ranges and captures the cohort effect

f2 <- function(x, ages) mean(ages) - x

constPlat <- function(ax, bx, kt, b0x, gc, wxt, ages){
 nYears <- dim(wxt)[2]
 x <- ages
 t <- 1:nYears
 c <- (1 - tail(ages, 1)):(nYears - ages[1])
 xbar <- mean(x)
 phiReg <- lm(gc ~ 1 + c + I(c ^ 2), na.action = na.omit)
 phi <- coef(phiReg)
 gc <- gc - phi[1] - phi[2] * c - phi[3] * c ^ 2
 kt[2, ] <- kt[2, ] + 2 * phi[3] * t
 kt[1, ] <- kt[1, ] + phi[2] * t + phi[3] * (t ^ 2 - 2 * xbar * t)
 ax <- ax + phi[1] - phi[2] * x + phi[3] * x ^ 2
 ci <- rowMeans(kt, na.rm = TRUE)
 ax <- ax + ci[1] + ci[2] * (xbar - x)
 kt[1, ] <- kt[1, ] - ci[1]
 kt[2, ] <- kt[2, ] - ci[2]
 list(ax = ax, bx = bx, kt = kt, b0x = b0x, gc = gc)
 }

PLAT <- StMoMo(link = "logit", staticAgeFun = TRUE, periodAgeFun = c("1", f2), cohortAgeFun = "1", constFun = constPlat)

#=========== Formulas of Models

library(gnm)
LC$gnmFormula
RH$gnmFormula
CBD$gnmFormula
M7$gnmFormula
PLAT$gnmFormula

#==================================================================

ages.fit <- 45:95 #  CBD model and the M7 model have been particularly designed to fit higher ages.

#==========================================================================================
#------------------------- | For Male
#==========================================================================================
# The first and last three cohorts years are excluded from the fitting via the argument
# wxt. The appropriate 0-1 weighting matrix, wxt, is constructed using the utility function genWeightMat

wxt <- genWeightMat(ages = ages.fit, years = SP_male_1$years,clip = 3)
LCfit_M <- fit(LC, data = SP_male_1, ages.fit = ages.fit, wxt = wxt)
RHfit_M <- fit(RH, data = SP_male_1, ages.fit = ages.fit, wxt = wxt,start.ax = LCfit_M$ax, start.bx = LCfit_M$bx, start.kt = LCfit_M$kt)
#APCfit_M <- fit(APC, data = SP_male_1, ages.fit = ages.fit, wxt = wxt)
CBDfit_M <- fit(CBD, data = SP_male_1, ages.fit = ages.fit, wxt = wxt)
M7fit_M <- fit(M7, data = SP_male_1, ages.fit = ages.fit, wxt = wxt)
PLATfit_M <- fit(PLAT, data = SP_male_1, ages.fit = ages.fit, wxt = wxt)


result_M<-data.frame(model=c('LC','RH','CBD','M7','PLAT'),
                   'AIC'=c(AIC(LCfit_M),AIC(RHfit_M),AIC(CBDfit_M),AIC(M7fit_M),AIC(PLATfit_M)),
                   'BIC'=c(BIC(LCfit_M),BIC(RHfit_M),BIC(CBDfit_M),BIC(M7fit_M),BIC(PLATfit_M)),
                   'No_Par'=c(LCfit_M$npar,RHfit_M$npar,CBDfit_M$npar,M7fit_M$npar,PLATfit_M$npar))

# Rank the result by the 'BIC' 
result_M <- result_M[order(result_M$BIC), ]
result_M$Rank <- rank(result_M$BIC, ties.method = "min")
result_M
# 
# model      AIC      BIC No_Par Rank
# 5  PLAT 11080.62 11817.98    150    1
# 2    RH 11043.08 11942.65    183    2
# 4    M7 11464.19 12058.99    121    3
# 1    LC 11564.95 12154.84    120    4
# 3   CBD 26858.75 27055.38     40    5


par(mar = c(5, 4, 4, 2) + 0.1) 

plot(LCfit_M, nCol = 3)
plot(RHfit_M, parametricbx = FALSE,nCol = 3)
plot(PLATfit_M, parametricbx = FALSE,nCol = 3)

LCres <- residuals(LCfit_M)
RHres <- residuals(RHfit_M)
PLATres <- residuals(PLATfit_M)

plot(LCres, type = "colourmap", reslim = c(-3.5, 3.5))
plot(RHres, type = "colourmap", reslim = c(-3.5, 3.5))
plot(PLATres, type = "colourmap", reslim = c(-3.5, 3.5))

plot(LCres, type = "scatter", reslim = c(-3.5, 3.5))
plot(RHres, type = "scatter", reslim = c(-3.5, 3.5))
plot(PLATres, type = "scatter", reslim = c(-3.5, 3.5))


LCfor_M <- forecast(LCfit_M, h = 21)
RHfor_M <- forecast(RHfit_M, h = 21,gc.order = c(1, 1, 0))
PLATfor_M <- forecast(PLATfit_M, h = 21, gc.order = c(2, 0, 0))


actual_M_mxt<-actual$rate$male
LCfor_M_mxt<-LCfor_M$rates
RHfor_M_mxt<-RHfor_M$rates
PLATfor_M_mxt<-PLATfor_M$rates

x <- c( "55","65","75")

matplot(actual$year, t(actual_M_mxt[x, ]), xlim = range(actual$year,LCfor_M$years), ylim = range(actual_M_mxt[x,],LCfor_M_mxt[x, ],RHfor_M_mxt[x,],PLATfor_M_mxt[x,]), type = "p", xlab = "years",
        ylab = "mortality rates (log scale)", log = "y", pch = 20, col = "maroon")

matlines(LCfor_M$year, t(LCfor_M_mxt[x, ]), type = "l", lty = 2, lwd = 1.5,col = "red")
matlines(RHfor_M$year, t(RHfor_M_mxt[x, ]), type = "l", lty = 2, lwd = 1.5,col = "blue")
matlines(PLATfor_M$year, t(PLATfor_M_mxt[x, ]), type = "l", lty = 2, lwd = 1.5,col = "darkgray")
# Adding legend
legend("topright", legend = c("LC", "RH","PLAT"), col = c("red", "blue","darkgray"), lty = 1, lwd = 1, cex = 0.4)





#==========================================================================================
#------------------------- | For Female
#==========================================================================================

wxt_f <- genWeightMat(ages = ages.fit, years = SP_female_1$years,clip = 3)
LCfit_F <- fit(LC, data = SP_female_1, ages.fit = ages.fit, wxt = wxt_f)
RHfit_F <- fit(RH, data = SP_female_1, ages.fit = ages.fit, wxt = wxt_f,start.ax = LCfit_F$ax, start.bx = LCfit_F$bx, start.kt = LCfit_F$kt)
CBDfit_F <- fit(CBD, data = SP_female_1, ages.fit = ages.fit, wxt = wxt_f)
M7fit_F <- fit(M7, data = SP_female_1, ages.fit = ages.fit, wxt = wxt_f)
PLATfit_F <- fit(PLAT, data = SP_female_1, ages.fit = ages.fit, wxt = wxt_f)

result_F<-data.frame(model=c('LC','RH','CBD','M7','PLAT'),
                     'AIC'=c(AIC(LCfit_F),AIC(RHfit_F),AIC(CBDfit_F),AIC(M7fit_F),AIC(PLATfit_F)),
                     'BIC'=c(BIC(LCfit_F),BIC(RHfit_F),BIC(CBDfit_F),BIC(M7fit_F),BIC(PLATfit_F)),
                     'No_Par'=c(LCfit_F$npar,RHfit_F$npar,CBDfit_F$npar,M7fit_F$npar,PLATfit_F$npar))

# Rank the result by the 'BIC' 
result_F <- result_F[order(result_F$BIC), ]
result_F$Rank <- rank(result_F$BIC, ties.method = "min")
result_F

# model      AIC      BIC No_Par Rank
# 5  PLAT 10765.54 11502.90    150    1
# 2    RH 10679.10 11578.68    183    2
# 1    LC 11288.14 11878.02    120    3
# 4    M7 12556.25 13151.05    121    4
# 3   CBD 59543.79 59740.42     40    5

# RH , PLAT and LC performing better 

plot(LCfit_F, nCol = 3)
plot(RHfit_F, parametricbx = FALSE,nCol = 3)
plot(PLATfit_F, parametricbx = FALSE,nCol = 3)

LCre_F <- residuals(LCfit_F)
RHre_F <- residuals(RHfit_F)
PLATre_F <- residuals(PLATfit_F)

plot(LCre_F, type = "colourmap", reslim = c(-3.5, 3.5))
plot(RHre_F, type = "colourmap", reslim = c(-3.5, 3.5))
plot(PLATre_F, type = "colourmap", reslim = c(-3.5, 3.5))

plot(LCre_F, type = "scatter", reslim = c(-3.5, 3.5))
plot(RHre_F, type = "scatter", reslim = c(-3.5, 3.5))
plot(PLATre_F, type = "scatter", reslim = c(-3.5, 3.5))


LCfor_F <- forecast(LCfit_F, h = 21)
RHfor_F <- forecast(RHfit_F, h = 21,gc.order = c(1, 1, 0))
PLATfor_F <- forecast(PLATfit_F, h = 21, gc.order = c(2, 0, 0))


actual_F_Fxt<-actual$rate$female
LCfor_F_Fxt<-LCfor_F$rates
RHfor_F_Fxt<-RHfor_F$rates
PLATfor_F_Fxt<-PLATfor_F$rates

x <- c( "55", "60","65","70","75")

matplot(actual$year, t(actual_F_Fxt[x, ]), xlim = range(actual$year,LCfor_F$years), ylim = range(actual_F_Fxt[x,],LCfor_F_Fxt[x, ],RHfor_F_Fxt[x,],PLATfor_F_Fxt[x,]), type = "p", xlab = "years",
        ylab = "mortality rates (log scale)", log = "y", pch = 20, col = "maroon")

matlines(LCfor_F$year, t(LCfor_F_Fxt[x, ]), type = "l", lty = 2, lwd = 1.5,col = "red")
matlines(RHfor_F$year, t(RHfor_F_Fxt[x, ]), type = "l", lty = 2, lwd = 1.5,col = "blue")
matlines(PLATfor_F$year, t(PLATfor_F_Fxt[x, ]), type = "l", lty = 2, lwd = 1.5,col = "darkgray")



#==============================| premiums for Male between 60-70


si_df<-read.csv("D:/SP/GitHub/P1/spanis_ insurance.csv")
dim(si_df)

View(si_df)
df_m<-si_df[si_df$Age >=60 & si_df$Age <=70 & si_df$Gender=='M',]
dim(df_m)
df_m$age_2010<-df_m$Age_Actuarial+1
df_m$age_2015<-df_m$Age_Actuarial+6
df_m$age_2020<-df_m$Age_Actuarial+11


# Desired age: 30.5 years
# Mortality rate at 30 years: 0.001
# Mortality rate at 31 years: 0.002
# Fraction of the year: 0.5
# mortality_fractional = 0.001 + (0.5 * (0.002 - 0.001))
# mortality_fractional = 0.0015

yr=c('2010','2015','2020')


create_df<-function(z,name){
  
  z1<-data.frame(z)
  names(z1) <- paste0(name, names(z1))
  z1$age=as.integer(rownames(z1))
  return(z1)
}

actual_M_df<-create_df(actual_M_mxt[,c('2009',yr)],'ACTUAL_')
LCfor_M_df<-create_df(LCfor_M_mxt[,yr],'LC_')
RHfor_M_df<-create_df(RHfor_M_mxt[,yr],'RH_')
PLATfor_M_df<-create_df(PLATfor_M_mxt[,yr],'PLAT_')


df_m<-merge(df_m,actual_M_df,by.x = 'Age_Actuarial', by.y = 'age',all.x = TRUE)
df_m<-merge(df_m,LCfor_M_df,by.x = 'Age_Actuarial', by.y = 'age',all.x = TRUE)
df_m<-merge(df_m,RHfor_M_df,by.x = 'Age_Actuarial', by.y = 'age',all.x = TRUE)
df_m<-merge(df_m,PLATfor_M_df,by.x = 'Age_Actuarial', by.y = 'age',all.x = TRUE)

View(df_m)


cal_p<-function(m,C){
  
  i=0.03 # interest rate 
  q=m/(1+0.5*m)
  A=(q)/(1+i)
  P=C*A
  return(P)
}


df_m$actual_P<-cal_p(df_m$ACTUAL_X2010,df_m$Capital)
df_m$LC_P<-cal_p(df_m$LC_X2010,df_m$Capital)
df_m$RH_P<-cal_p(df_m$RH_X2010,df_m$Capital)
df_m$PLAT_P<-cal_p(df_m$PLAT_X2010,df_m$Capital)


df_m$actual_P15<-cal_p(df_m$ACTUAL_X2015,df_m$Capital)
df_m$LC_P15<-cal_p(df_m$LC_X2015,df_m$Capital)
df_m$RH_P15<-cal_p(df_m$RH_X2015,df_m$Capital)
df_m$PLAT_P15<-cal_p(df_m$PLAT_X2015,df_m$Capital)

df_m$actual_P20<-cal_p(df_m$ACTUAL_X2020,df_m$Capital)
df_m$LC_P20<-cal_p(df_m$LC_X2010,df_m$Capital)
df_m$RH_P20<-cal_p(df_m$RH_X2020,df_m$Capital)
df_m$PLAT_P20<-cal_p(df_m$PLAT_X2010,df_m$Capital)

selected_columns<-c('actual_P','LC_P','RH_P','PLAT_P',
                    'actual_P15','LC_P15','RH_P15','PLAT_P15',
                    'actual_P20','LC_P20','RH_P20','PLAT_P20'
)

df<-df_m[,selected_columns]

# Calculate sum and median for each column
summary_df <- df %>%
  summarise(across(everything(), list(sum = sum, median = median)))

# Display the summary data frame
View(t(summary_df))



#=============================================================
# Mapping Mortality Rate to Premiums
#=============================================================

df_F<-si_df[si_df$Age >=60 & si_df$Age <=70 & si_df$Gender=='F',]
dim(df_F)
df_F$age_2010<-df_F$Age_Actuarial+1
df_F$age_2015<-df_F$Age_Actuarial+6
df_F$age_2020<-df_F$Age_Actuarial+11


yr=c('2010','2015','2020')

actual_F_df<-create_df(actual_F_Fxt[,c('2009',yr)],'ACTUAL_')
LCfor_F_df<-create_df(LCfor_F_Fxt[,yr],'LC_')
RHfor_F_df<-create_df(RHfor_F_Fxt[,yr],'RH_')
PLATfor_F_df<-create_df(PLATfor_F_Fxt[,yr],'PLAT_')


df_F<-merge(df_F,actual_F_df,by.x = 'Age_Actuarial', by.y = 'age',all.x = TRUE)
df_F<-merge(df_F,LCfor_F_df,by.x = 'Age_Actuarial', by.y = 'age',all.x = TRUE)
df_F<-merge(df_F,RHfor_F_df,by.x = 'Age_Actuarial', by.y = 'age',all.x = TRUE)
df_F<-merge(df_F,PLATfor_F_df,by.x = 'Age_Actuarial', by.y = 'age',all.x = TRUE)


df_F$actual_P<-cal_p(df_F$ACTUAL_X2010,df_F$Capital)
df_F$LC_P<-cal_p(df_F$LC_X2010,df_F$Capital)
df_F$RH_P<-cal_p(df_F$RH_X2010,df_F$Capital)
df_F$PLAT_P<-cal_p(df_F$PLAT_X2010,df_F$Capital)

summary(df_F)

df_F$actual_P15<-cal_p(df_F$ACTUAL_X2015,df_F$Capital)
df_F$LC_P15<-cal_p(df_F$LC_X2015,df_F$Capital)
df_F$RH_P15<-cal_p(df_F$RH_X2015,df_F$Capital)
df_F$PLAT_P15<-cal_p(df_F$PLAT_X2015,df_F$Capital)

df_F$actual_P20<-cal_p(df_F$ACTUAL_X2020,df_F$Capital)
df_F$LC_P20<-cal_p(df_F$LC_X2010,df_F$Capital)
df_F$RH_P20<-cal_p(df_F$RH_X2020,df_F$Capital)
df_F$PLAT_P20<-cal_p(df_F$PLAT_X2010,df_F$Capital)



df<-df_F[,selected_columns]

# Calculate sum and median for each column
summary_df <- df %>%
  summarise(across(everything(), list(sum = sum, median = median)))

# Display the summary data frame
View(t(summary_df))



#======================================================================================================



par(mfrow = c(1, 2))

matplot(actual$year, (actual_M_mxt["55", ]), xlim = range(actual$year,LCfor_M$years), ylim = range(actual_M_mxt["55",],LCfor_M_mxt["55", ],RHfor_M_mxt["55",],PLATfor_M_mxt["55",]), type = "p", xlab = "years",
        ylab = "mortality rates (log scale)", log = "y", pch = 20, col = "maroon")
# Adding legend
legend("topright", legend = c("LC", "RH","PLAT"), col = c("red", "blue","darkgray"), lty = 2, lwd = 1.5, cex = 0.4)

mtext("Male", side = 3, line = 1)
matlines(LCfor_M$year, (LCfor_M_mxt["55", ]), type = "l", lty = 2, lwd = 1.5,col = "red")
matlines(RHfor_M$year, (RHfor_M_mxt["55", ]), type = "l", lty = 2, lwd = 1.5,col = "blue")
matlines(PLATfor_M$year, (PLATfor_M_mxt["55", ]), type = "l", lty = 2, lwd = 1.5,col = "darkgray")


matplot(actual$year, (actual_F_Fxt["55", ]), xlim = range(actual$year,LCfor_F$years), ylim = range(actual_F_Fxt["55",],LCfor_F_Fxt["55", ],RHfor_F_Fxt["55",],PLATfor_F_Fxt["55",]), type = "p", xlab = "years",
        ylab = "mortality rates (log scale)", log = "y", pch = 20, col = "maroon")
mtext("Female", side = 3, line = 2)
matlines(LCfor_F$year, (LCfor_F_Fxt["55", ]), type = "l", lty = 2, lwd = 1.5,col = "red")
matlines(RHfor_F$year, (RHfor_F_Fxt["55", ]), type = "l", lty = 2, lwd = 1.5,col = "blue")
matlines(PLATfor_F$year, (PLATfor_F_Fxt["55", ]), type = "l", lty = 2, lwd = 1.5,col = "lightgreen")


par(mfrow = c(1, 2))

matplot(actual$year, (actual_M_mxt["65", ]), xlim = range(actual$year,LCfor_M$years), ylim = range(actual_M_mxt["65",],LCfor_M_mxt["65", ],RHfor_M_mxt["65",],PLATfor_M_mxt["65",]), type = "p", xlab = "years",
        ylab = "mortality rates (log scale)", log = "y", pch = 20, col = "maroon")
# Adding legend
legend("topright", legend = c("LC", "RH","PLAT"), col = c("red", "blue","darkgray"), lty = 2, lwd = 1.5, cex = 0.4)
mtext("Male", side = 3, line = 1)
matlines(LCfor_M$year, (LCfor_M_mxt["65", ]), type = "l", lty = 2, lwd = 1.5,col = "red")
matlines(RHfor_M$year, (RHfor_M_mxt["65", ]), type = "l", lty = 2, lwd = 1.5,col = "blue")
matlines(PLATfor_M$year, (PLATfor_M_mxt["65", ]), type = "l", lty = 2, lwd = 1.5,col = "darkgray")

matplot(actual$year, (actual_F_Fxt["65", ]), xlim = range(actual$year,LCfor_F$years), ylim = range(actual_F_Fxt["65",],LCfor_F_Fxt["65", ],RHfor_F_Fxt["65",],PLATfor_F_Fxt["65",]), type = "p", xlab = "years",
        ylab = "mortality rates (log scale)", log = "y", pch = 20, col = "maroon")
# Adding legend
legend("topright", legend = c("LC", "RH","PLAT"), col = c("red", "blue","darkgray"), lty = 2, lwd = 1.5, cex = 0.4)

mtext("Male", side = 3, line = 2)
matlines(LCfor_F$year, (LCfor_F_Fxt["65", ]), type = "l", lty = 2, lwd = 1.5,col = "red")
matlines(RHfor_F$year, (RHfor_F_Fxt["65", ]), type = "l", lty = 2, lwd = 1.5,col = "blue")
matlines(PLATfor_F$year, (PLATfor_F_Fxt["65", ]), type = "l", lty = 2, lwd = 1.5,col = "lightgreen")


par(mfrow = c(1, 2))

matplot(actual$year, (actual_M_mxt["75", ]), xlim = range(actual$year,LCfor_M$years), ylim = range(actual_M_mxt["75",],LCfor_M_mxt["75", ],RHfor_M_mxt["75",],PLATfor_M_mxt["75",]), type = "p", xlab = "years",
        ylab = "mortality rates (log scale)", log = "y", pch = 20, col = "maroon")
# Adding legend
legend("topright", legend = c("LC", "RH","PLAT"), col = c("red", "blue","darkgray"), lty = 2, lwd = 1.5, cex = 0.4)
mtext("Male", side = 3, line = 1)
matlines(LCfor_M$year, (LCfor_M_mxt["75", ]), type = "l", lty = 2, lwd = 1.5,col = "red")
matlines(RHfor_M$year, (RHfor_M_mxt["75", ]), type = "l", lty = 2, lwd = 1.5,col = "blue")
matlines(PLATfor_M$year, (PLATfor_M_mxt["75", ]), type = "l", lty = 2, lwd = 1.5,col = "darkgray")

matplot(actual$year, (actual_F_Fxt["75", ]), xlim = range(actual$year,LCfor_F$years), ylim = range(actual_F_Fxt["75",],LCfor_F_Fxt["75", ],RHfor_F_Fxt["75",],PLATfor_F_Fxt["75",]), type = "p", xlab = "years",
        ylab = "mortality rates (log scale)", log = "y", pch = 20, col = "maroon")
matlines(LCfor_F$year, (LCfor_F_Fxt["75", ]), type = "l", lty = 2, lwd = 1.5,col = "red")
matlines(RHfor_F$year, (RHfor_F_Fxt["75", ]), type = "l", lty = 2, lwd = 1.5,col = "blue")
matlines(PLATfor_F$year, (PLATfor_F_Fxt["75", ]), type = "l", lty = 2, lwd = 1.5,col = "lightgreen")

#==========================================================================


