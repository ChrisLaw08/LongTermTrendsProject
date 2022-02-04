library(lubridate)
library(tidyverse)
library(rstatix)
library(readxl)
library(tibbletime)
library(ggpubr)

ALSCLWC<-read_excel("WFC.2017.Met.R.xlsx", sheet = "2017 WFC METEOROLOGICAL DATA", skip = 7 , col_names = TRUE)
ALSCLWC<-filter(ALSCLWC, month(Date_Time) == 6| month(Date_Time) == 7)
#ALSCLWC = ALSCLWC[-1,]

MinuteData<-vroom("Minute2017MetData.csv", col_names = TRUE)

Minute2017Data<-MinuteData%>%
  mutate(date = as.POSIXct(paste(Date, Time), format = "%d/%m/%Y %H:%M"))
Minute2017Data<-as_tbl_time(MinuteData, index = date)

ALSCLWC<-ALSCLWC%>%
  mutate(date = Date_Time)

HourData<-MinuteData%>%
  collapse_by(period = "1 hourly")%>%
  group_by(date)%>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))%>%
  mutate(date = date +minutes(1))



ggplot(HourData)+
  geom_point(aes(x = date, y = LWC...Value, color = "Envidas"))+
  geom_point(aes(x = date, y = `LWC (g/m3)`, color = "ALSC"), data = ALSCLWC)


HourData$date<-as.character(HourData$date)
ALSCLWC$date<-as.character(ALSCLWC$date)

AllLWCData<-left_join(ALSCLWC, HourData, by = "date")
AllLWCData<-AllLWCData%>%
  mutate(LWCALSC = `LWC (g/m3)`, LWCAWI = LWC...Value)





ggplot(subset(AllLWCData), aes(x = LWCALSC, y = LWCAWI ))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y~x+0
  )

JustLWC<-select(AllLWCData, c("date","LWCALSC", "LWCAWI"))

##New Calculations to see if weighting by cloud time works ####
#Start with just getting the averagers with minute data with Confirm > 0.25

HourDataAdj<-MinuteData%>%
  collapse_by(period = "1 hourly")%>%
  filter(`CONFIRM : Value` > 0.25)%>%
  group_by(date)%>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))%>%
  mutate(date = date +minutes(1))

HourDataAprox<-MinuteData%>%
  collapse_by(period = "1 hourly")%>%
  group_by(date)%>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))%>%
  mutate(date = date +minutes(1))%>%
  filter(`CONFIRM : Value` > 0.25)%>%
  mutate(LWCConfirm = `LWC : Value`/`CONFIRM : Value`)


HourDataAdj$date<-as.character(HourDataAdj$date)
HourDataAprox$date<-as.character(HourDataAprox$date)


HourDataJoin<-inner_join(HourDataAdj,HourDataAprox, by = "date")

HourDataJoin<-HourDataJoin%>%
  mutate(LWCAdj = `LWC : Status.x`, LWCAprox = `LWC : Status.y`)


ggplot(HourDataJoin, aes(x = LWCAdj, y = LWCAprox ))+
  geom_point()+
  geom_smooth(method = lm)+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y~x
  )



####Now we want to rerun the LWC calculations and see how they do. ####

AllCloudData<-read.csv("Alldata.csv", header = TRUE)
LWCData3Hour<-AllCloudData%>%
  filter(Year > 2008 &  Year < 2014)%>%
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y %H:%M"))

LWCData12Hour<-AllCloudData%>%
  filter(Year >= 2014 &  Year < 2021)%>%
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y %H:%M"))



Hour3Data<-read_csv("3HourData.csv", col_names = TRUE)
Hour3Data<-Hour3Data%>%
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y %H:%M"))%>%
  mutate(LabID = LABNO)
#Hour3Data<-as_tbl_time(Hour3Data, index = date)

Hour12Data<-read_csv("WFC12HourMet.csv", col_names = TRUE)
Hour12Data<-Hour12Data%>%
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y %H:%M"))
Hour12Data<-as_tbl_time(Hour12Data, index = date)

Hour3Calc<-Hour3Data%>%
 # mutate(LWCCalc = LWC/(CLOUD_HR/60))%>%
 # filter(!is.infinite(LWCCalc))%>%
 # collapse_by(period = "3 hour")%>%
  group_by(LabID)%>%
  filter(CLOUD_HR/60 > 0.25)%>%
  mutate(LWCCalc = LWC/(CLOUD_HR/60))%>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))
  #mutate(date = date - hours(2))


ggplot(Hour3Data)+
  geom_density(aes(x =  RAIN_HR+CLOUD_HR))

LWC3Hour<-select(Hour3Calc, c("LabID","LWCCalc", "Year", "RAIN_HR", "PSA", "WSP", "WDR", "Temperature"))

Hour12Calc<-Hour12Data%>%
 # filter(LWC < 2.0)%>%
  collapse_by("12 hourly", start_date = "2014-06-10 18:00:00")%>%
  group_by(date)%>%
  filter(CloudCount/60 > 0.25 & RainCount/60 < 0.25)%>%
  mutate(LWCCalc = LWC/(CloudCount/60))%>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE))%>%
  mutate(date = date+hours(1))%>%
  mutate(RAIN_HR=RainCount)

LWC12Hour<-select(Hour12Calc, c('date', "LWCCalc", "Year", "RAIN_HR", "PSA", "WSP", "WDR", "Temperature"))


LWCData3HourAll<-inner_join(LWCData3Hour, LWC3Hour, by = "LabID")

LWCData12Hour$date<-as.character(LWCData12Hour$date)
LWC12Hour$date<-as.character(LWC12Hour$date)
LWCData12HourAll<-inner_join(LWCData12Hour, LWC12Hour, by = "date")

LWCCalcAll<-rbind(LWCData3HourAll,LWCData12HourAll)




#LWCCalcAll$date<-as.character(LWCCalcAll$date)
#LWCData$date<-as.character(LWCData$date)

NaLWC<-filter(LWCCalcAll, is.na(LWCCalc))
IncludedLWC<-filter(LWCCalcAll, !is.na(LWCCalc))

LWCDataAll%>%
  mutate(NaLWC = ifelse(is.na(LWCCalc), "Yes", "No"))%>%
  ggplot()+
  geom_boxplot(aes(x = "NH4", y = NH4, fill = NaLWC))+
  geom_boxplot(aes(x = "SO4", y = SO4, fill = NaLWC))+
  geom_boxplot(aes(x = "NO3", y = NO3, fill = NaLWC))


ggplot(LWCCalcAll)+
  geom_point(aes(x = LWC, y = LWCCalc, color = as.factor(Year.x)))+
  geom_abline(intercept = 0, slope = 1)


LWCDataAll%>%
  mutate(NaLWC = ifelse(is.na(LWCCalc), "Yes", "No"))%>%
  wilcox_test(WSOC~NaLWC)
  


ggplot(LWCCalcAll)+
  geom_boxplot(aes(x = "SO4Calc", y = NH4, group = "SO4", fill = "LWC Calc Data"))+
  geom_boxplot(aes(x = "SO4All", y = NH4, group = "SO4", fill = "All Cloud Data"), data = AllCloudData)
  
LWCCalcTest<-select(LWCCalcAll, c(date,CA:WSOC))
LWCCalcAll$date<-as.character(LWCCalcAll$date)
LWCAllTest<-select(LWCDataAll, c(date,CA:WSOC))
LWCAllTest$date<-as.character(LWCAllTest$date)

TTestAll<-right_join(LWCCalcTest,LWCAllTest, by = 'date')


TTestAll%>%
  


LWCLoadingData<-LWCCalcAll%>%
  filter(LWCCalc < 2 & Year.x < 2021)%>%
  mutate(Year = Year.x)%>%
  select(-c(Year.x, Year.y))%>%
  group_by(Year)%>%
  #mutate(TOCCorrected = case_when(Year.x == 2018 | Year == 2019 ~ WSOC*(1/.8486), 
            #                      TRUE~WSOC))%>%
  mutate(NH4Mass = NH4*LWCCalc, SO4Mass = SO4*LWCCalc, NO3Mass= NO3*LWCCalc, CAMass =  CA*LWCCalc, WSOCMass =WSOC*LWCCalc, KMass = K*LWCCalc, TOCMass = WSOC*LWCCalc,
         CLMass = CL*LWCCalc, NaMass = Sodium*LWCCalc, MGMass = MG*LWCCalc)
write.csv(LWCLoadingData, "CloudLWCLoading.csv")

ggplot(subset(LWCCalcAll))+
  geom_point(aes(x= LWC, y = WSOC, color= as.factor(Year.x)))+
  facet_wrap(~Year.x)+xlim(c(0,3))




###This function can calculate fraction of cateogorical variables (freaking wierd)
LWCDataAll%>%
  mutate(NaLWC = ifelse(is.na(LWCCalc), "Yes", "No"))%>%
  group_by(Year.x)%>%
  summarise(pct.Yes = mean(NaLWC =="Yes"))%>%
  ggplot()+
  geom_line(aes(x = Year.x, y = pct.Yes, color = "Percent of Samples Removed"))+
  geom_point(aes(x = Year.x, y = pct.Yes, color = "Percent of Samples Removed"))
  #geom_line(aes(x = Year.x, y = NH4, color = "NH4"))+
  #geom_point(aes(x = Year.x, y = NH4, color = "NH4"))+
  #geom_line(aes(x = Year.x, y = NO3, color = "NO3"))+
  #geom_point(aes(x = Year.x, y = NO3, color = "NO3"))+
  facet_wrap(~NaLWC)
  
LWCDataAll%>%
  mutate(NaLWC = ifelse(is.na(LWCCalc), "Yes", "No"))%>%
  group_by(Year.x, NaLWC)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  filter(NaLWC == "No")%>%
  ggplot()+
  geom_line(aes(x = Year.x, y = WSOC, color = "SO4"))+
  geom_point(aes(x = Year.x, y = WSOC, color = "SO4"))+
 # geom_line(aes(x = Year.x, y = NH4, color = "NH4"))+
 # geom_point(aes(x = Year.x, y = NH4, color = "NH4"))+
 # geom_line(aes(x = Year.x, y = NO3, color = "NO3"))+
#  geom_point(aes(x = Year.x, y = NO3, color = "NO3"))+
  geom_line(aes(x = Year, y = WSOC, color = "SO4 No Filter"), data = LWCDataSummary)+
  geom_point(aes(x = Year, y = WSOC, color = "SO4 No Filter"),data = LWCDataSummary)+
#  geom_line(aes(x = Year, y = NH4, color = "NH4 No Filter"), data = LWCDataSummary)+
 # geom_point(aes(x = Year, y = NH4, color = "NH4 No Filter"),data = LWCDataSummary)+
#  geom_line(aes(x = Year, y = NO3, color = "NO3 No Filter"), data = LWCDataSummary)+
#  geom_point(aes(x = Year, y = NO3, color = "NO3 No Filter"),data = LWCDataSummary)
 # facet_wrap(~NaLWC)  

  


LWCDataSummary<-LWCData%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))

state <- rep(c(rep("Idaho", 10), rep("Maine", 10)), 2)
student.id <- sample(1:1000,8,replace=T)
gender <- rep( c("Male","Female"), 100*c(0.25,0.75) )  
gender <- sample(gender, 40)
school.data <- data.frame(student.id, state, gender)

school.data %>% 
  group_by(state) %>%
  mutate(pct.female = mean(gender == "Female"))








ggplot(LWCDataAll, aes(x= LWC, y = LWCCalc))+
  geom_point(aes(x= LWC, y = LWCCalc, color = RAIN_HR))+
  geom_smooth(method = lm, formula = y~x+0)+
  scale_color_viridis_c(limits = c(0,15))+
  stat_regline_equation(
    aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    formula = y~x+0)
  


LWCDataAll%>%
  mutate(TOCMass = WSOC*LWCCalc, NH4Mass = NH4*LWCCalc)%>%
  group_by(Year.x)%>%
  summarise(across(where(is.numeric), median, na.rm=TRUE))%>%
  ggplot()+
  geom_line(aes(x = Year.x, y = Temperature, color = "Year Avg"), size = 3)+
  geom_point(aes(x= Year.x, y = Temperature, color = "Year Avg"), size= 4)
  
  
sdf

ggplot(subset(LWCDataAll, Year.x < 2017 & CA > 0 & SO4 > 0),
       aes(x= LWCCalc/PSA, y = NH4))+
  geom_point(aes(x= LWCCalc/PSA, y = NH4, color = Year.x))+scale_color_gradientn(colors = rainbow(6))+
  geom_smooth(aes(x = LWCCalc/PSA, y = NH4), method = "lm", formula = y~exp(-x))+
  xlim(c(0,0.002))+scale_x_log10()+ylim(c(0,1500))
  stat_regline_equation()
  
  

geom_smooth(aes(x= LWCCalc, y = LWCCalc))


