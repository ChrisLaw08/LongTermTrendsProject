library(tidyverse)
library(ggpubr)
library(tidymodels)
library(rstatix)
library(tibbletime)
library(gridExtra)
library(lubridate)
library(openair)
library(trend)
library(mblm)
library(deming)
library(Kendall)
library(latex2exp)
library(gt)

#### Reading in the Data and make relavent calculations ####
AllCloudData<-read.csv("Alldata.csv", header = TRUE)
AllCloudData<-AllCloudData%>%
  dplyr::rename(TOC = WSOC)%>%
  filter(Year > 1993 & Year<2021)%>%
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y %H:%M"))%>%
  mutate(TOCCorrected = case_when(Year == 2018 | Year == 2019 ~ TOC*(1/.8486), 
                                  TRUE~TOC))%>%
  mutate(Cations = 1e6*10^(-LABPH)+CA+MG+K+Sodium+NH4,
         Anions = NO3+SO4+CL)%>%### Create Columns with Cations and Anions ueq/L
  mutate(RPD = 200*(Cations-Anions)/(Cations+Anions))%>% ## RPD Calculations
  mutate(Ratio = Cations/Anions)%>%
  mutate(HCond = (Hplus*349.5)/1000, Regime = ifelse(HCond/SPCOND < 0.35, "Non-Linear", "Linear"))%>%
  mutate(Class = ifelse((Cations < 100 & Anions < 100 & abs(RPD) < 100) | (Cations > 100 & abs(RPD) < 25) |(Anions > 100 & abs(RPD)<25), "Valid", "Invalid"))%>%
  mutate(Class = as.factor(Class))%>%
  mutate(Hplus = 10^(-LABPH))%>%
  mutate(SurplusNH4 = NH4-NO3-SO4)%>%
  mutate(HCO3Gas = (3.4*10^(-2)*410*10^(-6)*10^(-6.36))/Hplus)%>%##Gas phase HCO3
  mutate(HCO3Fraction = (Hplus*10^(-6.3))/(Hplus^2 + Hplus*10^(-6.36) + 10^(-6.36-10.36)))%>%##Ionizations Fraction of HCO3
  mutate(HCO3Total = (HCO3Fraction*(MG+CA)/2) + 1e6*HCO3Gas)## Both Gas Phase Additions 
  
##Cloud Water Loading Data:
LWCLoadingData<-read.csv("CloudLWCLoading.csv", header = TRUE)
LWCLoadingData<-LWCLoadingData%>%
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y %H:%M"))



ggplot(subset(AllCloudData, Year > 2009))+
  geom_point(aes(x= 1e6*10^(-LABPH), y = TOC*12/(SO4*96), color = SO4))+ylim(c(0,10))+scale_color_viridis_c(trans = "log10")


#### Median IQR Function Used for Trend Plots ####
median_IQR <- function(x) {
    data.frame(y = median(x),# Median
               ymin = quantile(x)[2], # 1st quartile
               ymax = quantile(x)[4])  # 3rd quartile
}

#### Valid Plot for Trend Analysis: ####
ValidpH<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Valid"))+
  stat_summary(aes(x = Year, group = Year, y =LABPH), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = LABPH), fun.y = "median",
               size = 3 ,color = "red", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = LABPH),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "red", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{pH}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('pH')

ValidSPCOND<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Valid"))+
  stat_summary(aes(x = Year, group = Year, y =SPCOND), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = SPCOND), fun.y = "median",
               size = 3 ,color = "blue", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = SPCOND),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "blue", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Conductivity ($\\mu S$ cm^{-1}$})"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('Conductivity')

ValidSO4<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Valid"))+
  stat_summary(aes(x = Year, group = Year, y =SO4), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = SO4), fun.y = "median",
               size = 3 ,color = "Yellow", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = SO4),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Yellow", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{SO4 Concentration ($\\mu eq$ L^{-1}$})"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('SO4')
ValidNO3<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Valid"))+
  stat_summary(aes(x = Year, group = Year, y = NO3), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = NO3), fun.y = "median",
               size = 3 ,color = "Purple", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NO3),fun ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Purple", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{NO3 Concentration ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('NO3')

ValidNH4<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Valid"))+
  stat_summary(aes(x = Year, group = Year, y = NH4), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = NH4), fun.y = "median",
               size = 3 ,color = "Orange", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NH4),fun ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Orange", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{NH4 Concentration ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('NH4')

ValidTOC<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Valid"))+
  stat_summary(aes(x = Year, group = Year, y = TOC), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = TOC), fun.y = "median",
               size = 3 ,color = "Forest Green", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = TOC),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Forest Green", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 15, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{TOC Concentration ($\\mu molC $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2018,2))+
  ggtitle('TOC')

grid.arrange(ValidpH,ValidSPCOND,ValidSO4,ValidNH4,ValidNO3,ValidTOC)
#### Invalid Plots ####
InvalidpH<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Invalid"))+
  stat_summary(aes(x = Year, group = Year, y =LABPH), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = LABPH), fun.y = "median",
               size = 3 ,color = "red", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = LABPH),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "red", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{pH}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('pH')

InvalidSPCOND<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Invalid"))+
  stat_summary(aes(x = Year, group = Year, y =SPCOND), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = SPCOND), fun.y = "median",
               size = 3 ,color = "blue", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = SPCOND),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "blue", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Conductivity ($\\mu S$ cm^{-1}$})"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('Conductivity')
InvalidSO4<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Invalid"))+
  stat_summary(aes(x = Year, group = Year, y =SO4), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = SO4), fun.y = "median",
               size = 3 ,color = "Yellow", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = SO4),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Yellow", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{SO4 Concentration ($\\mu eq$L^{-1}$})"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('SO4')
InvalidNO3<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Invalid"))+
  stat_summary(aes(x = Year, group = Year, y = NO3), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = NO3), fun.y = "median",
               size = 3 ,color = "Purple", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NO3),fun ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Purple", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{NO3 Concentration ($\\mu eq $L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('NO3')

InvalidNH4<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Invalid"))+
  stat_summary(aes(x = Year, group = Year, y = NH4), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = NH4), fun.y = "median",
               size = 3 ,color = "Orange", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NH4),fun ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Orange", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{NH4 Concentration ($\\mu eq $L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('NH4')

InvalidTOC<-ggplot(data = subset(AllCloudData, Year < 2018 & Class == "Invalid"))+
  stat_summary(aes(x = Year, group = Year, y = TOC), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = TOC), fun.y = "median",
               size = 3 ,color = "Forest Green", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = TOC),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Forest Green", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{TOC Concentration ($\\mu molC $L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2018,2))+
  ggtitle('TOC')

grid.arrange(InvalidpH,InvalidSPCOND,InvalidSO4,InvalidNH4,InvalidNO3,InvalidTOC)



#### Complete Dataset ####

AllpH<-ggplot(data = subset(AllCloudData, Year < 2018))+
  stat_summary(aes(x = Year, group = Year, y =LABPH), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = LABPH), fun.y = "median",
               size = 3 ,color = "red", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = LABPH),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "red", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{pH}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('pH')

AllSPCOND<-ggplot(data = subset(AllCloudData, Year < 2018))+
  stat_summary(aes(x = Year, group = Year, y =SPCOND), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = SPCOND), fun.y = "median",
               size = 3 ,color = "blue", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = SPCOND),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "blue", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Conductivity ($\\mu S$ cm^{-1}$})"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('Conductivity')
AllSO4<-ggplot(data = subset(AllCloudData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y =K), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = K), fun.y = "median",
               size = 3 ,color = "Yellow", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = K),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Yellow", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{SO4 Concentration ($\\mu eq$L^{-1}$})"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('SO4')
AllNO3<-ggplot(data = subset(AllCloudData, Year < 2018))+
  stat_summary(aes(x = Year, group = Year, y = NO3), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = NO3), fun.y = "median",
               size = 3 ,color = "Purple", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NO3),fun ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Purple", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{NO3 Concentration ($\\mu eq $L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('NO3')

AllNH4<-ggplot(data = subset(AllCloudData, Year < 2018))+
  stat_summary(aes(x = Year, group = Year, y = NH4), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = NH4), fun.y = "median",
               size = 3 ,color = "Orange", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NH4),fun ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Orange", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{NH4 Concentration ($\\mu eq $L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('NH4')

AllTOC<-ggplot(data = subset(AllCloudData, Year < 2018))+
  stat_summary(aes(x = Year, group = Year, y = TOCCorrected), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = TOCCorrected), fun.y = "median",
               size = 3 ,color = "Forest Green", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = TOCCorrected),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Forest Green", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{TOC Concentration ($\\mu molC $L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2018,2))+
  ggtitle('TOC')

grid.arrange(InvalidpH,InvalidSPCOND,InvalidSO4,InvalidNH4,InvalidNO3,InvalidTOC)


ggplot(data = subset(AllCloudData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y = Sodium), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = Sodium), fun.y = "median",
               size = 3 ,color = "Blue", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = Sodium),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Blue", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Na Concentration ($\\mu molC $L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2018,2))+
  ggtitle('Na')


#### TheilSen Function and Results ####


TheilSenFunction<-function(df, analytes ,fun){
  df<-df%>%
    mutate(Hplus = Hplus*1e6)%>%
    group_by(Year)%>%
    summarise(across(where(is.numeric), fun, na.rm = TRUE))
 df<-as.data.frame(df)
 for (i in analytes){
 print(i)
   print(theilsen(df[,i]~Year, data = df))
   Ken<-MannKendall(df[,i])
   print("Mann Kendall Result")
   print(paste("tau:", Ken$tau,"p-value:", Ken$sl))
 }
}

ValidSlope<-TheilSenFunction(subset(AllCloudData, Class == "Valid"), analytes = c("LABPH", "SPCOND", "SO4", "NO3", "NH4", "TOC", "CA", "MG", "K", 
                                                                              "Sodium", "CL", "Ratio"), fun = median)### Loops through all the analytes in a given list.
InvalidSlope<-TheilSenFunction(subset(AllCloudData, Class == "Invalid"), analytes = c("LABPH", "SPCOND", "SO4", "NO3", "NH4", "TOC", "CA", "MG", "K", 
                                                                                    "Sodium", "CL", "Ratio"), fun = median)
AllSlope<-TheilSenFunction(subset(AllCloudData, Year >2008), analytes = c("LABPH", "SPCOND", "SO4", "NO3", "NH4", "TOC", "CA", "MG", "K", 
                                                                                  "Sodium", "CL", "Ratio", "SurplusNH4"), fun = median)

CloudLoadingSlope<-TheilSenFunction(subset(LWCLoadingData, Year > 2008), analytes =c("LWCCalc","SO4Mass","NO3Mass", "NH4Mass", "TOCMass", "CAMass", "MGMass", "KMass", "NaMass", "CLMass"), fun = median)


#### Valid vs Invalid Data Statistical Tests ####

AllCloudData%>%
  dplyr::select(c(Class,LABPH:TOC))%>%
  pivot_longer(cols = !Class, names_to = "Species", values_to = "Conc")%>%
  group_by(Species)%>%
  kruskal_test(Conc~Class)


DunnTestbyClass<-AllCloudData%>%
  dplyr::select(c(Class,LABPH:TOC))%>%
  dplyr::rename(pH = LABPH, Conductivity = SPCOND, Na = Sodium, Ca = CA, Mg = MG, Cl = CL)%>%
  pivot_longer(cols = !Class, names_to = "Species", values_to = "Conc")%>%
  group_by(Species)%>%
  dunn_test(Conc~Class,p.adjust.method = "holm", detailed = TRUE)%>%
  dplyr::select(c("Species", "estimate", "n1", "n2", "p.adj", "p.adj.signif"))

MedianValid<-AllCloudData%>%
  dplyr::select(c(Class,LABPH:TOC))%>%
  dplyr::rename(pH = LABPH, Conductivity = SPCOND, Na = Sodium, Ca = CA, Mg = MG, Cl = CL)%>%
  pivot_longer(cols = !Class, names_to = "Species", values_to = "Conc")%>%
  filter(Class == "Valid")%>%
  group_by(Species)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  dplyr::rename(`Valid Median` = Conc)
MedianInvalid<-AllCloudData%>%
  dplyr::select(c(Class,LABPH:TOC))%>%
  dplyr::rename(pH = LABPH, Conductivity = SPCOND, Na = Sodium, Ca = CA, Mg = MG, Cl = CL)%>%
  pivot_longer(cols = !Class, names_to = "Species", values_to = "Conc")%>%
  filter(Class == "Invalid")%>%
  group_by(Species)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  dplyr::rename(`Invalid Median` = Conc)

DunnTestbyClass<-inner_join(DunnTestbyClass, MedianValid, by = "Species")
DunnTestbyClass<-inner_join(DunnTestbyClass, MedianInvalid, by = "Species")
  
DunnTestbyClass%>%
  dplyr::select(c("Species","n1", "n2", "Invalid Median", "Valid Median", "p.adj", "p.adj.signif"))%>%
  arrange(desc(`Invalid Median`) )%>%
  dplyr::rename(`Valid n` = n2, `Invalid n` = n1)%>%
  mutate(Difference = `Invalid Median` - `Valid Median`)%>%
  dplyr::select(c(Species, `Invalid Median`,`Valid Median`, Difference, p.adj, `Invalid n`, `Valid n`))%>%
  mutate(across(where(is.numeric), round, digits = 3))%>%
  gt()
  



#### Regime Statistical Tests ####
AllCloudData%>%
  mutate(`Ion Balance` = Cations - Anions)%>%
  dplyr::select(c(Regime,LABPH:TOC, `Ion Balance`))%>%
  pivot_longer(cols = !Regime, names_to = "Species", values_to = "Conc")%>%
  group_by(Species)%>%
  kruskal_test(Conc~Regime)

DunnTestbyRegime<-AllCloudData%>%
  mutate(`Ion Balance` = Cations - Anions)%>%
  dplyr::select(c(Regime,LABPH:TOC, `Ion Balance`))%>%
  dplyr::rename(pH = LABPH, Conductivity = SPCOND, Na = Sodium, Ca = CA, Mg = MG, Cl = CL)%>%
  pivot_longer(cols = !Regime, names_to = "Species", values_to = "Conc")%>%
  group_by(Species)%>%
  dunn_test(Conc~Regime,p.adjust.method = "holm", detailed = TRUE)%>%
  dplyr::select(c("Species", "estimate", "n1", "n2", "p.adj", "p.adj.signif"))


MedianOld<-AllCloudData%>%
  mutate(`Ion Balance` = Cations - Anions)%>%
  dplyr::select(c(Regime,LABPH:TOC, `Ion Balance`))%>%
  dplyr::rename(pH = LABPH, Conductivity = SPCOND, Na = Sodium, Ca = CA, Mg = MG, Cl = CL)%>%
  pivot_longer(cols = !Regime, names_to = "Species", values_to = "Conc")%>%
  filter(Regime == "Linear")%>%
  group_by(Species)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  dplyr::rename(`Old Regime Median` = Conc)
MedianNew<-AllCloudData%>%
  mutate(`Ion Balance` = Cations - Anions)%>%
  dplyr::select(c(Regime,LABPH:TOC, `Ion Balance`))%>%
  dplyr::rename(pH = LABPH, Conductivity = SPCOND, Na = Sodium, Ca = CA, Mg = MG, Cl = CL)%>%
  pivot_longer(cols = !Regime, names_to = "Species", values_to = "Conc")%>%
  filter(Regime == "Non-Linear")%>%
  group_by(Species)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  dplyr::rename(`New Regime Median` = Conc)


DunnTestbyRegime<-inner_join(DunnTestbyRegime, MedianOld, by = "Species")
DunnTestbyRegime<-inner_join(DunnTestbyRegime, MedianNew, by = "Species")

DunnTestbyRegime%>%
  dplyr::select(c("Species","n1", "n2", "New Regime Median", "Old Regime Median", "p.adj", "p.adj.signif"))%>%
  arrange(desc(`New Regime Median`) )%>%
  dplyr::rename(`Old Regime n` = n1, `New Regime n` = n2)%>%
  mutate(Difference = `New Regime Median` - `Old Regime Median`)%>%
  dplyr::select(c(Species, `New Regime Median`,`Old Regime Median`, Difference, p.adj, `New Regime n`, `Old Regime n`))%>%
  mutate(across(where(is.numeric), round, digits = 3))%>%
  gt()
  
#### Surplus NH4 and Cation/Anion Ratio #### 

SurplusNH4<-ggplot(data = subset(AllCloudData, Year < 2018))+
#  stat_summary(aes(x = Year, group = Year, y =LABPH), fun.data = median_IQR,
              # size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = 1e6*10^(-LABPH)), fun.y = "median",
               size = 3 ,color = "red", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = 1e6*10^(-LABPH)),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "red", shape = 21)+
  stat_summary(aes(x = Year, y = NH4-SO4-NO3), fun.y = "median",
               size = 3 ,color = "blue", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NH4-SO4-NO3),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "blue", shape = 21)+
  geom_hline(yintercept = 0, size =3)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{H+ or Suplus NH4 Concentrations ($\\mu eq $L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('H+ and Surplus NH4')

CationAnion<-ggplot(data = subset(AllCloudData, Year < 2018))+
  #  stat_summary(aes(x = Year, group = Year, y =LABPH), fun.data = median_IQR,
  # size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = Cations/Anions), fun.y = "median",
               size = 3 ,color = "grey", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = Cations/Anions),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "grey", shape = 21)+
  #geom_hline(yintercept = 0, size =3)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Cations/Anions}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('Cation/Anion Ratio')

grid.arrange(SurplusNH4,CationAnion)

  
#### Cloud Water Loading Plots ####


LoadingLWC<-ggplot(data = subset(LWCLoadingData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y = LWCCalc), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = LWCCalc), fun.y = "median",
               size = 3 ,color = "Red", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = LWCCalc),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Red", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 14, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Liquid Water Content \\ (g $m^{-3}$)}"))+scale_x_continuous(breaks = seq(2010,2018,2))+
  ggtitle('Liquid Water Content')

LoadingSO4<-ggplot(data = subset(LWCLoadingData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y =SO4Mass), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = SO4Mass), fun.y = "median",
               size = 3 ,color = "Yellow", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = SO4Mass),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Yellow", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{SO4 CWL \\ (neq $m^{-3}$})"))+scale_x_continuous(breaks = seq(2010,2018,2))+
  ggtitle('SO4')
LoadingNO3<-ggplot(data = subset(LWCLoadingData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y = NO3Mass), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = NO3Mass), fun.y = "median",
               size = 3 ,color = "Purple", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NO3Mass),fun ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Purple", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{NO3 CWL \\ (neq $m^{-3}$)}"))+scale_x_continuous(breaks = seq(2010,2018,2))+
  ggtitle('NO3')

LoadingNH4<-ggplot(data = subset(LWCLoadingData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y = NH4Mass), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = NH4Mass), fun.y = "median",
               size = 3 ,color = "Orange", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NH4Mass),fun ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Orange", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{NH4 CWL \\ (neq $m^{-3}$})"))+scale_x_continuous(breaks = seq(2010,2018,2))+
  ggtitle('NH4')

LoadingTOC<-ggplot(data = subset(LWCLoadingData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y = TOCMass), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = TOCMass), fun.y = "median",
               size = 3 ,color = "Forest Green", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = TOCMass),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Forest Green", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{TOC CWL \\ (nmolC $m^{-3}$)}"))+scale_x_continuous(breaks = seq(2010,2018,2))+
  ggtitle('TOC')

LoadingCa<-ggplot(data = subset(LWCLoadingData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y = CAMass), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = CAMass), fun.y = "median",
               size = 3 ,color = "Blue", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = CAMass),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Blue", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Ca CWL \\ (neq $m^{-3}$)}"))+scale_x_continuous(breaks = seq(2010,2018,2))+
  ggtitle('CA')



LoadingCL<-ggplot(data = subset(LWCLoadingData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y = CLMass), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = CLMass), fun.y = "median",
               size = 3 ,color = "Blue", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = CLMass),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Blue", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Ca CWL \\ (neq $m^{-3}$)}"))+scale_x_continuous(breaks = seq(2010,2018,2))+
  ggtitle('Cl')


LoadingNa<-ggplot(data = subset(LWCLoadingData, Year < 2021))+
  stat_summary(aes(x = Year, group = Year, y = NaMass), fun.data = median_IQR,
               size = 1 ,color = "black")+
  stat_summary(aes(x = Year, y = NaMass), fun.y = "median",
               size = 3 ,color = "Blue", geom = "line")+
  stat_summary(aes(x = Year,group = Year, y = NaMass),fun.y ="median" ,geom = "point", 
               size = 5, color = "black", fill = "Blue", shape = 21)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Na CWL \\ (neq $m^{-3}$)}"))+scale_x_continuous(breaks = seq(2010,2018,2))+
  ggtitle('Na')

grid.arrange(LoadingCL, LoadingNa, ncol = 2)
grid.arrange(LoadingLWC,LoadingSO4,LoadingNO3,LoadingNH4,LoadingTOC, LoadingCa)


#### 


ggplot(subset(LWCLoadingData, Year < 2018))+
  geom_point(aes(x = NH4Mass, y = TOCMass, color = "Mass Loading"))+
  geom_point(aes(x= NH4, y = WSOC, color = "Cloud Concentration"))+
  stat_regline_equation(aes(x = NH4Mass, y = TOCMass, color = "Mass Loading",
                            label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")))
  scale_color_viridis_c()+xlim(c(0,1))


  
  
Loading2016<-LWCLoadingData%>%
  filter(Year == 2016)%>%
  ggplot()+
  geom_point(aes(x= date, y = TOCMass, color = "TOC"))+
  geom_point(aes(x= date, y = KMass*80, color = "K"))

Loading2015<-LWCLoadingData%>%
  filter(Year == 2017)%>%
  ggplot()+
  geom_point(aes(x= date, y = TOCMass, color = "TOC"))+
  geom_point(aes(x= date, y = KMass*80, color = "K"))

LWCLoadingData%>%
#  filter(K/CA > 0)%>%
  filter(Year == 2018)%>%
  ggplot()+
  geom_point(aes(x= date, y = K/WSOC, color = "TOC"))
  #geom_point(aes(x= date, y = (K/CA)*1000, color = "K"))




