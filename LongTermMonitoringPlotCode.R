#install.packages("ggedit")
library(extrafont)
library(gridExtra)
library(tidymodels)
library(tidyverse)
library(latex2exp)
library(gridExtra)
library(ggpubr)
library(gtsummary)
library(deming)
library(Kendall)
library(scales)
library(patchwork)
library(zyp)
library(car)
library(vip)

##
show_col(hue_pal()(6))
colsgg<-(hue_pal()(6))

### Data Clean Up ####
AllCloudData<-read.csv("Alldata.csv", header = TRUE)

PaperQuality<-function(...){
  theme_bw()%+replace%
    theme(axis.title.x = element_text(size = 17), axis.title.y = element_text(size = 17, angle = 90),
          axis.text = element_text(size = 17), title = element_text(size = 17, face = "bold"),
          legend.text = element_text(size = 15), legend.title = element_blank(),
          legend.position = "bottom")
}

AllCloudData<-AllCloudData%>%
  dplyr::rename(TOC = WSOC)%>%
  filter(Year > 1993 & Year<2021)%>%
  mutate(date = as.POSIXct(date, format = "%m/%d/%Y %H:%M"))%>%
#  mutate(TOCCorrected = case_when(Year == 2018 | Year == 2019 ~ TOC*(1/.8486), 
 #                                 TRUE~TOC))%>%
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
  mutate(HCO3Total = (HCO3Fraction*(MG+CA)/2) + 1e6*HCO3Gas,
         pHBin = cut(LABPH, 6))%>%
  mutate(IonBalance = Cations-Anions,
         IonBalanceHCO3 = Cations - Anions -HCO3Total)

### Create Median Trends Separated by Valid and All Data ####

ValidMedians<-AllCloudData%>%
  filter(Class == "Valid")%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))

InvalidMedians<-AllCloudData%>%
  filter(Class == "Invalid")%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))

AllMedians<-AllCloudData%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))

AllMedians<-cbind(AllMedians, ValidMedians)
names(AllMedians)<-make.unique(names(AllMedians))
names(AllMedians)<-gsub(".1", "Valid", x = names(AllMedians))

AllMedians<-cbind(AllMedians, InvalidMedians)
names(AllMedians)<-make.unique(names(AllMedians))
names(AllMedians)<-gsub(".1", "Invalid", x = names(AllMedians))


## Run the Sen Slope Functions for Valid and Total Data Sets and Create Datatables####

sen <- function(..., weights = NULL) {
  mblm::mblm(...)
} 

TheilSenFunction<-function(df, analytes ,fun){
  df<-df%>%
    mutate(Hplus = Hplus*1e6)%>%
    group_by(Year)%>%
    summarise(across(where(is.numeric), fun, na.rm = TRUE))
  df<-as.data.frame(df)
  SenSlope<-list()
  SenInter<-list()
  #Coeff<-list()
  Ken<-list()
  for (i in analytes){
    #print(i)
    # print(theilsen(df[,i]~Year, data = df))
    #SenSlope[[i]]<-theilsen(df[,i]~Year, data = df)
    Slope<-theilsen(df[,i]~Year, data = df)
    SenSlope[[i]]<-Slope$coefficients[[2]]
    SenInter[[i]]<-Slope$coefficients[[1]]
    pvalue<-MannKendall(df[,i])
    Ken[[i]]<-pvalue$sl
    #print(Sen$coefficients[[2]])
    #print("Mann Kendall Result")
    #print(paste("tau:", Ken$tau,"p-value:", Ken$sl, sep = " "))
  }
  SenDF<-tibble(Analyte = unlist(analytes), Slope = unlist(SenSlope), Intercept =unlist(SenInter) , `P-Value` = unlist(Ken))
  SenDF<-SenDF%>%
    mutate(`P-Value` = scientific(`P-Value`, digits = 3),
           Slope = round(Slope, digits = 2),
           Intercept = round(Intercept, digits = 3))
  return(SenDF)
  }


ValidSlope<-TheilSenFunction(df = subset(AllCloudData, Class == "Valid"), analytes = list("LABPH", "SPCOND", "SO4", "NO3", "NH4", "TOC", "CA", "MG", "K", 
                                                                                     "Sodium", "CL"), fun = median)
AllSlope<-TheilSenFunction(AllCloudData, analytes = list("LABPH", "SPCOND", "SO4", "NO3", "NH4", "TOC", "CA", "MG", "K", 
                                                                                   "Sodium", "CL", "Ratio"), fun = median)
ValidSlope$Analyte<-c("pH", "Conductivity", "SO4", "NO3", "NH4", "TOC", "Ca", "Mg", "K", "Na", "Cl")
AllSlope$Analyte<-c("pH", "Conductivity", "SO4", "NO3", "NH4", "TOC", "Ca", "Mg", "K", "Na", "Cl", "Ratio")

ValidSlopeNoInter<-ValidSlope%>%
  select(-Intercept)
ValidSlopeTable<-ggtexttable(ValidSlopeNoInter, theme = ttheme("light", base_size = 16, colnames.style = colnames_style(face = "italic", size = 16)), rows = NULL)%>%
  tab_add_title(text = paste("Theil-Sen Slope and",
                              "Mann Kendall P-Value", sep = "\n"), size = 17, face = "bold",padding = unit(1, "line"),
                just = "left")

AllSlopeNoInter<-AllSlope%>%
  select(-Intercept)%>%
  filter(Analyte != "Ratio")
AllSlopeTable<-ggtexttable(AllSlopeNoInter, theme = ttheme("light", base_size = 16, colnames.style = colnames_style(face = "italic", size = 16)), rows = NULL)%>%
  tab_add_title(text = paste("Theil-Sen Slope and",
                             "Mann Kendall P-Value", sep = "\n"), size = 17, face = "bold",padding = unit(1, "line"),
                just = "left")



## Figure 2 Valid Plots ####



ValidpHPlot<-ggplot(ValidMedians)+
  geom_line(aes(x = Year, y = LABPH), color = "purple",size = 4)+
  geom_point(aes(x = Year, y = LABPH, fill = "pH"), size = 4, shape = 21)+
  geom_smooth(aes(x= Year, y = LABPH), color = "purple", se = FALSE, method = sen, size = 2)+
  scale_fill_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"), axis.ticks = element_line(color = "forest green"), 
                       axis.title.y.left = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
  ylab("pH")+scale_x_continuous(breaks = seq(1994, 2018, by = 4))+
  theme(legend.position = "none")+
  annotate(geom = 'text',x = 2008, y = 5, size = 6,
           label = paste("y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="pH"][[1]],digits = 3),"x +", signif(ValidSlope$Intercept[ValidSlope$Analyte=="pH"][[1]],digits = 3)), color = "purple")+
  ggtitle("pH")


ValidConductivityPlot<-ggplot(ValidMedians)+
  geom_line(aes(x = Year, y = SPCOND), color = "black",size = 4)+
  geom_point(aes(x = Year, y = SPCOND, fill = "Conductivity"), size = 4, shape = 21)+
  geom_smooth(aes(x= Year, y = SPCOND), color = "black", se = FALSE, method = sen, size = 2)+
  scale_fill_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  PaperQuality()+theme(legend.position = "none")+
  ylab(TeX("\\textbf{Conductivity ($\\mu S $ cm^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2017, by = 4))+
  annotate(geom = 'text',x = 2010, y = 90, size = 6,
           label = paste("y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="Conductivity"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="Conductivity"][[1]], digits = 3)), color = "black")+
  ggtitle("Conductivity")

ValidConcPlots<-ggplot(ValidMedians)+
  geom_line(aes(x = Year, y = SO4), color = "red",size = 4)+
  geom_point(aes(x = Year, y = SO4, fill = "SO4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = SO4, color = "SO4"), method = sen, size =2, se = FALSE)+
  geom_line(aes(x = Year, y = NO3),color = "blue", size = 4)+
  geom_point(aes(x = Year, y = NO3, fill = "NO3"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = NO3, color = "NO3"), method = sen, size = 2, se = FALSE)+
  geom_line(aes(x = Year, y = NH4), color = "orange", size = 4)+
  geom_point(aes(x = Year, y = NH4, fill = "NH4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = NH4, color = "NH4"), method = sen, size = 2, se = FALSE)+
  geom_line(aes(x = Year, y = TOC/4),color = "forest green",size = 4)+
  geom_point(aes(x = Year, y = TOC/4, fill = "TOC"), size = 4, shape =21)+
  geom_smooth(aes(x = Year, y = TOC/4, color = "TOC"), method = sen, size = 2, se = FALSE)+
  scale_x_continuous(breaks = seq(1994, 2018, by = 4))+scale_y_continuous(guide = guide_axis(check.overlap = TRUE),sec.axis = sec_axis(trans = ~.x*4, name = TeX("\\textbf{TOC Concentration ($\\mu molC $ L^{-1}$)}")),
                                                                          name =TeX("\\textbf{$Ion Concentration$ ($\\mu eq $ L^{-1}$)}"))+
  scale_fill_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  scale_color_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"), guide = FALSE)+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"),
                       axis.ticks = element_line(color = "forest green"), axis.title.y.right = element_text(color = "forest green", size = 15))+
  annotate(geom = 'text',x = 2012, y = 200, size = 6,
           label = paste("y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="SO4"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="SO4"][[1]], digits = 3)), color = "red")+
  annotate(geom = 'text',x = 2012, y = 175, size = 6,
           label = paste("y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="NH4"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="NH4"][[1]], digits = 3)), color = "orange")+
  annotate(geom = 'text',x = 2012, y = 150, size = 6,
           label = paste("y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="NO3"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="NO3"][[1]], digits = 3)), color = "blue")+
  annotate(geom = 'text',x = 2012, y = 125, size = 6,
           label = paste("y=",signif(ValidSlope$Slope[ValidSlope$Analyte=="TOC"][[1]],digits = 3),"x +", 
                         signif(ValidSlope$Intercept[ValidSlope$Analyte=="TOC"][[1]], digits = 3)), color = "forest green")+
  ggtitle("Major Analytes")

 
ValidPlots<-ggarrange(ValidConcPlots,ValidpHPlot, ValidConductivityPlot, ValidSlopeTable, nrow = 2, ncol = 2)


ggsave(ValidPlots, filename = "ValidPlots.png",width = 14, height = 10)




## Figure 3 Percent Valid Plots ####

# Create Precent Valid and Invalid
Percent<-AllCloudData%>%
  select(-TOC)%>% ##Include whole dataset, as TOC start at 2009
  group_by(Year,Class)%>%
  summarise(n = n())%>%
  group_by(Year)%>%
  na.omit()%>%
  summarise(Percent = 100*n/(sum(n)))
Percent<-Percent[seq(1, nrow(Percent), by =2),] ##Filters out data to only look at percent invalid


PercentPlot<-ggplot(Percent)+
  geom_line(aes(Year, y = Percent), color = 'red', size = 2)+
  geom_point(aes(Year, y = Percent), fill = "red", shape = 21, size = 4)+
  PaperQuality()+scale_y_continuous(TeX('\\textbf{Percent Invalid (%)}'))+scale_x_continuous(breaks = seq(1994, 2017, 4))

ggsave(PercentPlot, filename = "PercentPlot.png", width =5, height = 5)

##Figure 4 Comparing Valid and Invalid Medians####

ValidInvalidpH<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = LABPHInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = LABPHInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = LABPHValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = LABPHValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  #scale_fill_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"), axis.ticks = element_line(color = "forest green"), 
                       axis.title.y.left = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)))+
  ggtitle("pH")+
  ylab("pH")+scale_x_continuous(breaks = seq(1994, 2018, by = 4))

ValidInvalidCond<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = SPCONDInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = SPCONDInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = SPCONDValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = SPCONDValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{Conductivity ($\\mu S $ cm^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2017, by = 4))+
  ggtitle("Conductivity")

ValidInvalidSO4<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = SO4Invalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = SO4Invalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = SO4Valid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = SO4Valid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{SO_4^{-2} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2017, by = 4))+
  ggtitle(TeX("\\textbf{SO_4^{ 2-}}"))

ValidInvalidNH4<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = NH4Invalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = NH4Invalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = NH4Valid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = NH4Valid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{NH_4^{ +} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2017, by = 4))+
  ggtitle(TeX("\\textbf{NH_4^{ +}}"))

ValidInvalidNO3<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = NO3Invalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = NO3Invalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = NO3Valid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = NO3Valid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{NO_3^{ -} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2017, by = 4))+
  ggtitle(TeX("\\textbf{NO_3^{ -}}"))

ValidInvalidTOC<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = TOCInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = TOCInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = TOCValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = TOCValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{TOC ($\\mu mol $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(2008, 2018, by = 2), limits = c(2009,2018))+
  ggtitle(TeX("\\textbf{TOC"))

ValidInvalidCA<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = CAInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = CAInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = CAValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = CAValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{Ca^{2+} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2018, by = 4), limits = c(1994,2018))+
  ggtitle(TeX("\\textbf{Ca^{2+}}"))


ValidInvalidK<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = KInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = KInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = KValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = KValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{K^{+} ($\\mu eq $ L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2018, by = 4), limits = c(1994,2018))+
  ggtitle(TeX("\\textbf{K^{+}}"))



ValidInvalidLWC<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = LWCInvalid), color = '#F8766D',size = 4)+
  geom_point(aes(x = Year, y = LWCInvalid, fill = "Invalid"), color = 'black', size = 4, shape = 21)+
  geom_line(aes(x = Year, y = LWCValid),color ='#00BFC4' ,size = 4)+
  geom_point(aes(x = Year, y = LWCValid, fill = "Valid"),color = 'black', size = 4, shape = 21)+
  PaperQuality()+
  ylab(TeX("\\textbf{LWC (g m^{-3}$)}"))+scale_x_continuous(breaks = seq(1994, 2018, by = 4), limits = c(1994,2018))+
  ggtitle(TeX("\\textbf{LWC}"))




ValidInvalidAll<-ValidInvalidpH+ValidInvalidCond+ValidInvalidSO4+ValidInvalidNH4+ValidInvalidNO3+ValidInvalidTOC+
  ValidInvalidCA+ValidInvalidK+ValidInvalidLWC&theme(legend.position = "bottom")
ValidInvalidAll<-ValidInvalidAll+plot_layout(guides = "collect")
ggsave(filename = "ValidInvalidAll.png", ValidInvalidAll,height = 10,width = 16)

## Figure 5 Complete Dataset Trends####
AllpHPlot<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = LABPH), color = "purple",size = 4)+
  geom_point(aes(x = Year, y = LABPH, fill = "pH"), size = 4, shape = 21)+
  geom_smooth(aes(x= Year, y = LABPH), color = "purple", se = FALSE, method = sen, size = 2)+
  scale_fill_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"), axis.ticks = element_line(color = "forest green"), 
                       axis.title.y.left = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)), legend.position = "none")+
  ylab("pH")+scale_x_continuous(breaks = seq(1994, 2018, by = 4))+
  annotate(geom = 'text',x = 2010, y = 5, size = 6,
           label = paste("y=",signif(AllSlope$Slope[AllSlope$Analyte=="pH"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="pH"][[1]], digits = 3)), color = "purple")+
  ggtitle("pH")



AllConductivityPlot<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = SPCOND), color = "black",size = 4)+
  geom_point(aes(x = Year, y = SPCOND, fill = "Conductivity"), size = 4, shape = 21)+
  geom_smooth(aes(x= Year, y = SPCOND), color = "black", se = FALSE, method = sen, size = 2)+
  scale_fill_manual(values = c("pH" = 'purple', 'Conductivity' = "black","SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  PaperQuality()+
  ylab(TeX("\\textbf{Conductivity ($\\mu S $ cm^{-1}$)}"))+scale_x_continuous(breaks = seq(1994, 2017, by = 4))+
  annotate(geom = 'text',x = 2010, y = 60, size = 6,
           label = paste("y=",signif(AllSlope$Slope[AllSlope$Analyte=="Conductivity"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="Conductivity"][[1]], digits = 3)), color = "black")+
  ggtitle("Conductivity")+theme(legend.position = 'none')

AllConcPlots<-ggplot(AllMedians)+
  geom_line(aes(x = Year, y = SO4), color = "red",size = 4)+
  geom_point(aes(x = Year, y = SO4, fill = "SO4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = SO4, color = "SO4"), method = sen, size =2, se = FALSE)+
  geom_line(aes(x = Year, y = NO3),color = "blue", size = 4)+
  geom_point(aes(x = Year, y = NO3, fill = "NO3"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = NO3, color = "NO3"), method = sen, size = 2, se = FALSE)+
  geom_line(aes(x = Year, y = NH4), color = "orange", size = 4)+
  geom_point(aes(x = Year, y = NH4, fill = "NH4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = NH4, color = "NH4"), method = sen, size = 2, se = FALSE)+
  geom_line(aes(x = Year, y = TOC/4),color = "forest green",size = 4)+
  geom_point(aes(x = Year, y = TOC/4, fill = "TOC"), size = 4, shape =21)+
  geom_smooth(aes(x = Year, y = TOC/4, color = "TOC"), method = sen, size = 2, se = FALSE)+
  scale_x_continuous(breaks = seq(1994, 2018, by = 4))+scale_y_continuous(guide = guide_axis(check.overlap = TRUE),sec.axis = sec_axis(trans = ~.x*4, name = TeX("\\textbf{TOC Concentration ($\\mu molC $ L^{-1}$)}")),
                                                                          name =TeX("\\textbf{$Ion Concentration$ ($\\mu eq $ L^{-1}$)}"))+
  scale_fill_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  scale_color_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"), guide = FALSE)+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"),
                       axis.ticks = element_line(color = "forest green"), axis.title.y.right = element_text(color = "forest green", size = 13.4))+
  annotate(geom = 'text',x = 2010, y = 200, size = 6,
           label = paste("y=",signif(AllSlope$Slope[AllSlope$Analyte=="SO4"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="SO4"][[1]], digits = 3)), color = "red")+
  annotate(geom = 'text',x = 2010, y = 175, size = 6,
           label = paste("y=",signif(AllSlope$Slope[AllSlope$Analyte=="NH4"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="NH4"][[1]], digits = 3)), color = "orange")+
  annotate(geom = 'text',x = 2010, y = 150, size = 6,
           label = paste("y=",signif(AllSlope$Slope[AllSlope$Analyte=="NO3"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="NO3"][[1]], digits = 3)), color = "blue")+
  annotate(geom = 'text',x = 2010, y = 125, size = 6,
           label = paste("y=",signif(AllSlope$Slope[AllSlope$Analyte=="TOC"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="TOC"][[1]], digits = 3)), color = "forest green")+
  ggtitle("Major Analytes")


AllPlots<-ggarrange(AllConcPlots,AllpHPlot, AllConductivityPlot, AllSlopeTable, nrow = 2, ncol = 2)

ggsave(AllPlots, filename = "AllPlots.png",width = 14, height = 10)

## Figure 6 Surplus NH4 and Cation/Anion Ratio #### 
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
  scale_y_continuous(name = TeX("\\textbf{Concentrations ($\\mu eq $L^{-1}$)}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('H+ and Surplus NH4')+
  annotate(x = 2010, y = -75, geom = "text", label = TeX("\\textbf{NH_4^+ - SO_4^{-2} - NO_3^-}"), color = "blue", size = 6)+
  annotate(x = 2010, y = 75, geom = "text", label = TeX("\\textbf{H^+}"), color = 'red', size = 6)

CationAnion<-ggplot(data = subset(AllMedians, Year < 2018))+
  #  stat_summary(aes(x = Year, group = Year, y =LABPH), fun.data = median_IQR,
  # size = 1 ,color = "black")+
  geom_line(aes(x = Year, y = Ratio), 
               size = 3 ,color = "grey", geom = "line")+
  geom_point(aes(x = Year, y = Ratio), 
               size = 5, color = "black", fill = "grey", shape = 21)+
  #geom_hline(yintercept = 0, size =3)+
  theme_bw()+theme(axis.text.y = element_text(size = 18),axis.text.x = element_text(size = 18),axis.title = element_text(size = 16, face = "bold"))+
  theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), title = element_text(face = "bold", size = 18))+
  scale_y_continuous(name = TeX("\\textbf{Cations/Anions}"))+scale_x_continuous(breaks = seq(1994,2017,4))+
  ggtitle('Cation/Anion Ratio')+
  annotate(geom = 'text',x = 1998, y = 1.4, size = 6,
           label = paste("y=",signif(AllSlope$Slope[AllSlope$Analyte=="Ratio"][[1]],digits = 3),"x +", 
                         signif(AllSlope$Intercept[AllSlope$Analyte=="Ratio"][[1]], digits = 3)))+
  geom_smooth(aes(x = Year, y = Ratio),color = "grey", method = sen, size = 2, se = FALSE)
  

NH4CationRatio<-grid.arrange(SurplusNH4,CationAnion, layout_matrix = rbind(c(1,NA,1,NA), c(2,NA,2,NA)))

ggsave(filename = "SurplusNH4_CationAnion.png", NH4CationRatio, width = 8, height = 8)


## Figure 7 HCO3 Estimates affecting Ion Balance ####

HCO3IonBalanceAll<-ggplot(AllCloudData)+
  geom_point(aes(x = Cations, y = Anions, color = "No HCO3"), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions, color = "No HCO3"), method = lm , formula = y~x+0, size =2)+
  geom_point(aes(x = Cations, y = Anions+HCO3Total, color = "HCO3 Included"), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions+HCO3Total, color = "HCO3 Included"), method = lm , formula = y~x+0, size = 2)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed" , size = 2)+
  PaperQuality()+scale_x_continuous(name = TeX("\\textbf{Cations ($\\mu eq $ L^{-1}$)}"))+scale_y_continuous(name = TeX("\\textbf{Anions ($\\mu eq $ L^{-1}$)}"))+
  stat_regline_equation(aes(x = Cations, y = Anions, color = "No HCO3", label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), 
                        label.x = 250, label.y = 4750,formula = y~x+0, size = 6)+
  stat_regline_equation(aes(x= Cations, y = Anions, color = "HCO3 Included",label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                        label.x = 250, label.y = 3750,formula = y~x+0, size = 6)+ggtitle("Ion Balance All Data")

HCO3IonBalanceHighpH<-ggplot(subset(AllCloudData, LABPH > 5))+
  geom_point(aes(x = Cations, y = Anions, color = "No HCO3"), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions, color = "No HCO3"), method = lm , formula = y~x+0, size = 2)+
  geom_point(aes(x = Cations, y = Anions+HCO3Total, color = "HCO3 Included"), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions+HCO3Total, color = "HCO3 Included"), method = lm , formula = y~x+0, size = 2)+
  PaperQuality()+scale_x_continuous(name = TeX("\\textbf{Cations ($\\mu eq $ L^{-1}$)}"))+scale_y_continuous(name = TeX("\\textbf{Anions ($\\mu eq $ L^{-1}$)}"))+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed" , size = 2)+
  stat_regline_equation(aes(x = Cations, y = Anions, color = "No HCO3", label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), 
                        label.x = 250, label.y = 3000,formula = y~x+0, size = 6)+
  stat_regline_equation(aes(x= Cations, y = Anions+HCO3Total, color = "HCO3 Included",label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                        label.x = 250, label.y = 2250,formula = y~x+0, size = 6)+ggtitle("Ion Balance for pH > 5")



IonBalancePlots<-grid.arrange(HCO3IonBalanceAll, HCO3IonBalanceHighpH)

ggsave(filename = "CationAnionBalance.png", IonBalancePlots, width = 8, height = 8)


#Figure 8 Cations vs Anions as Function of pH ####

IonBalancebypH<-ggplot(subset(AllCloudData, !is.na(pHBin)))+
  geom_point(aes(x = Cations, y = Anions, color = LABPH), size = 2)+
  geom_smooth(aes(x= Cations,y = Anions), method = lm , formula = y~x+0, size =2)+
  PaperQuality()+scale_x_continuous(name = TeX("\\textbf{Cations ($\\mu eq $ L^{-1}$)}"))+scale_y_continuous(name = TeX("\\textbf{Anions ($\\mu eq $ L^{-1}$)}"))+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed" , size = 2)+
  facet_wrap(~pHBin)+scale_color_gradientn(colors = rainbow(6))+
  theme(legend.position = "right", legend.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"))+
  stat_regline_equation(aes(x = Cations, y = Anions, label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")), 
                      label.x = 250, label.y = 4500,formula = y~x+0, size = 5)
ggsave(plot = IonBalancebypH, file = "IonBalancebypH.png", height = 8, width = 12)


## Figure 9 TOC vs Ion Balance, All Data ####

IonBalanceSenAll<-zyp.sen(TOC~IonBalance, 
                          dataframe = na.omit(AllCloudData))


TOCBalanceAll<-ggplot(na.omit(AllCloudData),aes(x= Cations-Anions, y = TOC))+
  geom_point(aes(x= (Cations-Anions), y = TOC), size = 3, color = 'forest green')+
  #  geom_smooth(method = 'rlm', aes(color = "RLM"))+
  # geom_smooth(method = "rlm")+
  geom_smooth(method = lm, size = 3)+
 # geom_abline(aes(intercept = IonBalanceSenNonLinear$coefficients[[1]], 
    #              slope = IonBalanceSenNonLinear$coefficients[[2]]), 
   #           size = 2, color = "red", 
  #            data = AllCloudData%>%filter(Regime == "Non-Linear"))+
  scale_color_gradientn(colors = rainbow(6))+
  labs(color = "pH")+ylim(c(0,3500))+ylab(TeX("\\textbf{TOC ($\\mu molC $ L^{-1}$)}"))+
  xlab(TeX("\\textbf{Cations - Anions ($\\mu eq $ L^{-1}$)}"))+
 # geom_text(aes(x = 0, y = 3000), size = 6,
    #        data = AllCloudData%>%filter(Regime == "Non-Linear"),
     #       label = paste("y=",signif(IonBalanceSenAll$coefficients[[1]],digits = 3),"+", 
             #             signif(IonBalanceSenAll$coefficients[[2]],digits = 3),"x"), color = "red")+
  # stat_regline_equation(method = "rlm")
  stat_regline_equation(size =6, label.x = 0, label.y = 3000,
                        aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")))+
  PaperQuality()+theme(legend.position = 'right', 
                       legend.title = element_text(size =14),
                       strip.text.y = element_text(face = "bold", size = 14))

ggsave(plot = TOCBalanceAll, filename = "TOCBalanceAll.png", width = 8, height = 8)


## Figure 10 TOC vs LWC ####
#There is a seperate script for calculating the Average LWC Data

LWCLoadingData<-read.csv("CloudLWCLoading.csv", header = TRUE)%>%
  mutate(TOC = WSOC)%>%
  mutate()

TOCvsLWC<-ggplot(LWCLoadingData)+
  geom_point(aes(x = LWCCalc, y = TOC), color = "Forest Green",size = 3)+
  PaperQuality()+scale_y_continuous(name = TeX('\\textbf{TOC ($\\mu molC$ L^{-1})}'))+
  scale_x_continuous(name = TeX('\\textbf{LWC (g m^{-3})}'))



ggsave(plot = TOCvsLWC, filename = "TOCLWC.png", width = 8, height = 8)


## Figure 11 Cloud Water Loadings ####

CWLSlope<-TheilSenFunction(LWCLoadingData, analytes = list("LWCCalc","SO4Mass", "NO3Mass", "NH4Mass", "TOCMass", "CAMass", "MGMass", "KMass", 
                                                         "NaMass", "CLMass"), fun = median)
CWLSlope$Analyte<-c("LWC","SO4","NO3", "NH4", "TOC", "Ca", "Mg", "K", "Na", "Cl")
CWLNoInter<-CWLSlope%>%
  select(-Intercept)%>%
  filter(Analyte != "Ratio")
CWLTable<-ggtexttable(CWLNoInter, theme = ttheme("light", base_size = 16, colnames.style = colnames_style(face = "italic", size = 16)), rows = NULL)%>%
  tab_add_title(text = paste("Theil-Sen Slope and",
                             "Mann Kendall P-Value", sep = "\n"), size = 17, face = "bold",padding = unit(1, "line"),
                just = "left")

CWLMedians<-LWCLoadingData%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))


CWLMassPlots<-ggplot(CWLMedians)+
  geom_line(aes(x = Year, y = SO4Mass), color = "red",size = 4)+
  geom_point(aes(x = Year, y = SO4Mass, fill = "SO4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = SO4Mass, color = "SO4"), method = sen, size =2, se = FALSE)+
  geom_line(aes(x = Year, y = NO3Mass),color = "blue", size = 4)+
  geom_point(aes(x = Year, y = NO3Mass, fill = "NO3"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = NO3Mass, color = "NO3"), method = sen, size = 2, se = FALSE)+
  geom_line(aes(x = Year, y = NH4Mass), color = "orange", size = 4)+
  geom_point(aes(x = Year, y = NH4Mass, fill = "NH4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = NH4Mass, color = "NH4"), method = sen, size = 2, se = FALSE)+
  geom_line(aes(x = Year, y = TOCMass/4),color = "forest green",size = 4)+
  geom_point(aes(x = Year, y = TOCMass/4, fill = "TOC"), size = 4, shape =21)+
  geom_smooth(aes(x = Year, y = TOCMass/4, color = "TOC"), method = sen, size = 2, se = FALSE)+
  scale_x_continuous(breaks = seq(1994, 2018, by = 4))+scale_y_continuous(guide = guide_axis(check.overlap = TRUE),sec.axis = sec_axis(trans = ~.x*4, name = TeX("\\textbf{TOC Concentration (nmolC m^{-3})}")),
                                                                          name =TeX("\\textbf{Ion CWL (neq  m^{-3}$)}"))+
  scale_fill_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  scale_color_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"), guide = FALSE)+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"),
                       axis.ticks = element_line(color = "forest green"), axis.title.y.right = element_text(color = "forest green", size = 15))+
  annotate(geom = 'text',x = 2012, y = 80, size = 6,
           label = paste("y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="SO4"][[1]],digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="SO4"][[1]], digits = 3)), color = "red")+
  annotate(geom = 'text',x = 2012, y = 72, size = 6,
           label = paste("y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="NH4"][[1]],digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="NH4"][[1]], digits = 3)), color = "orange")+
  annotate(geom = 'text',x = 2012, y = 64, size = 6,
           label = paste("y=",signif(ValidSlope$Slope[CWLSlope$Analyte=="NO3"][[1]],digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="NO3"][[1]], digits = 3)), color = "blue")+
  annotate(geom = 'text',x = 2012, y = 56, size = 6,
           label = paste("y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="TOC"][[1]],digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="TOC"][[1]], digits = 3)), color = "forest green")+
  ggtitle("Major Analytes")



CWLMassPlots<-ggplot(CWLMedians)+
  geom_line(aes(x = Year, y = SO4Mass), color = "red",size = 4)+
  geom_point(aes(x = Year, y = SO4Mass, fill = "SO4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = SO4Mass, color = "SO4"), method = sen, size =2, se = FALSE)+
  geom_line(aes(x = Year, y = NO3Mass),color = "blue", size = 4)+
  geom_point(aes(x = Year, y = NO3Mass, fill = "NO3"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = NO3Mass, color = "NO3"), method = sen, size = 2, se = FALSE)+
  geom_line(aes(x = Year, y = NH4Mass), color = "orange", size = 4)+
  geom_point(aes(x = Year, y = NH4Mass, fill = "NH4"), size = 4, shape = 21)+
  geom_smooth(aes(x = Year, y = NH4Mass, color = "NH4"), method = sen, size = 2, se = FALSE)+
  geom_line(aes(x = Year, y = TOCMass/4),color = "forest green",size = 4)+
  geom_point(aes(x = Year, y = TOCMass/4, fill = "TOC"), size = 4, shape =21)+
  geom_smooth(aes(x = Year, y = TOCMass/4, color = "TOC"), method = sen, size = 2, se = FALSE)+
  scale_x_continuous(breaks = seq(1994, 2018, by = 4))+scale_y_continuous(guide = guide_axis(check.overlap = TRUE),sec.axis = sec_axis(trans = ~.x*4, name = TeX("\\textbf{TOC Concentration (nmolC m^{-3})}")),
                                                                          name =TeX("\\textbf{Ion CWL (neq  m^{-3}$)}"))+
  scale_fill_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"))+
  scale_color_manual(values = c("SO4" = "red", "NO3" = "blue", "NH4" = 'orange', "TOC" = "forest green"), guide = FALSE)+
  PaperQuality()+theme(axis.text.y.right = element_text(color = "forest green"),
                       axis.ticks = element_line(color = "forest green"), axis.title.y.right = element_text(color = "forest green", size = 15))+
  annotate(geom = 'text',x = 2012, y = 80, size = 6,
           label = paste("y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="SO4"][[1]],digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="SO4"][[1]], digits = 3)), color = "red")+
  annotate(geom = 'text',x = 2012, y = 72, size = 6,
           label = paste("y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="NH4"][[1]],digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="NH4"][[1]], digits = 3)), color = "orange")+
  annotate(geom = 'text',x = 2012, y = 64, size = 6,
           label = paste("y=",signif(ValidSlope$Slope[CWLSlope$Analyte=="NO3"][[1]],digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="NO3"][[1]], digits = 3)), color = "blue")+
  annotate(geom = 'text',x = 2012, y = 56, size = 6,
           label = paste("y=",signif(CWLSlope$Slope[CWLSlope$Analyte=="TOC"][[1]],digits = 3),"x +", 
                         signif(CWLSlope$Intercept[CWLSlope$Analyte=="TOC"][[1]], digits = 3)), color = "forest green")+
  ggtitle("Major Analytes")





CWLPlots<-ggarrange(CWLMassPlots, CWLTable, nrow = 1, ncol = 2)


ggsave(plot = CWLPlots, filename = "CWLPlots.png",width = 12, height = 6)

## Figure 12 Conductivity vs pH ####

HCondPlot<-ggplot(AllCloudData)+
  geom_point(aes(x= LABPH, y = SPCOND, color = HCond/SPCOND),
             size = 3)+
  scale_color_gradientn(colors = rainbow(6), limits = c(0,1))+
  scale_x_continuous(name = "pH")+
  scale_y_log10(limits =c(1,1000),
                name = TeX('\\textbf{Conductivity ($\\mu S$ cm^{-1})}'))+
  labs(color = TeX("\\textbf{HCond Fraction}"))+
  PaperQuality()+theme(legend.position = "right",
                       legend.title = element_text(size = 14))
  

ggsave(plot = HCondPlot, filename = "HCondPlot.png", 
       height = 6, width = 10)


## Figure 13 Percent in the New Regime #### 

PercentRegime<-AllCloudData%>%
  select(-TOC)%>% ##Include whole dataset, as TOC start at 2009
  group_by(Year,Regime)%>%
  summarise(n = n())%>%
  filter(!is.na(Regime))%>%
  group_by(Year)%>%
  summarise(Percent = 100*n/(sum(n)))%>%
  ungroup()%>%
  add_row(Year = 2004, Percent = 0, .before = 22)



PercentRegime<-PercentRegime[seq(2, nrow(PercentRegime), by =2),] ##Filters out data to only look at percent invalid


PercentPlotRegime<-ggplot(PercentRegime)+
  geom_line(aes(Year, y = Percent), color = 'red', size = 2)+
  geom_point(aes(Year, y = Percent), fill = "red", shape = 21, size = 4)+
  PaperQuality()+scale_y_continuous(TeX('\\textbf{Percent Non-Linear (%)}'))+scale_x_continuous(breaks = seq(1994, 2017, 4))

ggsave(PercentPlotRegime, filename = "PercentPlotRegime.png", width =5, height = 5)


# Table 1. Regime Dunn Test ####

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

DunnResults<-DunnTestbyRegime%>%
  dplyr::select(c("Species","n1", "n2", "New Regime Median", "Old Regime Median", "p.adj", "p.adj.signif"))%>%
  arrange(desc(`New Regime Median`) )%>%
  dplyr::rename(`Old Regime n` = n1, `New Regime n` = n2)%>%
  mutate(Difference = `New Regime Median` - `Old Regime Median`)%>%
  dplyr::select(c(Species, `New Regime Median`,`Old Regime Median`, Difference, p.adj, `New Regime n`, `Old Regime n`))%>%
  mutate(across(where(is.numeric), signif, digits = 4))
  
  
RegimeTable<-ggtexttable(DunnResults, theme = ttheme("light", base_size = 16, colnames.style = colnames_style(face = "italic", size = 16)), rows = NULL)%>%
  tab_add_title(text = paste("Dunn-Test Comparison of Regimes", sep = "\n"), size = 17, face = "bold",padding = unit(1, "line"),
                just = "left")


ggsave(filename = "RegimeTable.png", plot = RegimeTable, width = 12, height = 8)




## Figure 14. Measured H+ vs Predicted H+ ####


MeasuredpHRegime<-ggplot(subset(AllCloudData, !is.na(Regime),
                                !is.na(LABPH)))+
  geom_point(aes(x = Anions - (Cations - Hplus*1e6),
                 y = Hplus*1e6, color = Regime),
             size = 3)+
  geom_smooth(aes(x = Anions - (Cations - Hplus*1e6),
                  y = Hplus*1e6), method = "lm",
              size  =2)+
  geom_abline(intercept = 0, slope = 1)+
  stat_regline_equation(size =6 ,aes(x = Anions - (Cations - Hplus*1e6),
                            y = Hplus*1e6, color = TOC,
                            label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
                        formula = y~x+0)+
  xlab(TeX("\\textbf{Predicted H^{+} ($\\mu$eq L^{-1})}"))+
  ylab(TeX("\\textbf{Measured H^{+} ($\\mu$eq L^{-1})}"))+
  facet_grid(Regime~.)+
  PaperQuality()+
  theme(strip.text = element_text(face = "bold", size = 14),
        legend.position = "none")

ggsave(filename = "MeasuredpHRegime.png", plot = MeasuredpHRegime,
       width = 8, height = 8)

#Figure 15 RPD vs TOC ####


## Regresion Lines

IonBalancelmLinear<-lm(TOC~IonBalance, 
                          data = subset(na.omit(AllCloudData), 
                                        Regime == "Linear"))
summary(IonBalancelmLinear)
AIC(IonBalancelmLinear)
BIC(IonBalancelmLinear)


IonBalanceSenLinear<-theilsen(TOC~IonBalance, 
                              data = subset(na.omit(AllCloudData), 
                              Regime == "Linear"), nboot = 600)
AIC(IonBalanceSenLinear)
BIC(IonBalanceSenLinear)


IonBalancelmNonLinear<-lm(TOC~IonBalance, data = subset(na.omit(AllCloudData), Regime == "Non-Linear"))
summary(IonBalancelmNonLinear)
AIC(IonBalancelmNonLinear)
BIC(IonBalancelmNonLinear)


IonBalanceSenNonLinear<-zyp.sen(TOC~IonBalance, 
                             data = subset(na.omit(AllCloudData), 
                                           Regime == "Non-Linear"))


TOCBalance<-ggplot(na.omit(AllCloudData),aes(x= Cations-Anions, y = TOC))+
  geom_point(aes(x= (Cations-Anions), y = TOC, color = LABPH))+
#  geom_smooth(method = 'rlm', aes(color = "RLM"))+
 # geom_smooth(method = "rlm")+
  geom_smooth(method = lm)+facet_grid(Regime~.)+
  geom_abline(aes(intercept = IonBalanceSenLinear$coefficients[[1]], 
              slope = IonBalanceSenLinear$coefficients[[2]]), 
              size = 2, color = "red", 
              data = AllCloudData%>%filter(Regime == "Linear"))+
  geom_ribbon(aes(ymax = IonBalanceSen$ci))
  geom_abline(aes(intercept = IonBalanceSenNonLinear$coefficients[[1]], 
                  slope = IonBalanceSenNonLinear$coefficients[[2]]), 
              size = 2, color = "red", 
              data = AllCloudData%>%filter(Regime == "Non-Linear"))+
  scale_color_gradientn(colors = rainbow(6))+
  labs(color = "pH")+ylim(c(0,4000))+ylab(TeX("\\textbf{TOC ($\\mu molC $ L^{-1}$)}"))+
  xlab(TeX("\\textbf{Cations - Anions ($\\mu eq $ L^{-1}$)}"))+
  geom_text(aes(x = 0, y = 3000), size = 6,
           data = AllCloudData%>%filter(Regime == "Non-Linear"),
  label = paste("y=",signif(IonBalanceSenNonLinear$coefficients[[1]],digits = 3),"+", 
                signif(IonBalanceSenNonLinear$coefficients[[2]],digits = 3),"x"), color = "red")+
  geom_text(aes(x = 0, y = 3000), size = 6,
            data = AllCloudData%>%filter(Regime == "Linear"),
            label = paste("y=",signif(IonBalanceSenLinear$coefficients[[1]],digits = 3),"+", 
                          signif(IonBalanceSenLinear$coefficients[[2]],digits = 3),"x"), color = "red")+
  #geom_smooth(method = "rlm")+
 # stat_regline_equation(method = "rlm")
  stat_regline_equation(size =6, label.x = -150, label.y = 2500,
                        aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")))+
  PaperQuality()+theme(legend.position = 'right', 
                       legend.title = element_text(size =14),
                       strip.text.y = element_text(face = "bold", size = 14))


ggsave(plot = TOCBalance,filename = "TOCBalance.png", width = 8, height = 8)

### With HCO3

IonBalancelmLinearHCO3<-lm(TOC~IonBalanceHCO3, 
                       data = subset(na.omit(AllCloudData), 
                                     Regime == "Linear"))
summary(IonBalancelmLinearHCO3)
AIC(IonBalancelmLinearHCO3)
BIC(IonBalancelmLinearHCO3)


IonBalanceSenLinearHCO3<-theilsen(TOC~IonBalanceHCO3, 
                             data = subset(na.omit(AllCloudData), 
                                           Regime == "Linear"), nboot = 10)
AIC(IonBalanceSenLinearHCO3)
BIC(IonBalanceSenLinearHCO3)


IonBalancelmNonLinear<-lm(TOC~IonBalance, data = subset(na.omit(AllCloudData), Regime == "Non-Linear"))
summary(IonBalancelmNonLinear)
AIC(IonBalancelmNonLinear)
BIC(IonBalancelmNonLinear)


IonBalanceSenNonLinearHCO3<-zyp.sen(TOC~IonBalanceHCO3, 
                             data = subset(na.omit(AllCloudData), 
                                           Regime == "Non-Linear"))
AIC(IonBalanceSenNonLinearHCO3)
BIC(IonBalanceSenNonLinearHCO3)


TOCBalanceHCO3<-ggplot(na.omit(AllCloudData),aes(x= Cations-Anions-HCO3Total, y = TOC))+
  geom_point(aes(x= (Cations-Anions-HCO3Total), y = TOC, color = LABPH))+
  #  geom_smooth(method = 'rlm', aes(color = "RLM"))+
  # geom_smooth(method = "rlm")+
  geom_smooth(method = lm)+facet_grid(Regime~.)+
  geom_abline(aes(intercept = IonBalanceSenLinearHCO3$coefficients[[1]], 
                  slope = IonBalanceSenLinearHCO3$coefficients[[2]]), 
              size = 2, color = "red", 
              data = AllCloudData%>%filter(Regime == "Linear"))+
  geom_abline(aes(intercept = IonBalanceSenNonLinearHCO3$coefficients[[1]], 
                  slope = IonBalanceSenNonLinearHCO3$coefficients[[2]]), 
              size = 2, color = "red", 
              data = AllCloudData%>%filter(Regime == "Non-Linear"))+
  scale_color_gradientn(colors = rainbow(6))+
  labs(color = "pH")+ylim(c(0,4000))+ylab(TeX("\\textbf{TOC ($\\mu molC $ L^{-1}$)}"))+
  xlab(TeX("\\textbf{Cations - Anions ($\\mu eq $ L^{-1}$)}"))+
  geom_text(aes(x = 0, y = 3000), size = 6,
            data = AllCloudData%>%filter(Regime == "Non-Linear"),
            label = paste("y=",signif(IonBalanceSenNonLinearHCO3$coefficients[[1]],digits = 3),"+", 
                          signif(IonBalanceSenNonLinearHCO3$coefficients[[2]],digits = 3),"x"), color = "red")+
  geom_text(aes(x = 0, y = 3000), size = 6,
            data = AllCloudData%>%filter(Regime == "Linear"),
            label = paste("y=",signif(IonBalanceSenLinearHCO3$coefficients[[1]],digits = 3),"+", 
                          signif(IonBalanceSenLinearHCO3$coefficients[[2]],digits = 3),"x"), color = "red")+
  #geom_smooth(method = "rlm")+
  # stat_regline_equation(method = "rlm")
  stat_regline_equation(size =6, label.x = -150, label.y = 2500,
                        aes(label = paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")))+
  PaperQuality()+theme(legend.position = 'right', 
                       legend.title = element_text(size =14),
                       strip.text.y = element_text(face = "bold", size = 14))


ggsave(plot = TOCBalanceHCO3,filename = "TOCBalanceHCO3.png", width = 8, height = 8)

## Figure 16. Inferred pH vs Measured pH ####

InferredpH<-AllCloudData%>%
  mutate(InferredHplus = CA+MG+1e6*10^(-LABPH))%>%
  group_by(Year)%>%
  summarise(across(where(is.numeric), median, na.rm = TRUE))%>%
  mutate(InferredpH = -log10(InferredHplus/1e6))%>% ## Convert from ueq
  ggplot()+
  geom_line(aes(x= Year, y = InferredpH, fill = "Inferred pH"), color = "#F8766D", size = 4)+
  geom_point(aes(x = Year, y = InferredpH, fill = "Inferred pH"), size = 4, shape = 21)+
  geom_line(aes(x= Year, y = -log10(Hplus), fill = "Measured pH"), size = 4, color = "#00C0B7")+
  geom_point(aes(x = Year, y = -log10(Hplus), fill = "Measured pH"), size = 4, shape = 21)+
  PaperQuality()+scale_x_continuous(breaks = seq(1994, 2017, by = 4))+
  ylab('pH')
  
  
ggsave(plot = InferredpH, filename = "InferredpH.png",
       width = 8, height = 8)



# Figure 17 Acidity from TOC ####


InferredpH<-ggplot(subset(AllCloudData, !is.na(TOC) & !is.na(Regime)))+
  geom_point(aes(y = SO4+NO3-NH4,
                 x = -log10((1e6*10^(-LABPH)+CA+MG)/1e6),
                 color = TOC))+
  scale_color_gradientn(colors = rainbow(6))+
  geom_line(aes(x= -log10((1e6*10^(-LABPH)+CA+MG)/1e6),
                y = ((1e6*10^(-LABPH)+CA+MG))),
            linetype = "dashed", size = 2)+
  facet_wrap(~Regime)+xlab("Inferred pH")+
  ylab(TeX("\\textbf{(SO_4^{-2}+NO_3^-) - NH_4^+ ($\\mu eq$ L^{-1} )}"))+
  labs(color = TeX("\\textbf{TOC $\\mu molC$ L^{-1}}"))+
  ggtitle(TeX("\\textbf{Inferred pH vs (SO_4^{-2}+NO_3^-) - NH_4^+}"))+PaperQuality()+
  theme(legend.position = "right", legend.title = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"))



MeasuredpH<-ggplot(subset(AllCloudData, !is.na(TOC) & !is.na(Regime)))+
  geom_point(aes(y = SO4+NO3-NH4,
                 x = LABPH,
                 color = TOC))+
  scale_color_gradientn(colors = rainbow(6))+
  geom_line(aes(x= LABPH,
                y = ((1e6*10^(-LABPH)))),
            linetype = "dashed", size = 2)+
  facet_wrap(~Regime)+xlab("Measured pH")+
  ylab(TeX("\\textbf{(SO_4^{-2}+NO_3^-) - NH_4^+ ($\\mu eq$ L^{-1}})"))+
  labs(color = TeX("\\textbf{TOC $\\mu molC$ L^{-1}}"))+
  ggtitle(TeX("\\textbf{Measured pH vs (SO_4^{-2}+NO_3^-) - NH_4^+}"))+PaperQuality()+
  theme(legend.position = "right", legend.title = element_text(size = 14),
        strip.text = element_text(size = 14, face = "bold"))


pHPlotsAll<-MeasuredpH+InferredpH&theme(legend.position = "right")
pHPlotsAll<-pHPlotsAll+plot_layout(guides = "collect",ncol = 1)

ggsave(plot = pHPlotsAll, filename = "TOCpHPlots.png", height = 10, width =10)
