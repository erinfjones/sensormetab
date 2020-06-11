### GAMUT Aquatic Sensor analysis
### Using data downloaded from commands in GAMUT_allsites_download.R
### Erin Fleming Jones, erinfjones3@gmail.com
### Created June 2018, updated Oct 2019

############# Setup ####
setwd("~/Box/Aanderud Lab/Projects/iUTAH/PR_CH Metabolism 2018/GAMUT Paper 2018/GAMUT 2018 streamMetabolizer_R/Sensor_downloads")
library(zoo)
library(corrgram)
library(readr)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(xts)
library(ggthemes)
library(gapminder)
library(magrittr)
library(data.table)
library(vegan)
theme_set(theme_few())

### Load functions
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
range01 <- function(x){(x-min(x, na.rm=TRUE))/(max(x, na.rm=TRUE)-min(x, na.rm=TRUE))}
fun1 <- function(x) rollmean(x, 2, na.pad=TRUE, align="left")
funfill<- function(x) na.approx(x, maxgap = 12, na.rm = FALSE)
funfillall<- function(x) na.approx(x, maxgap = 30, na.rm = FALSE)
# getSeason <- function(DATES) {
#   WS <- as.Date("2012-12-15", format = "%Y-%m-%d") # Winter Solstice
#   SE <- as.Date("2012-3-15",  format = "%Y-%m-%d") # Spring Equinox
#   SS <- as.Date("2012-6-15",  format = "%Y-%m-%d") # Summer Solstice
#   FE <- as.Date("2012-9-15",  format = "%Y-%m-%d") # Fall Equinox
#   
#   # Convert dates from any year to 2012 dates
#   d <- as.Date(strftime(DATES, format="2012-%m-%d"))
#   
#   ifelse (d >= WS | d < SE, "Winter",
#           ifelse (d >= SE & d < SS, "Spring",
#                   ifelse (d >= SS & d < FE, "Summer", "Fall")))
# }
getHydroperiod <- function(DATES) {
  Low <- as.Date("2012-10-1", format = "%Y-%m-%d") # Start of low
  Up <- as.Date("2012-3-31",  format = "%Y-%m-%d") # Start of increase
  Down <- as.Date("2012-6-30",  format = "%Y-%m-%d") # Start of decrease

  # Convert dates from any year to 2012 dates
  d <- as.Date(strftime(DATES, format="2012-%m-%d"))
  
  ifelse (d >= Low | d < Up, "Winter",
                  ifelse (d >= Up & d < Down, "SpringRunoff", "Summer"))
}

############# Load data ###########
### Sensor downloads by site (made in GAMUT_allsites_downloadXXXXXXX.R, edited in excel with precip/storm data)
Mendon<- read_csv("Mendon.csv", col_types = cols(NitrateMendon = col_double(), time = col_datetime(format = "%m/%d/%y %H:%M"), StormEvent= col_character()))
WaterLab<- read_csv("WaterLab.csv", col_types = cols( NitrateWaterLab = col_double(), time = col_datetime(format = "%m/%d/%y %H:%M")))
CH <- read_csv("CH.csv", col_types = cols(NitrateCH = col_double(),time = col_datetime(format = "%m/%d/%y %H:%M"), StormEvent= col_character()))
BJ <- read_csv("BJ.csv", col_types = cols(NitrateBJ = col_double(),time = col_datetime(format = "%m/%d/%y %H:%M")))
ARBR<- read_csv("ARBR.csv", col_types = cols( NitrateARBR = col_double(),time = col_datetime(format = "%m/%d/%Y %H:%M")))
FD <- read_csv("FD.csv", col_types = cols( NitrateFD = col_double(), time = col_datetime(format = "%m/%d/%y %H:%M"), StormEvent= col_character()))

### fill gaps of fewer than 12 obs (3 hours)
CH_fill=CH %>%
  mutate_each (funs(funfill), EXWTempCH, ODOCH, pHCH, QCH, NitrateCH, SpConCH)
Mendon_fill=Mendon %>%
  mutate_each (funs(funfill), EXWTempMendon, ODOMendon, pHMendon, QMendon, NitrateMendon, SpConMendon)
FD_fill=FD %>%
  mutate_each (funs(funfill), EXWTempFD, ODOFD, pHFD, QFD, NitrateFD, SpConFD)
ARBR_fill=ARBR %>%
  mutate_each (funs(funfill), EXWTempARBR, ODOARBR, pHARBR, QARBR, NitrateARBR, SpConARBR)
BJ_fill=BJ %>%
  mutate_each (funs(funfill), EXWTempBJ, ODOBJ, pHBJ, QBJ, NitrateBJ, SpConBJ)
WaterLab_fill=WaterLab %>%
  mutate_each (funs(funfill), EXWTempWaterLab, ODOWaterLab, pHWaterLab, QWaterLab, NitrateWaterLab, SpConWaterLab)

# #### Calculate Season
# Mendon_fill$season=getSeason(Mendon_fill$time)
# CH_fill$season=getSeason(CH_fill$time)
# FD_fill$season=getSeason(FD_fill$time)
# ARBR_fill$season=getSeason(ARBR_fill$time)
# BJ_fill$season=getSeason(BJ_fill$time)
# WaterLab_fill$season=getSeason(WaterLab_fill$time)

#calculate hydroperiod
Mendon_fill$Hydroperiod=getHydroperiod(Mendon_fill$time)
CH_fill$Hydroperiod=getHydroperiod(CH_fill$time)
FD_fill$Hydroperiod=getHydroperiod(FD_fill$time)
ARBR_fill$Hydroperiod=getHydroperiod(ARBR_fill$time)
BJ_fill$Hydroperiod=getHydroperiod(BJ_fill$time)
WaterLab_fill$Hydroperiod=getHydroperiod(WaterLab_fill$time)



############# hi-freq C-Q graphs ############
#### by hydroperiod/season
ggplot(Mendon_fill, aes(QMendon, NitrateMendon, color=as.integer(time))) + 
  facet_wrap(.~Hydroperiod, ncol = 4) + 
  geom_point(alpha=0.6, shape=1)+  guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  ylab("Nitrate (mg/L)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/hydroperiodCQ_Mendon.png", width=12, height=4)

ggplot(CH_fill, aes(QCH, NitrateCH, color=as.integer(time))) + 
  facet_wrap(.~Hydroperiod, ncol = 4) + 
  geom_point(alpha=0.6, shape=1)+  guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  ylab("Nitrate (mg/L)") + xlab("Discharge (m^3/s)" )ggsave("./figures/hydroperiodCQ_CH.png", width=12, height=4 )

ggplot(FD_fill, aes(QFD, NitrateFD, color=as.integer(time))) + 
  facet_wrap(.~Hydroperiod, ncol = 4) + 
  geom_point(alpha=0.6, shape=1)+  guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  ylab("Nitrate (mg/L)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/hydroperiodCQ_FD.png", width=12, height=4 )

ggplot(ARBR_fill, aes(QARBR, NitrateARBR, color=as.integer(time))) + 
  facet_wrap(.~Hydroperiod, ncol = 4) + 
  geom_point(alpha=0.6, shape=1)+  guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  ylab("Nitrate (mg/L)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/hydroperiodCQ_ARBR.png", width=12, height=4)

ggplot(WaterLab_fill, aes(QWaterLab, NitrateWaterLab, color=as.integer(time))) + 
  facet_wrap(.~Hydroperiod, ncol = 4) + 
  geom_point(alpha=0.6, shape=1)+  guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  ylab("Nitrate (mg/L)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/hydroperiodCQ_WaterLab.png", width=12, height=4 )

ggplot(BJ_fill, aes(QBJ, NitrateBJ, color=as.integer(time))) + 
  facet_wrap(.~Hydroperiod, ncol = 4) + 
  geom_point(alpha=0.6, shape=1)+  guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  ylab("Nitrate (mg/L)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/hydroperiodCQ_BJ.png", width=12, height=4 )


####### Cq graphs relativized ####
loganhydrocq=Mendon_fill%>%
  group_by(Hydroperiod) %>%
  mutate_each(funs(fun1), QMendon, NitrateMendon)%>%
  mutate_each(funs(range01), QMendon, NitrateMendon) 
ggplot(loganhydrocq, aes(QMendon, NitrateMendon, color=as.integer(time))) + 
  facet_wrap(.~Hydroperiod, ncol = 4) +
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="open"), size=0.5)  +
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Logan\n\n") + xlab("Relativized Discharge" )+
  theme(axis.title.x=element_blank())
ggsave("./figures/CQhydrop_M.png",width=6.5, height=2.5)

provohydrocq=CH_fill%>%
  group_by(Hydroperiod) %>%
  mutate_each(funs(fun1), QCH, NitrateCH)%>%
  mutate_each(funs(range01), QCH, NitrateCH) 
ggplot(provohydrocq, aes(QCH, NitrateCH, color=as.integer(time))) + 
  facet_wrap(.~Hydroperiod, ncol = 4) +
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="open"), size=0.5)  + 
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Provo\n\n") + xlab("Relativized Discharge" )
ggsave("./figures/CQhydrop_CH.png",width=6.5, height=2.5)

rbhydrocq=FD_fill%>%
  group_by(Hydroperiod) %>%
  mutate_each(funs(fun1), QFD, NitrateFD)%>%
  mutate_each(funs(range01), QFD, NitrateFD) 
ggplot(rbhydrocq, aes(QFD, NitrateFD, color=as.integer(time))) + 
  facet_wrap(.~Hydroperiod, ncol = 4) +  
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="open"), size=0.5)  + 
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Red Butte\nRelativized Nitrate") + xlab("Relativized Discharge" )+
  theme(axis.title.x=element_blank())
ggsave("./figures/CQhydrop_FD.png",width=6.5, height=2.5)


####### Storm C-Q graphs ################
Mqtile=Mendon_fill %>%
  filter(!is.na(StormEvent)) %>%
  group_by(StormEvent) %>%
  mutate_each(funs(fun1), QMendon, NitrateMendon)%>%
  mutate_each(funs(range01), QMendon, NitrateMendon) 
ggplot(Mqtile, aes(QMendon, NitrateMendon, color=as.integer(time))) + 
  facet_wrap(.~StormEvent, ncol = 4) + 
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="closed"))  + 
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Relativized Nitrate") + xlab("Relativized Discharge" )
ggsave("./figures/CQstorm_M.png",width=12, height=4)

FDqtile=FD_fill %>%
  filter(!is.na(StormEvent)) %>%
  group_by(StormEvent) %>%
  mutate_each(funs(fun1), QFD, NitrateFD) %>%
  mutate_each(funs(range01), QFD, NitrateFD)
ggplot(FDqtile, aes(QFD, NitrateFD, color=as.integer(time))) +
  facet_wrap(.~StormEvent, ncol = 4)+ 
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="closed"))  + 
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Relativized Nitrate") + xlab("Relativized Discharge" )
ggsave("./figures/CQstorm_FD.png",width=12, height=4)

CHqtile=CH_fill %>%
  filter(!is.na(StormEvent)) %>%
  group_by(StormEvent) %>%
  mutate_each(funs(fun1), QCH, NitrateCH) %>%
  mutate_each(funs(range01), QCH, NitrateCH) 
ggplot(CHqtile, aes(QCH, NitrateCH, color=as.integer(time))) + 
  facet_wrap(.~StormEvent, ncol = 4)+ 
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="closed"))  +
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Relativized Nitrate") + xlab("Relativized Discharge" )
ggsave("./figures/CQstorm_CH.png",width=12, height=4)

####### C-Q graphs by N Load and N conc ############
CH_fill$Nload=CH_fill$NitrateCH*CH_fill$QCH*1000 ##mg/s
FD_fill$Nload=FD_fill$NitrateFD*FD_fill$QFD*1000
Mendon_fill$Nload=Mendon_fill$NitrateMendon*Mendon_fill$QMendon*1000

ggplot(CH_fill, aes(rollmean(QCH,2,na.pad = TRUE, align="left"), rollmean(Nload,2,na.pad = TRUE, align="left"), color=as.integer(time))) +
  scale_color_gradientn(colours = rainbow(6)) + guides(color = "none")+
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="open"), alpha=0.7)+
  ylab("Nitrate load (mg/s)") + xlab("Discharge (m^3/s)" )
 ggsave("./figures/Nflux-QCH.png", width=10, height=6)
ggplot(FD_fill, aes(rollmean(QFD,2,na.pad = TRUE, align="left"), rollmean(Nload,2,na.pad = TRUE, align="left"), color=as.integer(time))) +
  scale_color_gradientn(colours = rainbow(6)) + guides(color = "none")+
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="open"), alpha=0.7)+
  ylab("Nitrate load (mg/s)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/Nflux-QFD.png", width=10, height=6)
ggplot(Mendon_fill, aes(rollmean(QMendon,2,na.pad = TRUE, align="left"), rollmean(Nload,2,na.pad = TRUE, align="left"), color=as.integer(time))) +
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="open"), alpha=0.7)+
  ylab("Nitrate load (mg/s)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/Nflux-QMendon.png", width=10, height=6)
ggplot(CH_fill, aes(rollmean(QCH,2,na.pad = TRUE, align="left"), rollmean(NitrateCH,2,na.pad = TRUE, align="left"), color=as.integer(time))) +
  scale_color_gradientn(colours = rainbow(6)) + guides(color = "none")+
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="open"), alpha=0.7)+
  ylab("Nitrate concentration (mg/L)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/N-QCH.png", width=10, height=6)
ggplot(FD_fill, aes(rollmean(QFD,2,na.pad = TRUE, align="left"), rollmean(NitrateFD,2,na.pad = TRUE, align="left"), color=as.integer(time))) +
  scale_color_gradientn(colours = rainbow(6)) + guides(color = "none")+
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="open"), alpha=0.7)+
  ylab("Nitrate concentration (mg/L)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/N-QFD.png", width=10, height=6)
ggplot(Mendon_fill, aes(rollmean(QMendon,2,na.pad = TRUE, align="left"), rollmean(NitrateMendon,2,na.pad = TRUE, align="left"), color=as.integer(time))) +
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  geom_path(arrow=arrow(length = unit(0.05, "inches"), type="open"), alpha=0.7)+
  ylab("Nitrate concentration (mg/L)") + xlab("Discharge (m^3/s)" )
ggsave("./figures/N-QMendon.png", width=10, height=6)



############# Diel NO3 ################
attributes(Mendon_fill$time)$tzone <- "America/Denver"  
attributes(Mendon_fill$time)
Mendon_fill$justtime <- as.POSIXct(as.ITime(Mendon_fill$time))
Mendon_fill$justtime=  format(as.POSIXct(Mendon_fill$justtime) ,format = "%H:%M") 
ggplot(Mendon_fill, aes(justtime, NitrateMendon, color=as.integer(time))) + 
  geom_point(alpha=0.6, shape=1,size=0.8)  + facet_grid(Hydroperiod~.)+
  scale_color_gradientn(colours = rainbow(6)) + 
  ylab("Nitrate (mg/L)") + theme(axis.title.x=element_blank()) + theme(axis.title.y=element_blank()) + 
  scale_x_discrete(breaks=c("00:00","06:00","12:00", "18:00"), labels=c("00:00","06:00","12:00", "18:00"))
ggsave("figures/dielNMendon.png",width=4, height=4)

attributes(CH_fill$time)$tzone <- "America/Denver"  
attributes(CH_fill$time)
CH_fill$justtime <-as.POSIXct( as.ITime(CH_fill$time))
CH_fill$justtime=  format(as.POSIXct(CH_fill$justtime) ,format = "%H:%M") 
ggplot(CH_fill, aes(justtime, NitrateCH, color=as.integer(time))) + 
  geom_point(alpha=0.6, shape=1,size=0.8)  + facet_grid(Hydroperiod~.)+
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Nitrate (mg/L)") + theme(axis.title.x=element_blank()) +theme(axis.title.y=element_blank()) + 
  scale_x_discrete(breaks=c("00:00","06:00","12:00", "18:00"), labels=c("00:00","06:00","12:00", "18:00"))
ggsave("figures/dielNCH.png",width=4, height=4)

attributes(ARBR_fill$time)$tzone <- "America/Denver"  
attributes(ARBR_fill$time)
ARBR_fill$justtime <- as.POSIXct( as.ITime(ARBR_fill$time))
ARBR_fill$justtime=  format(as.POSIXct(ARBR_fill$justtime) ,format = "%H:%M") 
ggplot(ARBR_fill, aes(justtime, NitrateARBR, color=as.integer(time))) + 
  geom_point(alpha=0.6, shape=1,size=0.8)  + facet_grid(Hydroperiod~.)+
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Nitrate (mg/L)") + theme(axis.title.x=element_blank()) +
  theme(strip.text.y = element_blank())+
  scale_x_discrete(breaks=c("00:00","06:00","12:00", "18:00"), labels=c("00:00","06:00","12:00", "18:00"))
ggsave("figures/dielNARBR.png",width=4, height=4)

attributes(WaterLab_fill$time)$tzone <- "America/Denver"  
attributes(WaterLab_fill$time)
WaterLab_fill$justtime <- as.POSIXct( as.ITime(WaterLab_fill$time))
WaterLab_fill$justtime=  format(as.POSIXct(WaterLab_fill$justtime) ,format = "%H:%M") 
ggplot(WaterLab_fill, aes(justtime, NitrateWaterLab, color=as.integer(time))) + 
  geom_point(alpha=0.6, shape=1,size=0.8)  + facet_grid(Hydroperiod~.)+
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Nitrate (mg/L)") + theme(axis.title.x=element_blank()) +
  theme(strip.text.y = element_blank())+
  scale_x_discrete(breaks=c("00:00","06:00","12:00", "18:00"), labels=c("00:00","06:00","12:00", "18:00"))
ggsave("figures/dielNWaterLab.png",width=4, height=4)

attributes(BJ_fill$time)$tzone <- "America/Denver"  
attributes(BJ_fill$time)
BJ_fill$justtime <- as.POSIXct( as.ITime(BJ_fill$time))
BJ_fill$justtime=  format(as.POSIXct(BJ_fill$justtime) ,format = "%H:%M") 
ggplot(BJ_fill, aes(justtime, NitrateBJ, color=as.integer(time))) + 
  geom_point(alpha=0.6, shape=1,size=0.8)  + facet_grid(Hydroperiod~.)+
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Nitrate (mg/L)") + theme(axis.title.x=element_blank()) +
  theme(strip.text.y = element_blank())+
  scale_x_discrete(breaks=c("00:00","06:00","12:00", "18:00"), labels=c("00:00","06:00","12:00", "18:00"))
ggsave("figures/dielNBJ.png",width=4, height=4)

attributes(FD_fill$time)$tzone <- "America/Denver"  
attributes(FD_fill$time)
FD_fill$justtime <- as.POSIXct( as.ITime(FD_fill$time))
FD_fill$justtime=  format(as.POSIXct(FD_fill$justtime) ,format = "%H:%M") 
ggplot(FD_fill, aes(justtime, NitrateFD, color=as.integer(time))) + 
  geom_point(alpha=0.6, shape=1,size=0.8)  + facet_grid(Hydroperiod~.)+
  scale_color_gradientn(colours = rainbow(6)) +  guides(color = "none")+
  ylab("Nitrate (mg/L)") + theme(axis.title.x=element_blank()) +theme(axis.title.y=element_blank()) + 
  scale_x_discrete(breaks=c("00:00","06:00","12:00", "18:00"), labels=c("00:00","06:00","12:00", "18:00"))
ggsave("figures/dielNFD.png",width=4, height=4)

####### Calculate and graph daily min/max ####
###WaterLab
means=as.data.frame(WaterLab_fill)
means$hour= cut(means$time, breaks="hour")
means=aggregate(NitrateWaterLab ~ hour, means, mean)
means$day <- as.Date(means$hour)
means$Time <- format(as.POSIXct(means$hour) ,format = "%H:%M:%S") 
midnight= means%>%
  filter(Time == '02:00:00')
noon= means%>%
  filter( Time== '14:00:00')
daily_NO3_WaterLab=merge(midnight[,c("day","NitrateWaterLab")], noon[,c("day","NitrateWaterLab")], by="day")
daily_NO3_WaterLab$dif=daily_NO3_WaterLab$NitrateWaterLab.x-daily_NO3_WaterLab$NitrateWaterLab.y
#daily_NO3_WaterLab$season=getSeason(daily_NO3_WaterLab$day)
daily_NO3_WaterLab$hydro=getHydroperiod(daily_NO3_WaterLab$day)
write.csv(daily_NO3_WaterLab,file="./output/dailyNO3WaterLab.csv")
hydroNO3_WaterLab= daily_NO3_WaterLab %>%
  group_by(hydro) %>%
  summarise(avgNO3=mean(dif, na.rm=TRUE), s=sd(dif, na.rm=TRUE), n=sum(!is.na(dif)))
write.csv(hydroNO3_WaterLab, file="./output/hydroNO3_WaterLab.csv")
png('./figures/NhistWaterLab.png', height=5, width=5, units="in", res=120)
hist( daily_NO3_WaterLab$dif,breaks=40, main="WaterLab", xlab="Daily change in [NO3]")
dev.off()
ts1=ggplot(daily_NO3_WaterLab, aes(day, dif))+ geom_hline(yintercept=0, linetype="dashed")+ geom_point()+ ylab("Logan Up") +theme(axis.title.x=element_blank())+ xlim(as.Date(c("2016-02-01","2017-02-01")))

###Mendon
means=as.data.frame(Mendon_fill)
means$hour= cut(means$time, breaks="hour")
means=aggregate(NitrateMendon ~ hour, means, mean)
means$day <- as.Date(means$hour)
means$Time <- format(as.POSIXct(means$hour) ,format = "%H:%M:%S") 
midnight= means%>%
  filter(Time == '02:00:00')
noon= means%>%
  filter( Time== '14:00:00')
daily_NO3_Mendon=merge(midnight[,c("day","NitrateMendon")], noon[,c("day","NitrateMendon")], by="day")
daily_NO3_Mendon$dif=daily_NO3_Mendon$NitrateMendon.x-daily_NO3_Mendon$NitrateMendon.y
daily_NO3_Mendon$season=getSeason(daily_NO3_Mendon$day)
daily_NO3_Mendon$hydro=getHydroperiod(daily_NO3_Mendon$day)
write.csv(daily_NO3_Mendon,file="./output/dailyNO3Mendon.csv")
hydroNO3_Mendon= daily_NO3_Mendon %>%
  group_by(hydro) %>%
  summarise(avgNO3=mean(dif, na.rm=TRUE), s=sd(dif, na.rm=TRUE), n=sum(!is.na(dif)))
write.csv(hydroNO3_Mendon, file="./output/hydroNO3_Mendon.csv")
png('./figures/NhistMendon.png', height=5, width=5, units="in", res=120)
hist( daily_NO3_Mendon$dif,breaks=100, main="Mendon", xlab="Daily change in [NO3]", xlim=c(-0.2,0.4))
dev.off()
ts2=ggplot(daily_NO3_Mendon, aes(day, dif))+ geom_hline(yintercept=0, linetype="dashed")+ geom_point()+ ylab("Logan Down")+ ylim(-0.2,0.2) + theme(axis.title.x=element_blank())+ xlim(as.Date(c("2016-02-01","2017-02-01")))

###ARBR
means=as.data.frame(ARBR_fill)
means$hour= cut(means$time, breaks="hour")
means=aggregate(NitrateARBR ~ hour, means, mean, na.omit=TRUE)
means$day <- as.Date(means$hour)
means$Time <- format(as.POSIXct(means$hour) ,format = "%H:%M:%S") 
midnight= means%>%
  filter(Time == '02:00:00')
noon= means%>%
  filter( Time== '14:00:00')
daily_NO3_ARBR=merge(midnight[,c("day","NitrateARBR")], noon[,c("day","NitrateARBR")], by="day")
daily_NO3_ARBR$dif=daily_NO3_ARBR$NitrateARBR.x-daily_NO3_ARBR$NitrateARBR.y
daily_NO3_ARBR$season=getSeason(daily_NO3_ARBR$day)
daily_NO3_ARBR$hydro=getHydroperiod(daily_NO3_ARBR$day)
write.csv(daily_NO3_ARBR,file="./output/dailyNO3ARBR.csv")
hydroNO3_ARBR= daily_NO3_ARBR %>%
  group_by(hydro) %>%
  summarise(avgNO3=mean(dif, na.rm=TRUE), s=sd(dif, na.rm=TRUE), n=sum(!is.na(dif)))
write.csv(hydroNO3_ARBR, file="./output/hydroNO3_ARBR.csv")
png('./figures/NhistARBR.png', height=5, width=5, units="in", res=120)
hist( daily_NO3_ARBR$dif,breaks=40, main="ARBR", xlab="Daily change in [NO3]")
dev.off()
ts3=ggplot(daily_NO3_ARBR, aes(day, dif))+ geom_hline(yintercept=0, linetype="dashed")+ geom_point()+ ylab("RB Up") + theme(axis.title.x=element_blank())+ylim(-0.05,0.05) + xlim(as.Date(c("2015-02-01","2016-02-01")))

###FD
means=as.data.frame(FD_fill)
means$hour= cut(means$time, breaks="hour")
means=aggregate(NitrateFD ~ hour, means, mean, na.omit=TRUE)
means$day <- as.Date(means$hour)
means$Time <- format(as.POSIXct(means$hour) ,format = "%H:%M:%S") 
midnight= means%>%
  filter(Time == '02:00:00')
noon= means%>%
  filter( Time== '14:00:00')
daily_NO3_FD=merge(midnight[,c("day","NitrateFD")], noon[,c("day","NitrateFD")], by="day")
daily_NO3_FD$dif=daily_NO3_FD$NitrateFD.x-daily_NO3_FD$NitrateFD.y
daily_NO3_FD$season=getSeason(daily_NO3_FD$day)
daily_NO3_FD$hydro=getHydroperiod(daily_NO3_FD$day)
write.csv(daily_NO3_FD,file="./output/dailyNO3FD.csv")
hydroNO3_FD= daily_NO3_FD %>%
  group_by(hydro) %>%
  summarise(avgNO3=mean(dif, na.rm=TRUE), s=sd(dif, na.rm=TRUE), n=sum(!is.na(dif)))
write.csv(hydroNO3_FD, file="./output/hydroNO3_FD.csv")
png('./figures/NhistFD.png', height=5, width=5, units="in", res=120)
hist( daily_NO3_FD$dif,breaks=40, main="FD", xlab="Daily change in [NO3]")
dev.off()
ts4=ggplot(daily_NO3_FD, aes(day, dif))+ geom_hline(yintercept=0, linetype="dashed")+ geom_point()+ ylab("RB Down") + theme(axis.title.x=element_blank())+ylim(-0.3,0.3) +xlim(as.Date(c("2015-02-01","2016-02-01")))
###BJ
means=as.data.frame(BJ_fill)
means$hour= cut(means$time, breaks="hour")
means=aggregate(NitrateBJ ~ hour, means, mean)
means$day <- as.Date(means$hour)
means$Time <- format(as.POSIXct(means$hour) ,format = "%H:%M:%S") 
midnight= means%>%
  filter(Time == '02:00:00')
noon= means%>%
  filter( Time== '14:00:00')
daily_NO3_BJ=merge(midnight[,c("day","NitrateBJ")], noon[,c("day","NitrateBJ")], by="day")
daily_NO3_BJ$dif=daily_NO3_BJ$NitrateBJ.x-daily_NO3_BJ$NitrateBJ.y
daily_NO3_BJ$season=getSeason(daily_NO3_BJ$day)
daily_NO3_BJ$hydro=getHydroperiod(daily_NO3_BJ$day)
write.csv(daily_NO3_BJ,file="./output/dailyNO3BJ.csv")
hydroNO3_BJ= daily_NO3_BJ %>%
  group_by(hydro) %>%
  summarise(avgNO3=mean(dif, na.rm=TRUE), s=sd(dif, na.rm=TRUE), n=sum(!is.na(dif)))
write.csv(hydroNO3_BJ, file="./output/hydroNO3_BJ.csv")
png('./figures/NhistBJ.png', height=5, width=5, units="in", res=120)
hist( daily_NO3_BJ$dif,breaks=40, main="BJ", xlab="Daily change in [NO3]")
dev.off()
ts5=ggplot(daily_NO3_BJ, aes(day, dif))+ geom_hline(yintercept=0, linetype="dashed")+ geom_point()+ ylab("Provo Up") + theme(axis.title.x=element_blank())+xlim(as.Date(c("2015-02-01","2016-02-01")))

###CH
means=as.data.frame(CH_fill)
means$hour= cut(means$time, breaks="hour")
means=aggregate(NitrateCH ~ hour, means, mean)
means$day <- as.Date(means$hour)
means$Time <- format(as.POSIXct(means$hour) ,format = "%H:%M:%S") 
midnight= means%>%
  filter(Time == '02:00:00')
noon= means%>%
  filter( Time== '14:00:00')
daily_NO3_CH=merge(midnight[,c("day","NitrateCH")], noon[,c("day","NitrateCH")], by="day")
daily_NO3_CH$dif=daily_NO3_CH$NitrateCH.x-daily_NO3_CH$NitrateCH.y
daily_NO3_CH$season=getSeason(daily_NO3_CH$day)
daily_NO3_CH$hydro=getHydroperiod(daily_NO3_CH$day)
write.csv(daily_NO3_CH,file="./output/dailyNO3CH.csv")
hydroNO3_CH= daily_NO3_CH %>%
  group_by(hydro) %>%
  summarise(avgNO3=mean(dif, na.rm=TRUE), s=sd(dif, na.rm=TRUE), n=sum(!is.na(dif)))
write.csv(hydroNO3_CH, file="./output/hydroNO3_CH.csv")
png('./figures/NhistCH.png', height=5, width=5, units="in", res=120)
hist( daily_NO3_CH$dif,breaks=40, main="CH", xlab="Daily change in [NO3]")
dev.off()
ts6=ggplot(daily_NO3_CH, aes(day, dif))+ geom_hline(yintercept=0, linetype="dashed")+ geom_point()+ ylab("Provo Down") + theme(axis.title.x=element_blank()) +xlim(as.Date(c("2015-02-01","2016-02-01")))

png("./figures/NO3difftimeseries.png",width = 7, height = 6, units = "in", res=120)
multiplot( ts1,ts2,ts3,ts4,ts5,ts6, cols = 1)
dev.off()


daily_NO3_FD<- read_csv("./output/dailyNO3FD.csv", col_types = cols(day = col_date(format = "%Y-%m-%d")))
daily_NO3_ARBR<- read_csv("./output/dailyNO3ARBR.csv", col_types = cols(day = col_date(format = "%Y-%m-%d")))


####### Stats#####
Mendonhydro=aovSufficient(avgNO3 ~ hydro, data=hydroNO3_Mendon) 
summary(Mendonhydro)
TukeyHSD(Mendonhydro)
WaterLabhydro=aovSufficient(avgNO3 ~ hydro, data=hydroNO3_WaterLab) 
summary(WaterLabhydro)
TukeyHSD(WaterLabhydro)
ARBRhydro=aovSufficient(avgNO3 ~ hydro, data=hydroNO3_ARBR) 
summary(ARBRhydro)
TukeyHSD(ARBRhydro)
FDhydro=aovSufficient(avgNO3 ~ hydro, data=hydroNO3_FD) 
summary(FDhydro)
TukeyHSD(FDhydro)
BJhydro=aovSufficient(avgNO3 ~ hydro, data=hydroNO3_BJ) 
summary(BJhydro)
TukeyHSD(BJhydro)
CHhydro=aovSufficient(avgNO3 ~ hydro, data=hydroNO3_CH) 
summary(CHhydro)
TukeyHSD(CHhydro)


p4=ggplot(hydroNO3_Mendon, aes(x=hydro, y=avgNO3)) + 
  geom_errorbar(aes(ymin=avgNO3-s, ymax=avgNO3+s), width=.2) +
  geom_line() +  geom_point() + 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+ #ylim(-0.2,0.5) +
  ylab("Mendon")
p1=ggplot(hydroNO3_WaterLab, aes(x=hydro, y=avgNO3)) + 
  geom_errorbar(aes(ymin=avgNO3-s, ymax=avgNO3+s), width=.2) +
  geom_line() +  geom_point()+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+ #ylim(-0.2,0.5) +
  ylab("WaterLab")
p2=ggplot(hydroNO3_ARBR, aes(x=hydro, y=avgNO3)) + 
  geom_errorbar(aes(ymin=avgNO3-s, ymax=avgNO3+s), width=.2) +
  geom_line() +  geom_point()+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())+ #ylim(-0.2,0.5) +
  ylab("ARBR")
p5=ggplot(hydroNO3_FD, aes(x=hydro, y=avgNO3)) + 
  geom_errorbar(aes(ymin=avgNO3-s, ymax=avgNO3+s), width=.2) +
  geom_line() +  geom_point()+ 
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank() )+ #ylim(-0.2,0.5) +
  ylab("FD")
p3=ggplot(hydroNO3_BJ, aes(x=hydro, y=avgNO3)) + 
  geom_errorbar(aes(ymin=avgNO3-s, ymax=avgNO3+s), width=.2) +
  geom_line() +  geom_point()+ xlab("Hydroperiod")+ #ylim(-0.1,0.1) +
  ylab("BJ")
p6=ggplot(hydroNO3_CH, aes(x=hydro, y=avgNO3)) + 
  geom_errorbar(aes(ymin=avgNO3-s, ymax=avgNO3+s), width=.2) +
  geom_line() +  geom_point()+ xlab("Hydroperiod")+ #ylim(-0.1,0.1) +
  ylab("CH")
png("./figures/NO3difhydro.png", width=6, height=6, units = "in",res = 120)
multiplot(p1,p2,p3,p4,p5,p6, cols=2)
dev.off()


###Calculate upstream integration###
# take discharge, divide by depth and channel width
ARBR_fill$velocityARBR=ARBR_fill$QARBR/ ARBR_fill$DepthARBR

############# Merge to watershed, combine daily metab ####
logan= merge(Mendon_fill, WaterLab_fill,by="time")
provo= merge (CH_fill, BJ_fill, by="time")
rb= merge(ARBR_fill, FD_fill, by="time")

#loganday= apply.daily(read.zoo(merge(Mendon_fill, WaterLab_fill,by="time"),index.column = "time"),mean,na.rm=TRUE)
#provoday= apply.daily(read.zoo(merge (CH_fill, BJ_fill, by="time"),index.column="time"),mean,na.rm=TRUE)
#rbday= apply.daily(read.zoo(merge(ARBR_fill, FD_fill, by="time"), index.column="time"),mean,na.rm=TRUE)

# write.zoo(loganday, "./output/daily_logan.csv")
# write.zoo(provoday, "./output/daily_provo.csv")
# write.zoo(rbday, "./output/daily_rb.csv")

## Use excel to remove pesky time stamp
loganday <- read_csv("./output/daily_logan.csv", col_types = cols(date = col_date(format = "%m/%d/%y")))
provoday<- read_csv("./output/daily_provo.csv", col_types = cols(date = col_date(format = "%m/%d/%y")))
rbday <- read_csv("./output/daily_rb.csv", col_types = cols(date = col_date(format = "%m/%d/%y")))

### Daily metab data (made in GAMUT_Metab_XXXXXX.R)
metabLR= read_csv("metabolism/LR_M_param_Data.csv", col_types=cols(
  X1 = col_double(),
  date = col_date(format = "%m/%d/%y"),
  GPP.daily = col_double(),
  GPP.daily.sd = col_double(),
  ER.daily = col_double(),
  ER.daily.sd = col_double(),
  K600.daily = col_double(),
  K600.daily.sd = col_double(),
  warnings = col_character(),
  errors = col_character()
))
metabPR= read_csv("metabolism/PR_CH_param_Data.csv", col_types=cols(
  X1 = col_double(),
  date = col_date(format = ""),
  GPP.daily = col_double(),
  GPP.daily.sd = col_double(),
  ER.daily = col_double(),
  ER.daily.sd = col_double(),
  K600.daily = col_double(),
  K600.daily.sd = col_double(),
  warnings = col_character(),
  errors = col_character()
))
metabRB= read_csv("metabolism/RB_FD_param_Data191115.csv", col_types=cols(
  X1 = col_double(),
  date = col_date(format = ""),
  GPP.daily = col_double(),
  GPP.daily.sd = col_double(),
  ER.daily = col_double(),
  ER.daily.sd = col_double(),
  K600.daily = col_double(),
  K600.daily.sd = col_double(),
  warnings = col_character(),
  errors = col_character()
))

provoday_metab=merge(provoday,metabPR, by="date")
loganday_metab=merge(loganday,metabLR, by="date")
rbday_metab=merge(rbday,metabRB, by="date")

provoday_metab$hydro=getHydroperiod(provoday_metab$date)
loganday_metab$hydro=getHydroperiod(loganday_metab$date)
rbday_metab$hydro=getHydroperiod(rbday_metab$date)

#### Merge N dif
loganday_metab= merge(loganday_metab, daily_NO3_WaterLab, by.x="date", by.y="day" )
loganday_metab= merge(loganday_metab, daily_NO3_Mendon, by.x="date", by.y="day" )
rbday_metab= merge(rbday_metab, daily_NO3_ARBR, by.x="date", by.y="day" )
rbday_metab= merge(rbday_metab, daily_NO3_FD, by.x="date", by.y="day" )
provoday_metab= merge(provoday_metab, daily_NO3_BJ, by.x="date", by.y="day" )
provoday_metab= merge(provoday_metab, daily_NO3_CH, by.x="date", by.y="day" )


write_csv(loganday_metab, "output/loganday_metab.csv")
write_csv(rbday_metab, "output/rbday_metab191115.csv")
write_csv(provoday_metab, "output/provoday_metab.csv")


#### Plot daily C-Q
ggplot(provoday_metab, aes(QCH, NitrateCH, color=as.integer(date))) + 
  geom_point()+  facet_wrap(.~hydro, ncol=3) + 
  scale_color_gradientn(colours = rainbow(6)) +
  ylab("Provo Nitrate (mg/L)") + xlab("Discharge (m^3/s)" ) 
ggplot(rbday_metab, aes(QFD, NitrateFD, color=as.integer(date))) +
  geom_point()+   facet_wrap(.~hydro, ncol=3) + 
  scale_color_gradientn(colours = rainbow(6)) +
  ylab("Red Butte Nitrate (mg/L)") + xlab("Discharge (m^3/s)" )
ggplot(loganday_metab, aes(QMendon, NitrateMendon, color=as.integer(date))) + 
  geom_point()+   facet_wrap(.~hydro, ncol=3) + 
  scale_color_gradientn(colours = rainbow(6)) +
  ylab("Logan Nitrate (mg/L)") + xlab("Discharge (m^3/s)" )


############# N-N graphs #############
#### Daily, Concentrations
ggplot(rbday_metab, aes(NitrateARBR, NitrateFD, color=as.integer(date)))   +
  facet_wrap(.~hydro, ncol=3) + guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", alpha=0.6)+
  ylab("Downstream Nitrate (mg/L)") + xlab("Upstream Nitrate (mg/L)")+
  geom_point()
ggsave("./figures/NvN_RB_daily.png",width=6, height=3, units = "in")

ggplot(provoday_metab, aes(NitrateCH, NitrateBJ, color=as.integer(date))) +
  facet_wrap(.~hydro, ncol=3) +guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", alpha=0.6)+
  geom_point()
ggsave("./figures/NvN_provo_daily.png",width=6, height=3, units = "in")

ggplot(loganday_metab, aes(NitrateMendon, NitrateWaterLab, color=as.integer(date))) +
  facet_wrap(.~hydro, ncol=3) +guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", alpha=0.6)+
  geom_point()
ggsave("./figures/NvN_logan_daily.png",width=6, height=3, units = "in")



#### hi-freq, hydroperiod #
ggplot(rb, aes(NitrateARBR, NitrateFD, color=as.integer(time)))   +
  facet_wrap(.~Hydroperiod.x, ncol=3) + guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", alpha=0.6)+
  ylab("Red Butte \n \nDownstream Nitrate (mg/L)") + xlab("Upstream Nitrate (mg/L)")+
  geom_point(size=0.5, alpha=0.6)+
  theme(axis.title.x=element_blank())
ggsave("./figures/NvN_RB_hydro.png",width=6, height=2.5, units = "in")

ggplot(provo, aes(NitrateBJ, NitrateCH, color=as.integer(time))) +
  facet_wrap(.~Hydroperiod.x, ncol=3) + guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", alpha=0.6)+
  ylab(c("Provo\n \n")) + xlab("Upstream Nitrate (mg/L)")+
  geom_point(size=0.5, alpha=0.6)
ggsave("./figures/NvN_provo_hydro.png",width=6, height=2.5, units = "in")

ggplot(logan, aes(NitrateWaterLab, NitrateMendon, color=as.integer(time))) +
  facet_wrap(.~Hydroperiod.x, ncol=3) + guides(color = "none")+
  scale_color_gradientn(colours = rainbow(6)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", alpha=0.6)+
  ylab("Logan\n \n") + xlab("Upstream Nitrate (mg/L)")+
  geom_point(size=0.5, alpha=0.6)+
  theme(axis.title.x=element_blank())
ggsave("./figures/NvN_logan_hydro.png",width=6, height=2.5, units = "in")

loganreachnitrate= loganday_metab%>%
  mutate (reachloadchange=(NitrateMendon*QMendon*86.4)-(NitrateWaterLab*QWaterLab*86.4), reachconcchange=NitrateMendon-NitrateWaterLab) %>%
  group_by(hydro) %>%
  summarise(sumreachloadchange=sum(reachloadchange, na.rm=TRUE), avgreachloadchange=mean(reachloadchange, na.rm=TRUE),
            avgreachconcchange=mean(reachconcchange, na.rm=TRUE), downnitrate= mean(NitrateMendon, na.rm=TRUE)) %>%
  mutate(perconcchange=100*avgreachconcchange/downnitrate)
write.csv(loganreachnitrate, "output/loganreachnitrate.csv")

rbreachnitrate= rbday_metab%>%
  mutate (reachloadchange=(NitrateFD*QFD*86.4)-(NitrateARBR*QARBR*86.4),  reachconcchange=NitrateFD-NitrateARBR) %>%
  group_by(hydro) %>%
  summarise(sumreachloadchange=sum(reachloadchange, na.rm=TRUE), avgreachloadchange=mean(reachloadchange, na.rm=TRUE),
            avgreachconcchange=mean(reachconcchange, na.rm=TRUE), downnitrate= mean(NitrateFD, na.rm=TRUE)) %>%
  mutate(perconcchange=100*avgreachconcchange/downnitrate)
write.csv(rbreachnitrate, "output/rbreachnitrate.csv")

provoreachnitrate= provoday_metab%>%
  mutate (reachloadchange=(NitrateCH*QCH*86.4)-(NitrateBJ*QBJ*86.4), reachconcchange=NitrateCH-NitrateBJ) %>%
  group_by(hydro) %>%
  summarise(sumreachloadchange=sum(reachloadchange, na.rm=TRUE), avgreachloadchange=mean(reachloadchange, na.rm=TRUE),
            avgreachconcchange=mean(reachconcchange, na.rm=TRUE), downnitrate= mean(NitrateCH, na.rm=TRUE)) %>%
  mutate(perconcchange=100*avgreachconcchange/downnitrate)
write.csv(provoreachnitrate, "output/provoreachnitrate.csv")

#### hi-freq, loads

CH_fill$Nload=CH_fill$NitrateCH*CH_fill$QCH*1000 ##mg/s
FD_fill$Nload=FD_fill$NitrateFD*FD_fill$QFD*1000
Mendon_fill$Nload=Mendon_fill$NitrateMendon*Mendon_fill$QMendon*1000

g=ggplot(GAMUTload, aes(log(NloadARBR), log(NloadFD), color=as.integer(DateTime))) +geom_point(shape=1) +ggtitle("Red Butte") + 
  scale_color_gradientn(colours = rainbow(6)) + theme(legend.position="none") + geom_segment(aes(x = -9, y = -9, xend = 0, yend = 0), color="black")
ggsave(g, file="RedButte.png",width=4, height=4)
g=ggplot(GAMUTload, aes(NloadBJ, NloadCH, color=as.integer(DateTime))) +geom_point(shape=1) + ggtitle("Provo") + 
  scale_color_gradientn(colours = rainbow(6)) + theme(legend.position="none") + geom_segment(aes(x = 0.5, y = 0.5, xend = 2, yend = 2), color="black")
ggsave(g, file="Provo.png",width=4, height=4)
g=ggplot(GAMUTload, aes(NloadWaterLab, NloadMendon, color=as.integer(DateTime))) +ggtitle("Logan")+geom_point(shape=1) +
  scale_color_gradientn(colours = rainbow(6)) + theme(legend.position="right") + geom_segment(aes(x = 0, y = 0, xend = 2, yend = 2), color="black")
ggsave(g, file="Logan.png",width=6, height=4)

############# Site specific time series graphs ###########

#### Check data correlation with correlogram ####
png('./figures/corrgramFD.png',width = 6, height = 6, units = "in", res=120)
corrgram(FD_fill[,3:13],use="pairwise.complete.obs",main="FD", lower.panel=panel.shade, upper.panel=panel.pts, order=TRUE)
dev.off()

png('./figures/corrgramARBR.png',width = 6, height = 6, units = "in", res=120)
corrgram(ARBR_fill[,3:13],use="pairwise.complete.obs",main="ARBR", lower.panel=panel.shade, upper.panel=panel.pts, order=TRUE)
dev.off()

png('./figures/corrgramCH.png',width = 6, height = 6, units = "in", res=120)
corrgram(CH_fill[,3:13],use="pairwise.complete.obs",main="CH", lower.panel=panel.shade, upper.panel=panel.pts, order=TRUE)
dev.off()

png('./figures/corrgramBJ.png',width = 6, height = 6, units = "in", res=120)
corrgram(BJ_fill[,3:13],use="pairwise.complete.obs",main="BJ", lower.panel=panel.shade, upper.panel=panel.pts, order=TRUE)
dev.off()

png('./figures/corrgramWaterLab.png',width = 6, height = 6, units = "in", res=120)
corrgram(WaterLab_fill[,3:13],use="pairwise.complete.obs",main="WaterLab", lower.panel=panel.shade, upper.panel=panel.pts, order=TRUE)
dev.off()

png('./figures/corrgramMendon.png',width = 6, height = 6, units = "in", res=120)
corrgram(Mendon_fill[,3:13],use="pairwise.complete.obs",main="Mendon", lower.panel=panel.shade, upper.panel=panel.pts, order=TRUE)
dev.off()


#### Time Series ####
## Logan
p1=ggplot(logan, aes(time)) + 
  geom_line(aes(y=EXWTempWaterLab), color="grey55")+ 
  geom_line(aes(y=EXWTempMendon) ) + 
  xlab("") +ylab("Temp\n(C)")  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2=ggplot(logan, aes(time)) +   
  geom_line(aes(y=ODOWaterLab),  color="grey55")+ 
  geom_line(aes(y=ODOMendon))+ 
  xlab("") +ylab("DO\n(mg/L)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3=ggplot(logan, aes(time)) +   
  geom_line(aes(y=QWaterLab),  color="grey55")+ 
  geom_line(aes(y=QMendon))+ 
  xlab("") +ylab("Discharge\n(cms)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4=ggplot(logan, aes(time)) +   
  geom_line(aes(y=NitrateWaterLab),  color="grey55")+ 
  geom_line(aes(y=NitrateMendon))+ 
  xlab("") +ylab("NO3\n(mg N/L)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p5=ggplot(logan, aes(time)) +   
  geom_line(aes(y=SpConWaterLab),  color="grey55")+ 
  geom_line(aes(y=SpConMendon))+ 
  xlab("") +ylab("Sp. Cond.\n(uS/cm)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6=ggplot(logan, aes(time)) +   
  geom_line(aes(y=fDOMWaterLab),  color="grey55")+ 
  geom_line(aes(y=fDOMMendon))+ 
  xlab("") +ylab("fDOM\n(QSU)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p7=ggplot(logan, aes(time)) +   
  geom_line(aes(y=TurbWaterLab),  color="grey55")+ 
  geom_line(aes(y=TurbMendon))+ 
  xlab("") +ylab("Turbidity\n(NTU)")

png("./figures/TSlogan.png",width = 4, height = 8, units = "in", res=120)
multiplot(p1,p2, p3,p4,p5,p6,p7,cols = 1)
dev.off()

#rb
p1=ggplot(rb, aes(time)) + 
  geom_line(aes(y=EXWTempARBR), color="grey55")+ 
  geom_line(aes(y=EXWTempFD) ) + 
  xlab("") +ylab("Temp\n(C)")  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2=ggplot(rb, aes(time)) +   
  geom_line(aes(y=ODOARBR),  color="grey55")+ 
  geom_line(aes(y=ODOFD))+ 
  xlab("") +ylab("DO\n(mg/L)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3=ggplot(rb, aes(time)) +   
  geom_line(aes(y=QARBR),  color="grey55")+ 
  geom_line(aes(y=QFD))+ 
  xlab("") +ylab("Discharge\n(cms)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4=ggplot(rb, aes(time)) +   
  geom_line(aes(y=NitrateARBR),  color="grey55")+ 
  geom_line(aes(y=NitrateFD))+ 
  xlab("") +ylab("NO3\n(mg N/L)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p5=ggplot(rb, aes(time)) +   
  geom_line(aes(y=SpConARBR),  color="grey55")+ 
  geom_line(aes(y=SpConFD))+ 
  xlab("") +ylab("Sp. Cond.\n(uS/cm)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6=ggplot(rb, aes(time)) +   
  geom_line(aes(y=fDOMARBR),  color="grey55")+ 
  geom_line(aes(y=fDOMFD))+ 
  xlab("") +ylab("fDOM\n(QSU)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p7=ggplot(rb, aes(time)) +   
  geom_line(aes(y=TurbARBR),  color="grey55")+ 
  geom_line(aes(y=TurbFD))+ 
  xlab("") +ylab("Turbidity\n(NTU)")

png("./figures/TSrb.png",width = 4, height =8, units = "in", res=120)
multiplot(p1,p2, p3,p4,p5,p6,p7,cols = 1)
dev.off()

#provo
p1=ggplot(provo, aes(time)) + 
  geom_line(aes(y=EXWTempBJ), color="grey55")+ 
  geom_line(aes(y=EXWTempCH) ) + 
  xlab("") +ylab("Temp\n(C)")  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2=ggplot(provo, aes(time)) +   
  geom_line(aes(y=ODOBJ),  color="grey55")+ 
  geom_line(aes(y=ODOCH))+ 
  xlab("") +ylab("DO\n(mg/L)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3=ggplot(provo, aes(time)) +   
  geom_line(aes(y=QBJ),  color="grey55")+ 
  geom_line(aes(y=QCH))+ 
  xlab("") +ylab("Discharge\n(cms)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4=ggplot(provo, aes(time)) +   
  geom_line(aes(y=NitrateBJ),  color="grey55")+ 
  geom_line(aes(y=NitrateCH))+ 
  xlab("") +ylab("NO3\n(mg N/L)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p5=ggplot(provo, aes(time)) +   
  geom_line(aes(y=SpConBJ),  color="grey55")+ 
  geom_line(aes(y=SpConCH))+ 
  xlab("") +ylab("Sp. Cond.\n(uS/cm)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6=ggplot(provo, aes(time)) +   
  geom_line(aes(y=fDOMBJ),  color="grey55")+ 
  geom_line(aes(y=fDOMCH))+ 
  xlab("") +ylab("fDOM\n(QSU)")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p7=ggplot(provo, aes(time)) +   
  geom_line(aes(y=TurbBJ),  color="grey55")+ 
  geom_line(aes(y=TurbCH))+ 
  xlab("") +ylab("Turbidity\n(NTU)")

png("./figures/TSprovo.png",width = 4, height =8, units = "in", res=120)
multiplot(p1,p2, p3,p4,p5,p6,p7,cols = 1)
dev.off()


############# Statistics ###########
Mendonstats= Mendon_fill %>%
  group_by(Hydroperiod) %>%
  summarise( maxNO3=max(NitrateMendon, na.rm=TRUE), minNO3=min(NitrateMendon, na.rm=TRUE), avgNO3=mean(NitrateMendon, na.rm = TRUE),sd= sd(NitrateMendon,na.rm=TRUE ))  %>%  
  mutate(range=maxNO3-minNO3 )
write.csv(Mendonstats, "output/Mendonstats.csv")
write.csv(summary(Mendon), "output/Mendonsummary.csv")

WaterLabstats= WaterLab_fill %>%
  group_by(Hydroperiod) %>%
  summarise( maxNO3=max(NitrateWaterLab, na.rm=TRUE), minNO3=min(NitrateWaterLab, na.rm=TRUE), avgNO3=mean(NitrateWaterLab, na.rm = TRUE),sd= sd(NitrateWaterLab,na.rm=TRUE ))  %>%  
  mutate(range=maxNO3-minNO3 )
write.csv( WaterLabstats, "output/WaterLabstats.csv")
write.csv(summary(WaterLab), "output/WaterLabsummary.csv")

ARBRstats= ARBR_fill %>%
  group_by(Hydroperiod) %>%
  summarise( maxNO3=max(NitrateARBR, na.rm=TRUE), minNO3=min(NitrateARBR, na.rm=TRUE), avgNO3=mean(NitrateARBR, na.rm = TRUE),sd= sd(NitrateARBR,na.rm=TRUE ))  %>%  
  mutate(range=maxNO3-minNO3 )
write.csv( ARBRstats, "output/ARBRstats.csv")
write.csv(summary(ARBR), "output/ARBRsummary.csv")

FDstats= FD_fill %>%
  group_by(Hydroperiod) %>%
  summarise( maxNO3=max(NitrateFD, na.rm=TRUE), minNO3=min(NitrateFD, na.rm=TRUE), avgNO3=mean(NitrateFD, na.rm = TRUE),sd= sd(NitrateFD,na.rm=TRUE ))  %>%  
  mutate(range=maxNO3-minNO3 )
write.csv( FDstats, "output/FDstats.csv")
write.csv(summary(FD), "output/FDsummary.csv")

BJstats= BJ_fill %>%
  group_by(Hydroperiod) %>%
  summarise( maxNO3=max(NitrateBJ, na.rm=TRUE), minNO3=min(NitrateBJ, na.rm=TRUE), avgNO3=mean(NitrateBJ, na.rm = TRUE),sd= sd(NitrateBJ,na.rm=TRUE ))  %>%  
  mutate(range=maxNO3-minNO3 )
write.csv( BJstats, "output/BJstats.csv")
write.csv(summary(BJ), "output/BJsummary.csv")

CHstats= CH_fill %>%
  group_by(Hydroperiod) %>%
  summarise( maxNO3=max(NitrateCH, na.rm=TRUE), minNO3=min(NitrateCH, na.rm=TRUE), avgNO3=mean(NitrateCH, na.rm = TRUE),sd= sd(NitrateCH,na.rm=TRUE ))  %>%  
  mutate(range=maxNO3-minNO3 )
write.csv( CHstats, "output/CHstats.csv")
write.csv(summary(CH), "output/CHsummary.csv")


#### Time series Stats ####
CH_TSstats= CH_fill %>%
  group_by(Hydroperiod) %>%
  summarise(Site="CH" ,
            maxEXWTemp=max(EXWTempCH, na.rm=TRUE), minEXWTemp=min(EXWTempCH, na.rm=TRUE), avgEXWTemp=mean(EXWTempCH, na.rm = TRUE),sdEXWTemp= sd(EXWTempCH,na.rm=TRUE ),
             maxODO=max(ODOCH, na.rm=TRUE), minODO=min(ODOCH, na.rm=TRUE), avgODO=mean(ODOCH, na.rm = TRUE),sdODO= sd(ODOCH,na.rm=TRUE ),
             maxQ=max(QCH, na.rm=TRUE), minQ=min(QCH, na.rm=TRUE), avgQ=mean(QCH, na.rm = TRUE),sdQ= sd(QCH,na.rm=TRUE ),
             maxpH=max(pHCH, na.rm=TRUE), minpH=min(pHCH, na.rm=TRUE), avgpH=mean(pHCH, na.rm = TRUE),sdpH= sd(pHCH,na.rm=TRUE ),
             maxNitrate=max(NitrateCH, na.rm=TRUE), minNitrate=min(NitrateCH, na.rm=TRUE), avgNitrate=mean(NitrateCH, na.rm = TRUE),sdNitrate= sd(NitrateCH,na.rm=TRUE ),
             maxSpCon=max(SpConCH, na.rm=TRUE), minSpCon=min(SpConCH, na.rm=TRUE), avgSpCon=mean(SpConCH, na.rm = TRUE),sdSpCon= sd(SpConCH,na.rm=TRUE ),
             maxfDOM=max(fDOMCH, na.rm=TRUE), minfDOM=min(fDOMCH, na.rm=TRUE), avgfDOM=mean(fDOMCH, na.rm = TRUE),sdfDOM= sd(fDOMCH,na.rm=TRUE ),
             maxTurb=max(TurbCH, na.rm=TRUE), minTurb=min(TurbCH, na.rm=TRUE), avgTurb=mean(TurbCH, na.rm = TRUE),sdTurb= sd(TurbCH,na.rm=TRUE )
             )  
BJ_TSstats= BJ_fill %>%
  group_by(Hydroperiod) %>%
  summarise( Site="BJ" ,
            maxEXWTemp=max(EXWTempBJ, na.rm=TRUE), minEXWTemp=min(EXWTempBJ, na.rm=TRUE), avgEXWTemp=mean(EXWTempBJ, na.rm = TRUE),sdEXWTemp= sd(EXWTempBJ,na.rm=TRUE ),
             maxODO=max(ODOBJ, na.rm=TRUE), minODO=min(ODOBJ, na.rm=TRUE), avgODO=mean(ODOBJ, na.rm = TRUE),sdODO= sd(ODOBJ,na.rm=TRUE ),
             maxQ=max(QBJ, na.rm=TRUE), minQ=min(QBJ, na.rm=TRUE), avgQ=mean(QBJ, na.rm = TRUE),sdQ= sd(QBJ,na.rm=TRUE ),
             maxpH=max(pHBJ, na.rm=TRUE), minpH=min(pHBJ, na.rm=TRUE), avgpH=mean(pHBJ, na.rm = TRUE),sdpH= sd(pHBJ,na.rm=TRUE ),
             maxNitrate=max(NitrateBJ, na.rm=TRUE), minNitrate=min(NitrateBJ, na.rm=TRUE), avgNitrate=mean(NitrateBJ, na.rm = TRUE),sdNitrate= sd(NitrateBJ,na.rm=TRUE ),
             maxSpCon=max(SpConBJ, na.rm=TRUE), minSpCon=min(SpConBJ, na.rm=TRUE), avgSpCon=mean(SpConBJ, na.rm = TRUE),sdSpCon= sd(SpConBJ,na.rm=TRUE ),
             maxfDOM=max(fDOMBJ, na.rm=TRUE), minfDOM=min(fDOMBJ, na.rm=TRUE), avgfDOM=mean(fDOMBJ, na.rm = TRUE),sdfDOM= sd(fDOMBJ,na.rm=TRUE ),
             maxTurb=max(TurbBJ, na.rm=TRUE), minTurb=min(TurbBJ, na.rm=TRUE), avgTurb=mean(TurbBJ, na.rm = TRUE),sdTurb= sd(TurbBJ,na.rm=TRUE )
  )   

WaterLab_TSstats= WaterLab_fill %>%
  group_by(Hydroperiod) %>%
  summarise( Site="WaterLab" ,
             maxEXWTemp=max(EXWTempWaterLab, na.rm=TRUE), minEXWTemp=min(EXWTempWaterLab, na.rm=TRUE), avgEXWTemp=mean(EXWTempWaterLab, na.rm = TRUE),sdEXWTemp= sd(EXWTempWaterLab,na.rm=TRUE ),
             maxODO=max(ODOWaterLab, na.rm=TRUE), minODO=min(ODOWaterLab, na.rm=TRUE), avgODO=mean(ODOWaterLab, na.rm = TRUE),sdODO= sd(ODOWaterLab,na.rm=TRUE ),
             maxQ=max(QWaterLab, na.rm=TRUE), minQ=min(QWaterLab, na.rm=TRUE), avgQ=mean(QWaterLab, na.rm = TRUE),sdQ= sd(QWaterLab,na.rm=TRUE ),
             maxpH=max(pHWaterLab, na.rm=TRUE), minpH=min(pHWaterLab, na.rm=TRUE), avgpH=mean(pHWaterLab, na.rm = TRUE),sdpH= sd(pHWaterLab,na.rm=TRUE ),
             maxNitrate=max(NitrateWaterLab, na.rm=TRUE), minNitrate=min(NitrateWaterLab, na.rm=TRUE), avgNitrate=mean(NitrateWaterLab, na.rm = TRUE),sdNitrate= sd(NitrateWaterLab,na.rm=TRUE ),
             maxSpCon=max(SpConWaterLab, na.rm=TRUE), minSpCon=min(SpConWaterLab, na.rm=TRUE), avgSpCon=mean(SpConWaterLab, na.rm = TRUE),sdSpCon= sd(SpConWaterLab,na.rm=TRUE ),
             maxfDOM=max(fDOMWaterLab, na.rm=TRUE), minfDOM=min(fDOMWaterLab, na.rm=TRUE), avgfDOM=mean(fDOMWaterLab, na.rm = TRUE),sdfDOM= sd(fDOMWaterLab,na.rm=TRUE ),
             maxTurb=max(TurbWaterLab, na.rm=TRUE), minTurb=min(TurbWaterLab, na.rm=TRUE), avgTurb=mean(TurbWaterLab, na.rm = TRUE),sdTurb= sd(TurbWaterLab,na.rm=TRUE )
  )   
Mendon_TSstats= Mendon_fill %>%
  group_by(Hydroperiod) %>%
  summarise( Site="Mendon" ,
             maxEXWTemp=max(EXWTempMendon, na.rm=TRUE), minEXWTemp=min(EXWTempMendon, na.rm=TRUE), avgEXWTemp=mean(EXWTempMendon, na.rm = TRUE),sdEXWTemp= sd(EXWTempMendon,na.rm=TRUE ),
             maxODO=max(ODOMendon, na.rm=TRUE), minODO=min(ODOMendon, na.rm=TRUE), avgODO=mean(ODOMendon, na.rm = TRUE),sdODO= sd(ODOMendon,na.rm=TRUE ),
             maxQ=max(QMendon, na.rm=TRUE), minQ=min(QMendon, na.rm=TRUE), avgQ=mean(QMendon, na.rm = TRUE),sdQ= sd(QMendon,na.rm=TRUE ),
             maxpH=max(pHMendon, na.rm=TRUE), minpH=min(pHMendon, na.rm=TRUE), avgpH=mean(pHMendon, na.rm = TRUE),sdpH= sd(pHMendon,na.rm=TRUE ),
             maxNitrate=max(NitrateMendon, na.rm=TRUE), minNitrate=min(NitrateMendon, na.rm=TRUE), avgNitrate=mean(NitrateMendon, na.rm = TRUE),sdNitrate= sd(NitrateMendon,na.rm=TRUE ),
             maxSpCon=max(SpConMendon, na.rm=TRUE), minSpCon=min(SpConMendon, na.rm=TRUE), avgSpCon=mean(SpConMendon, na.rm = TRUE),sdSpCon= sd(SpConMendon,na.rm=TRUE ),
             maxfDOM=max(fDOMMendon, na.rm=TRUE), minfDOM=min(fDOMMendon, na.rm=TRUE), avgfDOM=mean(fDOMMendon, na.rm = TRUE),sdfDOM= sd(fDOMMendon,na.rm=TRUE ),
             maxTurb=max(TurbMendon, na.rm=TRUE), minTurb=min(TurbMendon, na.rm=TRUE), avgTurb=mean(TurbMendon, na.rm = TRUE),sdTurb= sd(TurbMendon,na.rm=TRUE )
  )   
ARBR_TSstats= ARBR_fill %>%
  group_by(Hydroperiod) %>%
  summarise( Site="ARBR" ,
             maxEXWTemp=max(EXWTempARBR, na.rm=TRUE), minEXWTemp=min(EXWTempARBR, na.rm=TRUE), avgEXWTemp=mean(EXWTempARBR, na.rm = TRUE),sdEXWTemp= sd(EXWTempARBR,na.rm=TRUE ),
             maxODO=max(ODOARBR, na.rm=TRUE), minODO=min(ODOARBR, na.rm=TRUE), avgODO=mean(ODOARBR, na.rm = TRUE),sdODO= sd(ODOARBR,na.rm=TRUE ),
             maxQ=max(QARBR, na.rm=TRUE), minQ=min(QARBR, na.rm=TRUE), avgQ=mean(QARBR, na.rm = TRUE),sdQ= sd(QARBR,na.rm=TRUE ),
             maxpH=max(pHARBR, na.rm=TRUE), minpH=min(pHARBR, na.rm=TRUE), avgpH=mean(pHARBR, na.rm = TRUE),sdpH= sd(pHARBR,na.rm=TRUE ),
             maxNitrate=max(NitrateARBR, na.rm=TRUE), minNitrate=min(NitrateARBR, na.rm=TRUE), avgNitrate=mean(NitrateARBR, na.rm = TRUE),sdNitrate= sd(NitrateARBR,na.rm=TRUE ),
             maxSpCon=max(SpConARBR, na.rm=TRUE), minSpCon=min(SpConARBR, na.rm=TRUE), avgSpCon=mean(SpConARBR, na.rm = TRUE),sdSpCon= sd(SpConARBR,na.rm=TRUE ),
             maxfDOM=max(fDOMARBR, na.rm=TRUE), minfDOM=min(fDOMARBR, na.rm=TRUE), avgfDOM=mean(fDOMARBR, na.rm = TRUE),sdfDOM= sd(fDOMARBR,na.rm=TRUE ),
             maxTurb=max(TurbARBR, na.rm=TRUE), minTurb=min(TurbARBR, na.rm=TRUE), avgTurb=mean(TurbARBR, na.rm = TRUE),sdTurb= sd(TurbARBR,na.rm=TRUE )
  )   
FD_TSstats= FD_fill %>%
  group_by(Hydroperiod) %>%
  summarise( Site="FD" ,
             maxEXWTemp=max(EXWTempFD, na.rm=TRUE), minEXWTemp=min(EXWTempFD, na.rm=TRUE), avgEXWTemp=mean(EXWTempFD, na.rm = TRUE),sdEXWTemp= sd(EXWTempFD,na.rm=TRUE ),
             maxODO=max(ODOFD, na.rm=TRUE), minODO=min(ODOFD, na.rm=TRUE), avgODO=mean(ODOFD, na.rm = TRUE),sdODO= sd(ODOFD,na.rm=TRUE ),
             maxQ=max(QFD, na.rm=TRUE), minQ=min(QFD, na.rm=TRUE), avgQ=mean(QFD, na.rm = TRUE),sdQ= sd(QFD,na.rm=TRUE ),
             maxpH=max(pHFD, na.rm=TRUE), minpH=min(pHFD, na.rm=TRUE), avgpH=mean(pHFD, na.rm = TRUE),sdpH= sd(pHFD,na.rm=TRUE ),
             maxNitrate=max(NitrateFD, na.rm=TRUE), minNitrate=min(NitrateFD, na.rm=TRUE), avgNitrate=mean(NitrateFD, na.rm = TRUE),sdNitrate= sd(NitrateFD,na.rm=TRUE ),
             maxSpCon=max(SpConFD, na.rm=TRUE), minSpCon=min(SpConFD, na.rm=TRUE), avgSpCon=mean(SpConFD, na.rm = TRUE),sdSpCon= sd(SpConFD,na.rm=TRUE ),
             maxfDOM=max(fDOMFD, na.rm=TRUE), minfDOM=min(fDOMFD, na.rm=TRUE), avgfDOM=mean(fDOMFD, na.rm = TRUE),sdfDOM= sd(fDOMFD,na.rm=TRUE ),
             maxTurb=max(TurbFD, na.rm=TRUE), minTurb=min(TurbFD, na.rm=TRUE), avgTurb=mean(TurbFD, na.rm = TRUE),sdTurb= sd(TurbFD,na.rm=TRUE )
  )   

mergestats=t(BJ_TSstats%>% bind_rows(CH_TSstats,WaterLab_TSstats,Mendon_TSstats,ARBR_TSstats,FD_TSstats ))

write.csv(mergestats,"output/timeseriesstats.csv")


####Daily stats ##
summary(loganday)
summary(provoday)
summary(rbday)
####### Linear models ######
###N dif
###logan by season
loganday_metabwinter=subset(loganday_metab,hydro.y=="Winter")
loganday_metabsummer=subset(loganday_metab,hydro.y=="Summer")
loganday_metabspring=subset(loganday_metab,hydro.y=="SpringRunoff")

lmLRspring=lm(dif.y~ GPP.daily+ER.daily+NitrateWaterLab, data=loganday_metabspring)
summary(lmLRspring)
hist(lmLRspring$residuals)
lmLRsummer=lm(dif.y~ GPP.daily+ER.daily+NitrateWaterLab, data=loganday_metabsummer)
summary(lmLRsummer)
hist(lmLRsummer$residuals)
lmLRwinter=lm(dif.y~ GPP.daily+ER.daily+NitrateWaterLab, data=loganday_metabwinter)
summary(lmLRwinter)
hist(lmLRwinter$residuals)

###rb by season
rbday_metabwinter=subset(rbday_metab,hydro.y=="Winter")
rbday_metabsummer=subset(rbday_metab,hydro.y=="Summer")
rbday_metabspring=subset(rbday_metab,hydro.y=="SpringRunoff")

lmRBspring=lm(dif.y~ GPP.daily+ER.daily+NitrateARBR, data=rbday_metabspring)
summary(lmRBspring)
lmRBsummer=lm(dif.y~ GPP.daily+ER.daily+NitrateARBR, data=rbday_metabsummer)
summary(lmRBsummer)
lmRBwinter=lm(dif.y~ GPP.daily+ER.daily+NitrateARBR, data=rbday_metabwinter)
summary(lmRBwinter)

###provo by season
provoday_metabwinter=subset(provoday_metab,hydro.y=="Winter")
provoday_metabsummer=subset(provoday_metab,hydro.y=="Summer")
provoday_metabspring=subset(provoday_metab,hydro.y=="SpringRunoff")

lmprovospring=lm(dif.y~ GPP.daily+ER.daily+NitrateBJ, data=provoday_metabspring)
summary(lmprovospring)
lmprovosummer=lm(dif.y~ GPP.daily+ER.daily+NitrateBJ, data=provoday_metabsummer)
summary(lmprovosummer)
lmprovowinter=lm(dif.y~ GPP.daily+ER.daily+NitrateBJ, data=provoday_metabwinter)
summary(lmprovowinter)

#### sensor
lmNprovospring=lm(NitrateCH~ EXWTempCH+chlaCH+NitrateBJ+QCH+SpConCH+newprecip, data=provoday_metabspring)
summary(lmNprovospring)
lmNprovosummer=lm(NitrateCH~EXWTempCH+chlaCH+NitrateBJ+QCH+SpConCH+newprecip, data=provoday_metabsummer)
summary(lmNprovosummer)
lmNprovowinter=lm(NitrateCH~ EXWTempCH+chlaCH+NitrateBJ+QCH+SpConCH+newprecip, data=provoday_metabwinter)
summary(lmNprovowinter)

lmNloganspring=lm(NitrateMendon~ EXWTempMendon+chlaMendon+NitrateWaterLab+QMendon+SpConMendon+newprecip, data=loganday_metabspring)
summary(lmNloganspring)
lmNlogansummer=lm(NitrateMendon~EXWTempMendon+ chlaMendon+NitrateWaterLab+QMendon+SpConMendon+newprecip, data=loganday_metabsummer)
summary(lmNlogansummer)
lmNloganwinter=lm(NitrateMendon~ EXWTempMendon+chlaMendon+NitrateWaterLab+QMendon+SpConMendon+newprecip, data=loganday_metabwinter)
summary(lmNloganwinter)

lmNrbspring=lm(NitrateFD~ EXWTempFD+chlaFD+NitrateARBR+QFD+SpConFD+newprecip, data=rbday_metabspring)
summary(lmNrbspring)
lmNrbsummer=lm(NitrateFD~EXWTempFD+ chlaFD+NitrateARBR+QFD+SpConFD+newprecip, data=rbday_metabsummer)
summary(lmNrbsummer)
lmNrbwinter=lm(NitrateFD~ EXWTempFD+chlaFD+NitrateARBR+QFD+SpConFD+newprecip, data=rbday_metabwinter)
summary(lmNrbwinter)


####### Corrgrams ########
loganmetab=loganday_metab%>%
  mutate(metabratio=GPP.daily/ER.daily)%>%
  transmute_at (c("GPP.daily","ER.daily","dif.y","dif.x","metabratio","EXWTempMendon","ODOMendon", "pHMendon", "QMendon", "NitrateMendon", "SpConMendon","BGAMendon", "TurbMendon", "chlaMendon","fDOMMendon","newprecip",
                  "EXWTempWaterLab","ODOWaterLab", "pHWaterLab", "QWaterLab", "NitrateWaterLab", "SpConWaterLab","BGAWaterLab", "TurbWaterLab", "fDOMWaterLab"), funfill)

png('./figures/corrgramLoganmetab.png',width = 15, height = 15, units = "in", res=120)
corrgram(loganmetab,main="Logan", lower.panel=panel.shade, upper.panel=panel.pts, order=TRUE)
dev.off()

rbmetab=rbday_metab%>%
  mutate(metabratio=GPP.daily/ER.daily)%>%
  transmute_at (c("GPP.daily","ER.daily","dif.y","dif.x","metabratio","EXWTempFD","ODOFD", "pHFD", "QFD", "NitrateFD", "SpConFD","BGAFD", "TurbFD", "chlaFD","fDOMFD","newprecip",
                  "EXWTempARBR","ODOARBR", "pHARBR", "QARBR", "NitrateARBR", "SpConARBR","BGAARBR", "TurbARBR","chlaARBR", "fDOMARBR"), funfill)
png('./figures/corrgramrbmetab.png',width = 15, height = 15, units = "in", res=120)
corrgram(rbmetab,main="Red Butte", lower.panel=panel.shade, upper.panel=panel.pts, order=TRUE)
dev.off()

provometab=provoday_metab%>%
  mutate(metabratio=GPP.daily/ER.daily)%>%
  transmute_at (c("GPP.daily","ER.daily","dif.y","dif.x","metabratio","EXWTempCH","ODOCH", "pHCH", "QCH", "NitrateCH", "SpConCH","BGACH", "TurbCH", "chlaCH","fDOMCH","newprecip",
                  "EXWTempBJ","ODOBJ", "pHBJ", "QBJ", "NitrateBJ", "SpConBJ","BGABJ", "TurbBJ","chlaBJ", "fDOMBJ"), funfill)
png('./figures/corrgramprovometab.png',width = 15, height = 15, units = "in", res=120)
corrgram(provometab,main="Provo", lower.panel=panel.shade, upper.panel=panel.pts, order=TRUE)
dev.off()

#make the plot
library(GGally)
p = ggpairs(data = WaterLab_fill, columns = 3:12, aes(color = Hydroperiod, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
# #set the colors
# for(i in 1:p$nrow) {
#   for(j in 1:p$ncol){
#     p[i,j] <- p[i,j] + 
#       scale_fill_manual(values=c("deepskyblue","darkorange", "black")) +
#       scale_color_manual(values=c("deepskyblue","darkorange", "black"))  
#   }
# }
ggsave("figures/corWaterLab.png", plot = p, width = 14, height = 14)

p = ggpairs(data = Mendon_fill, columns = 3:12, aes(color = Hydroperiod, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
ggsave("figures/corMendon.png", plot = p, width = 14, height = 14)

p = ggpairs(data = BJ_fill, columns = 3:12, aes(color = Hydroperiod, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
ggsave("figures/corBJ.png", plot = p, width = 14, height = 14)

p = ggpairs(data = CH_fill, columns = 3:12, aes(color = Hydroperiod, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
ggsave("figures/corCH.png", plot = p, width = 14, height = 14)

p = ggpairs(data = ARBR_fill, columns = 3:12, aes(color = Hydroperiod, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
ggsave("figures/corARBR.png", plot = p, width = 14, height = 14)

p = ggpairs(data = FD_fill, columns = 3:12, aes(color = Hydroperiod, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
ggsave("figures/corFD.png", plot = p, width = 14, height = 14)

#select the data columns for all the sites

#calculate the spearman correlation matrix and extract the r values
ARBR_cor = ARBR_fill[3:12]
ARBRcor = rcorr(as.matrix(ARBR_cor), type = "spearman")
ARBRcor = data.frame(ARBRcor$r)
write.csv(ARBRcor, file = "output/ARBRcor.csv")

FD_cor = FD_fill[3:12]
FDcor = rcorr(as.matrix(FD_cor), type = "spearman")
FDcor = data.frame(FDcor$r)
write.csv(FDcor, file = "output/FDcor.csv")

WaterLab_cor = WaterLab_fill[3:12]
WaterLabcor = rcorr(as.matrix(WaterLab_cor), type = "spearman")
WaterLabcor = data.frame(WaterLabcor$r)
write.csv(WaterLabcor, file = "output/WaterLabcor.csv")

Mendon_cor = Mendon_fill[3:12]
Mendoncor = rcorr(as.matrix(Mendon_cor), type = "spearman")
Mendoncor = data.frame(Mendoncor$r)
write.csv(Mendoncor, file = "output/Mendoncor.csv")

CH_cor = CH_fill[3:12]
CHcor = rcorr(as.matrix(CH_cor), type = "spearman")
CHcor = data.frame(CHcor$r)
write.csv(CHcor, file = "output/CHcor.csv")

BJ_cor = BJ_fill[3:12]
BJcor = rcorr(as.matrix(BJ_cor), type = "spearman")
BJcor = data.frame(BJcor$r)
write.csv(BJcor, file = "output/BJcor.csv")


##Metab re run #####
RB_FD_param_Data191115 <- read_csv("metabolism/RB_FD_param_Data191115.csv",    col_types = cols(date = col_date(format = "%Y-%m-%d")))
RB_FD_param_Data <- read_csv("metabolism/RB_FD_param_Data.csv",  col_types = cols(date = col_date(format = "%Y-%m-%d")))
daily_NO3_FD<- read_csv("./output/dailyNO3FD.csv", col_types = cols(day = col_date(format = "%Y-%m-%d")))
daily_NO3_ARBR<- read_csv("./output/dailyNO3ARBR.csv", col_types = cols(day = col_date(format = "%Y-%m-%d")))
daily_DO_FD<- read_csv("./output/dailyDOFD.csv", col_types = cols(day = col_date(format = "%Y-%m-%d")))
rbday_metab=read_csv("./output/rbday_metab191115.csv",col_types = cols(date = col_date(format = "%m/%d/%y")))
rbday_metab_osat=merge(rbday_metab,daily_DO_FD, by.x = "date", by.y="day" )
p = ggpairs(data = rbday_metab_osat, columns = 15:41, aes(color = hydro.y, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
p
ggsave("figures/corFD.png", plot = p, width = 14, height = 14)

####### Calculate sp discharge, mean median and standard deviation ####
### Calculate Specific Discharge ###
# GAMUTfill$SpDisFD<-GAMUTfill$QFD*1384.58602
# GAMUTfill$SpDisARBR<-GAMUTfill$QARBR*1683.35491
# GAMUTfill$SpDisBJ<-GAMUTfill$QBJ*46.8911408
# GAMUTfill$SpDisCH<-GAMUTfill$QCH*40.4528045
# GAMUTfill$SpDisWaterLab<-GAMUTfill$QWaterLab*56.6766615
# GAMUTfill$SpDisMendon<-GAMUTfill$QMendon*16.3855282
# write.csv(GAMUTfill, "GAMUTfill.csv")
for (i in 2:49){
  c=colnames(GAMUTfill)[i] 
  print(c)
}

for (i in GAMUTfill[2:49]){
  c=mean(i, na.rm = TRUE)
  print(c)
}

for (i in GAMUTfill[2:49]){
  c=median(i, na.rm = TRUE)
  print(c)
}

for (i in GAMUTfill[2:49]){
  c=sd(i, na.rm = TRUE)
  print(c)
}

for (i in GAMUTfill[2:49]){
  c=sum(i, na.rm = TRUE)
  print(c)
}

for (i in GAMUTfill[2:49]){
  c=sum(!is.na(i)) 
  print(c)
}

for (i in 3:96){
  c=colnames(GAMUTload)[i] 
  print(c)
}

for (i in GAMUTload[3:96]){
  c=mean(i, na.rm = TRUE)
  print(c)
}

for (i in GAMUTload[3:96]){
  c=median(i, na.rm = TRUE)
  print(c)
}

for (i in GAMUTload[3:96]){
  c=sd(i, na.rm = TRUE)
  print(c)
}

# Not used ####
daily_NO3 <- Mendon_fill %>%
  mutate(day = as.Date(time, format="%Y-%m-%d")) %>%
  group_by(day) %>% # group by the day column
  summarise( maxNO3=max(NitrateMendon, na.rm=TRUE), minNO3=min(NitrateMendon, na.rm=TRUE), avgNO3=mean(NitrateMendon, na.rm = TRUE))  %>%  
  mutate(maxNO3 = ifelse(is.infinite(maxNO3), NA, maxNO3), minNO3 = ifelse(is.infinite(minNO3), NA, minNO3), avgNO3=ifelse(is.infinite(avgNO3), NA, avgNO3))
daily_NO3$dif= daily_NO3$maxNO3 - daily_NO3$minNO3
daily_NO3$season=getSeason(daily_NO3$day)
daily_NO3$hydro=getHydroperiod(daily_NO3$day)

GAMUTload <- read_csv("GAMUTload.csv", col_types = cols(DateTime = col_datetime(format = "%m/%d/%Y %H:%M"),  NitrateARBR = col_double(), NitrateBJ = col_double(), NitrateCH = col_double(), NitrateFD = col_double(), NitrateMendon = col_double(),
                                                        NitrateWaterLab = col_double(),  NloadARBR = col_double(), NloadBJ = col_double(),  NloadCH = col_double(), NloadFD = col_double(), NloadMendon = col_double(), NloadWaterLab = col_double()))
GAMUTload$integer<-as.integer(GAMUTload$DateTime)

# ggplot(GAMUTload, aes(log(TurbFD), log(TurbARBR), color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(log(TurbCH), log(TurbBJ), color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(log(TurbMendon), log(TurbWaterLab), color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# 
# ggplot(GAMUTload, aes(SpConFD, SpConARBR, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(SpConCH, SpConBJ, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(SpConMendon, SpConWaterLab, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# 
# ggplot(GAMUTload, aes(NitrateFD, NitrateARBR,  color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(NitrateCH, NitrateBJ, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(NitrateMendon, NitrateWaterLab, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# 
# ggplot(GAMUTload, aes(fDOMFD, fDOMARBR, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(fDOMCH, fDOMBJ,color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(fDOMMendon, fDOMWaterLab, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# 
# p2=ggplot(GAMUTload, aes(log(PloadFD),log(PloadARBR),  color=as.integer(DateTime))) +geom_point(shape=1) + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE) + theme(legend.position="none")
# p4=ggplot(GAMUTload, aes(log(PloadCH), log(PloadBJ), color=as.integer(DateTime))) +geom_point(shape=1) + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE) + theme(legend.position="none")
# ggplot(GAMUTload, aes(log(PloadMendon), log(PloadWaterLab), color=as.integer(DateTime))) +geom_point(shape=1)  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE) + theme(legend.position="none")

g=ggplot(GAMUTload, aes(SpConARBR, SpConFD, color=as.integer(DateTime))) +geom_point(shape=1) +ggtitle("Red Butte") + 
  scale_color_gradientn(colours = rainbow(6)) + theme(legend.position="none") + geom_segment(aes(x = 300, y = 300, xend = 800, yend = 800), color="black")
ggsave(g, file="RedButteSp.png",width=4, height=4)
g=ggplot(GAMUTload, aes(SpConBJ, SpConCH, color=as.integer(DateTime))) +geom_point(shape=1) + ggtitle("Provo") + 
  scale_color_gradientn(colours = rainbow(6)) + theme(legend.position="none") + geom_segment(aes(x = 160, y = 160, xend = 225, yend = 225), color="black")
ggsave(g, file="ProvoSp.png",width=4, height=4)
g=ggplot(GAMUTload, aes(SpConWaterLab, SpConMendon, color=as.integer(DateTime))) +ggtitle("Logan")+geom_point(shape=1) +
  scale_color_gradientn(colours = rainbow(6)) + theme(legend.position="right") + geom_segment(aes(x = 300, y = 300, xend = 430, yend = 430), color="black")
ggsave(g, file="LoganSp.png",width=6, height=4)




# ggplot(GAMUTload, aes(CloadFD,CloadARBR,  color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(CloadCH, CloadBJ,color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(CloadMendon, CloadWaterLab, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# 
# ggplot(GAMUTload, aes(condloadFD, condloadARBR,  color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(condloadCH, condloadBJ, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(condloadMendon, condloadWaterLab, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# 
# ggplot(GAMUTload, aes(log(QFD), log(QARBR), color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(QCH, QBJ, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(6)) + geom_smooth(method = "lm", se = FALSE)
# ggplot(GAMUTload, aes(QMendon, QWaterLab, color=as.integer(DateTime))) +geom_point()  + scale_color_gradientn(colours = rainbow(12)) + geom_smooth(method = "lm", se = FALSE)

##### calculate load and flux ##

# GAMUTload$SpDisFD<-GAMUTload$QFD*1384.58602
# GAMUTload$SpDisARBR<-GAMUTload$QARBR*1683.35491
# GAMUTload$SpDisBJ<-GAMUTload$QBJ*46.8911408
# GAMUTload$SpDisCH<-GAMUTload$QCH*40.4528045
# GAMUTload$SpDisWaterLab<-GAMUTload$QWaterLab*56.6766615
# GAMUTload$SpDisMendon<-GAMUTload$QMendon*16.3855282


#### rolling mean time series###
ggplot(GAMUTload) + geom_point(aes(DateTime, rollmean(SpDisARBR, k = 4, fill = NA, align = "left"))) #+ geom_line(aes(y=SpDisARBR), color="red" )+ geom_line(aes(y=SpDisFD), color="pink" )+
  geom_point(aes(y=SpDisBJ), color="blue" )+ geom_point(aes(y=SpDisCH), color="purple" )+ geom_line(aes(y=SpDisWaterLab), color="green" ) + geom_line(aes(y=SpDisMendon), color="yellow" )

ggplot(GAMUTload, aes(DateTime)) + geom_line(aes(y=rollapply(SpDisARBR,FUN=mean,12,fill=NA,align="right")), color="red" ) + geom_line(aes(y=rollapply(SpDisFD,FUN=mean,12,fill=NA,align="right")), color="pink" ) +
  geom_point(aes(y=rollapply(SpDisBJ,FUN=mean,1,fill=NA,align="right")), color="blue" ) + geom_point(aes(y=rollapply(SpDisCH,FUN=mean,1,fill=NA,align="right")), color="purple" ) + 
  geom_line(aes(y=rollapply(SpDisWaterLab,FUN=mean,12,fill=NA,align="right")), color="green" ) + geom_line(aes(y=rollapply(SpDisMendon,FUN=mean,12,fill=NA,align="right")), color="yellow" )

# GAMUTload$rollSpDisFD<-rollapply(GAMUTload$SpDisFD,FUN=mean,4,fill = NA, align = "right")
# GAMUTload$rollSpDisARBR<-rollmean(GAMUTload$SpDisARBR,4)
# GAMUTload$rollSpDisFD<-rollmean(GAMUTload$SpDisFD,4)
# GAMUTload$rollSpDisFD<-rollmean(GAMUTload$SpDisFD,4)

p1=ggplot(GAMUTload, aes(DateTime)) + geom_line(aes(y=rollapply(NloadFD,FUN=mean,4,fill=NA,align="right")),color="red") + geom_line(aes(y=rollapply(NloadARBR,FUN=mean,4,fill=NA,align="right")), color="pink")
p2=ggplot(GAMUTload, aes(DateTime)) + geom_point(aes(y=rollapply(NloadBJ,FUN=mean,1,fill=NA,align="right")), color="blue",size=0.8) + geom_point(aes(y=rollapply(NloadCH,FUN=mean,1,fill=NA,align="right")), color="purple",size=0.8)
p3=ggplot(GAMUTload, aes(DateTime)) + geom_line(aes(y=rollapply(NloadWaterLab,FUN=mean,4,fill=NA,align="right")), color="darkgreen") + geom_line(aes(y=rollapply(NloadMendon,FUN=mean,4,fill=NA,align="right")), color="springgreen")
multiplot(p1, p2, p3, cols = 1)

ggplot(NitrateWaterLab, aes(time)) + geom_line(aes(y=DataValue))
ggplot(NitrateMendon, aes(time)) + geom_line(aes(y=DataValue))




#### Monthly Stats ###
GAMUTfill=read_csv("GAMUTfill.csv", col_types = cols(DateTime = col_datetime(format = "%m/%d/%y %H:%M"), NitrateARBR = col_double(), NitrateBJ = col_double(),  NitrateCH = col_double(), NitrateFD = col_double(), NitrateMendon = col_double(), NitrateWaterLab = col_double(), NloadARBR = col_double(), 
                                                     NloadBJ = col_double(), NloadCH = col_double(), NloadFD = col_double(), NloadMendon = col_double(), NloadWaterLab = col_double(), X1 = col_skip()))
GAMUTfillzoo=read.zoo(GAMUTfill,index.column="DateTime") 
GAMUTfillmonth<-apply.monthly(GAMUTfillzoo,mean)   ###calculate monthly mean
write.csv(GAMUTfillmonth,"GAMUTfillmonthmeans.csv")

## GAMUTstats: Edited headings in Excel on GAMUTfill to have parameter:watershed:position format
stats=read_csv("GAMUTstats.csv", col_types = cols(GAMUTload.DateTime = col_datetime(format = "%m/%d/%y %H:%M"), `Nitrate:LR:High` = col_double(),  `Nitrate:LR:Low` = col_double(),  `Nitrate:PR:High` = col_double(),  `Nitrate:PR:Low` = col_double(), `Nitrate:RB:High` = col_double(),  
                                                  `Nitrate:RB:Low` = col_double(), `Nload:LR:High` = col_double(), `Nload:LR:Low` = col_double(),  `Nload:PR:High` = col_double(), `Nload:PR:Low` = col_double(), `Nload:RB:High` = col_double(), `Nload:RB:Low` = col_double(), X1 = col_skip()))
stats=gather(stats, variable,datavalue, -GAMUTload.DateTime)
stats=separate(stats,col=variable, into=c("parameter","watershed","position"))
GAMUTdischarge=subset(stats, parameter=="Q")
anova(lm(datavalue~position*watershed, stats, na.action = na.omit, subset=(parameter=="Nload")))
ggplot(means, aes(x=Watershed, y=log(Mean), fill=Site)) + geom_bar(position=position_dodge(), stat="identity") + geom_errorbar(aes(ymin=log(Mean-sd), ymax=log(Mean+sd))) +facet_grid(means$Parameter~.)

GAMUTloadzoo=read.zoo(GAMUTload,index.column="DateTime") 
GAMUTloadmonth<-apply.monthly(GAMUTloadzoo,mean,na.rm=TRUE)
write.csv(GAMUTloadmonth,"GAMUTloadmonthmeans.csv")
## Edited headings to have parameter:watershed:position format
GAMUTloadmonthmeans <- read_csv("GAMUTloadmonthmeans.csv",  col_types = cols(Date = col_date(format = "%Y/%m")))
GAMUTloadmonthmeans=gather(GAMUTloadmonthmeans, variable, datavalue, -Date)
GAMUTloadmonthmeans=separate(GAMUTloadmonthmeans,col=variable, into=c("parameter","watershed","position"))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="Nitrate")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="SpCon")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="Turb")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="pH")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="EXWTemp")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="ODO")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="Q")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="BGA")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="chla")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="fDOM")))

anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="Nload")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="Pload")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="condload")))
anova(lm(datavalue~position*watershed+Date, GAMUTloadmonthmeans, na.action = na.omit, subset=(parameter=="Cload")))

ggplot(subset(GAMUTloadmonthmeans,parameter=="Nitrate"), aes(x=Date,y=datavalue, color=position)) + geom_line() +facet_grid(watershed~.) + ylab("Nitrate (mg/L)")
ggplot(subset(GAMUTloadmonthmeans,parameter=="SpCon"), aes(x=Date,y=datavalue, color=position)) + geom_line() +facet_grid(watershed~.) + ylab("Specific Conductivity (uS/cm)")
ggplot(subset(GAMUTloadmonthmeans,parameter=="Turb"), aes(x=Date,y=datavalue, color=position)) + geom_line() +facet_grid(watershed~.) + ylab("Turbidity (NTU)")
ggplot(subset(GAMUTloadmonthmeans,parameter=="pH"), aes(x=Date,y=datavalue, color=position)) + geom_line() +facet_grid(watershed~.) + ylab("pH")
ggplot(subset(GAMUTloadmonthmeans,parameter=="EXWTemp"), aes(x=Date,y=datavalue, color=position)) + geom_line() +facet_grid(watershed~.) + ylab("Temperature (C)")
ggplot(subset(GAMUTloadmonthmeans,parameter=="ODO"), aes(x=Date,y=datavalue, color=position)) + geom_line() +facet_grid(watershed~.) + ylab("Dissolved Oxygen (mg/L)")
ggplot(subset(GAMUTloadmonthmeans,parameter=="Q"), aes(x=Date,y=datavalue, color=position)) + geom_line() +facet_grid(watershed~.) + ylab("Discharge (cms)")
ggplot(subset(GAMUTloadmonthmeans,parameter=="SpDis"), aes(x=Date,y=datavalue, color=position)) + geom_line() +facet_grid(watershed~.) + ylab("Specific Discharge (mm/Y)")

###CH
p1=ggplot(CH_fill, aes(time)) +   geom_line(aes(y=EXWTempCH)) + 
  xlab("") +ylab("Temp")  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2=ggplot(CH_fill, aes(time)) +   geom_line(aes(y=ODOCH))+ 
  xlab("") +ylab("DO")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3=ggplot(CH_fill, aes(time)) +   geom_line(aes(y=QCH))+ 
  xlab("") +ylab("Discharge")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4=ggplot(CH_fill, aes(time)) +   geom_line(aes(y=NitrateCH))+ 
  xlab("") +ylab("NO3")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p5=ggplot(CH_fill, aes(time)) +   geom_line(aes(y=SpConCH))+ 
  xlab("") +ylab("Sp. Cond.")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6=ggplot(CH_fill, aes(time)) +   geom_line(aes(y=fDOMCH))+ 
  xlab("") +ylab("fDOM")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p7=ggplot(CH_fill, aes(time)) +   geom_line(aes(y=TurbCH))+ 
  xlab("") +ylab("Turb")

png("./figures/TSCH.png",width = 8, height = 10, units = "in", res=120)
multiplot(p1,p2, p3,p4,p5,p6,p7,cols = 1)
dev.off()
###BJ
p1=ggplot(BJ_fill, aes(time)) +   geom_line(aes(y=EXWTempBJ)) + 
  xlab("") +ylab("Temp")  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2=ggplot(BJ_fill, aes(time)) +   geom_line(aes(y=ODOBJ))+ 
  xlab("") +ylab("DO")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3=ggplot(BJ_fill, aes(time)) +   geom_line(aes(y=QBJ))+ 
  xlab("") +ylab("Discharge")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4=ggplot(BJ_fill, aes(time)) +   geom_line(aes(y=NitrateBJ))+ 
  xlab("") +ylab("NO3")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p5=ggplot(BJ_fill, aes(time)) +   geom_line(aes(y=SpConBJ))+ 
  xlab("") +ylab("Sp. Cond.")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6=ggplot(BJ_fill, aes(time)) +   geom_line(aes(y=fDOMBJ))+ 
  xlab("") +ylab("fDOM")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p7=ggplot(BJ_fill, aes(time)) +   geom_line(aes(y=TurbBJ))+ 
  xlab("") +ylab("Turb")

png("./figures/TSBJ.png",width = 8, height = 10, units = "in", res=120)
multiplot(p1,p2, p3,p4,p5,p6,p7,cols = 1)
dev.off()
###ARBR
p1=ggplot(ARBR_fill, aes(time)) +   geom_line(aes(y=EXWTempARBR)) + 
  xlab("") +ylab("Temp")  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2=ggplot(ARBR_fill, aes(time)) +   geom_line(aes(y=ODOARBR))+ 
  xlab("") +ylab("DO")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3=ggplot(ARBR_fill, aes(time)) +   geom_line(aes(y=QARBR))+ 
  xlab("") +ylab("Discharge")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4=ggplot(ARBR_fill, aes(time)) +   geom_line(aes(y=NitrateARBR))+ 
  xlab("") +ylab("NO3")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p5=ggplot(ARBR_fill, aes(time)) +   geom_line(aes(y=SpConARBR))+ 
  xlab("") +ylab("Sp. Cond.")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6=ggplot(ARBR_fill, aes(time)) +   geom_line(aes(y=fDOMARBR))+ 
  xlab("") +ylab("fDOM")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p7=ggplot(ARBR_fill, aes(time)) +   geom_line(aes(y=TurbARBR))+ 
  xlab("") +ylab("Turb")

png("./figures/TSARBR.png",width = 8, height = 10, units = "in", res=120)
multiplot(p1,p2, p3,p4,p5,p6,p7,cols = 1)
dev.off()
###FD
p1=ggplot(FD_fill, aes(time)) +   geom_line(aes(y=EXWTempFD)) + 
  ylab("Temp")  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2=ggplot(FD_fill, aes(time)) +   geom_line(aes(y=ODOFD))+ 
  ylab("DO")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3=ggplot(FD_fill, aes(time)) +   geom_line(aes(y=QFD))+ 
  ylab("Discharge")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4=ggplot(FD_fill, aes(time)) +   geom_line(aes(y=NitrateFD))+ 
  ylab("NO3")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p5=ggplot(FD_fill, aes(time)) +   geom_line(aes(y=SpConFD))+ 
  ylab("Sp. Cond.")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6=ggplot(FD_fill, aes(time)) +   geom_line(aes(y=fDOMFD))+ 
  +ylab("fDOM")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p7=ggplot(FD_fill, aes(time)) +   geom_line(aes(y=TurbFD))+ 
  xlab("") +ylab("Turb")

png("./figures/TSFD.png",width = 8, height = 10, units = "in", res=120)
multiplot(p1,p2, p3,p4,p5,p6,p7,cols = 1)
dev.off()
##Merge all the sites into a single object, GAMUTload
# GAMUTload<-merge(FD,ARBR, by="DateTime")
# GAMUTload<-merge(GAMUTload,CH, by="DateTime")
# GAMUTload<-merge(GAMUTload,BJ, by="DateTime")
# GAMUTload<-merge(GAMUTload,Mendon, by="DateTime")
# GAMUTload<-merge(GAMUTload,WaterLab, by="DateTime")
# GAMUTload$time.x<-NULL
# GAMUTload$time.y<-NULL
# GAMUTload$time<-NULL
# write.csv(GAMUTload, file="GAMUTload.csv")

###Mendon
p1=ggplot(Mendon_fill, aes(time)) +   geom_line(aes(y=EXWTempMendon)) + 
  xlab("") +ylab("Temp")  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2=ggplot(Mendon_fill, aes(time)) +   geom_line(aes(y=ODOMendon))+ 
  xlab("") +ylab("DO")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3=ggplot(Mendon_fill, aes(time)) +   geom_line(aes(y=QMendon))+ 
  xlab("") +ylab("Discharge")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4=ggplot(Mendon_fill, aes(time)) +   geom_line(aes(y=NitrateMendon))+ 
  xlab("") +ylab("NO3")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p5=ggplot(Mendon_fill, aes(time)) +   geom_line(aes(y=SpConMendon))+ 
  xlab("") +ylab("Sp. Cond.")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6=ggplot(Mendon_fill, aes(time)) +   geom_line(aes(y=fDOMMendon))+ 
  xlab("") +ylab("fDOM")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p7=ggplot(Mendon_fill, aes(time)) +   geom_line(aes(y=TurbMendon))+ 
  xlab("") +ylab("Turb")

png("./figures/TSMendon.png",width = 8, height = 10, units = "in", res=120)
multiplot(p1,p2, p3,p4,p5,p6,p7,cols = 1)
dev.off()

###WaterLab
p1=ggplot(WaterLab_fill, aes(time)) +   geom_line(aes(y=EXWTempWaterLab)) + 
  xlab("") +ylab("Temp")  + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2=ggplot(WaterLab_fill, aes(time)) +   geom_line(aes(y=ODOWaterLab))+ 
  xlab("") +ylab("DO")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p3=ggplot(WaterLab_fill, aes(time)) +   geom_line(aes(y=QWaterLab))+ 
  xlab("") +ylab("Discharge")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p4=ggplot(WaterLab_fill, aes(time)) +   geom_line(aes(y=NitrateWaterLab))+ 
  xlab("") +ylab("NO3")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p5=ggplot(WaterLab_fill, aes(time)) +   geom_line(aes(y=SpConWaterLab))+ 
  xlab("") +ylab("Sp. Cond.")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p6=ggplot(WaterLab_fill, aes(time)) +   geom_line(aes(y=fDOMWaterLab))+ 
  xlab("") +ylab("fDOM")+ 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p7=ggplot(WaterLab_fill, aes(time)) +   geom_line(aes(y=TurbWaterLab))+ 
  xlab("") +ylab("Turb")

png("./figures/TSWaterLab.png",width = 8, height = 10, units = "in", res=120)
multiplot(p1,p2, p3,p4,p5,p6,p7,cols = 1)
dev.off()
