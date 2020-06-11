#### Calculate GAMUT station metabolism
#### basic protocol from streamMetabolizer
#### Written by Dylan Dastrup and Erin Jones erinfjones3@gmail.com
#### Last updated 04 May 2019


#### Set Up ####
# install.packages("streamMetabolizer", dependencies=TRUE, 
#                  repos=c("https://owi.usgs.gov/R","https://cran.rstudio.com"))

library("dplyr")
require(cowplot)
library(streamMetabolizer)
library(readr)
library(zoo)
theme_set(theme_bw())                              ## I also like theme_calc, theme_classic, theme_gray
Sys.timezone()
setwd("~/Box Sync/Aanderud Lab/Projects/iUTAH/PR_CH Metabolism 2018/GAMUT Paper 2018/GAMUT 2018 streamMetabolizer_R/Sensor_downloads/metabolism") ##BYU
setwd("C:/Users/Erin/Box Sync/Aanderud Lab/Projects/iUTAH/PR_CH Metabolism 2018/GAMUT Paper 2018/GAMUT 2018 streamMetabolizer_R/Sensor_downloads/metabolism") ##Lappy
funfill<- function(x) na.approx(x, maxgap = 8, na.rm = FALSE)
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

#### Read in data (made using GAMUT_GPPdownloads.R)###

gamutCH <- read_csv("CHmetab190815_clean.csv",  col_types = cols(X1 = col_skip(), time= col_datetime(format = "%m/%d/%y %H:%M" ), EXWTempCH = col_number(), ODOCH = col_number(),  
                                              ODO_SatCH = col_number(), PARIn_AvgCH = col_number(), DepthCH= col_number() ))
gamutMendon = read_csv("Mendonmetab190815_clean.csv", col_types = cols(X1 = col_skip(), time= col_datetime(format = "%m/%d/%y %H:%M"), EXWTempMendon = col_number(), ODOMendon = col_number(), 
                                                ODO_SatMendon = col_number(), PARIn_AvgMendon = col_number(), DepthMendon= col_number() ))
gamutFD = read_csv("FDmetab190815_clean.csv", col_types = cols(X1 = col_skip(), time= col_datetime(format = "%m/%d/%y %H:%M"), EXWTempFD= col_number(), ODOFD = col_number(),
                                                ODO_SatFD = col_number(), PARIn_AvgFD = col_number(), DepthFD= col_number() ))

### Fill in gaps (up to 2 hours)

gamutCH_fill=gamutCH %>%
  mutate_each (funs(funfill), EXWTempCH, ODOCH, ODO_SatCH, PARIn_AvgCH, DepthCH)
gamutMendon_fill=gamutMendon %>%
  mutate_each (funs(funfill), EXWTempMendon, ODOMendon, ODO_SatMendon, PARIn_AvgMendon, DepthMendon)
gamutFD_fill=gamutFD %>%
  mutate_each (funs(funfill), EXWTempFD, ODOFD, ODO_SatFD, PARIn_AvgFD, DepthFD)

p1=ggplot(gamutCH_fill, aes(time)) + geom_line(aes(y=EXWTempCH))
p4=ggplot(gamutCH_fill, aes(time)) + geom_line(aes(y=ODOCH))
p2=ggplot(gamutCH_fill, aes(time)) + geom_line(aes(y=ODO_SatCH))
p3=ggplot(gamutCH_fill, aes(time)) + geom_line(aes(y=PARIn_AvgCH))
p5=ggplot(gamutCH_fill, aes(time)) + geom_line(aes(y=DepthCH))
multiplot( p1,p4,p2,p3,p5,cols = 1)


p1=ggplot(gamutMendon_fill, aes(time)) + geom_line(aes(y=EXWTempMendon))
p4=ggplot(gamutMendon_fill, aes(time)) + geom_line(aes(y=ODOMendon))
p2=ggplot(gamutMendon_fill, aes(time)) + geom_line(aes(y=ODO_SatMendon))
p3=ggplot(gamutMendon_fill, aes(time)) + geom_line(aes(y=PARIn_AvgMendon))
p5=ggplot(gamutMendon_fill, aes(time)) + geom_line(aes(y=DepthMendon))
multiplot( p1,p4,p2,p3,p5,cols = 1)

p1=ggplot(gamutFD_fill, aes(time)) + geom_line(aes(y=EXWTempFD))
p4=ggplot(gamutFD_fill, aes(time)) + geom_line(aes(y=ODOFD))
p2=ggplot(gamutFD_fill, aes(time)) + geom_line(aes(y=ODO_SatFD))
p3=ggplot(gamutFD_fill, aes(time)) + geom_line(aes(y=PARIn_AvgFD))
p5=ggplot(gamutFD_fill, aes(time)) + geom_line(aes(y=DepthFD))
multiplot( p1,p4,p2,p3,p5,cols = 1)

gamutCH_fill$light=convert_PAR_to_SW(gamutCH_fill$PARIn_AvgCH, coef = 0.473)
ggplot(gamutCH_fill, aes(gamutCH.time)) + geom_line(aes(y=light))

gamutMendon_fill$Hydroperiod=getHydroperiod(gamutMendon_fill$time)
gamutCH_fill$Hydroperiod=getHydroperiod(gamutCH_fill$time)
gamutFD_fill$Hydroperiod=getHydroperiod(gamutFD_fill$time)

p = ggpairs(data = gamutFD_fill, columns = 2:6, aes(color = Hydroperiod, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
ggsave("figures/corFD.png", plot = p, width = 10, height = 10)
p = ggpairs(data = gamutMendon_fill, columns = 2:6, aes(color = Hydroperiod, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
ggsave("figures/corMendon.png", plot = p, width = 10, height = 10)
p = ggpairs(data = gamutCH_fill, columns = 2:6, aes(color = Hydroperiod, alpha = 0.7), 
            upper = list(continuous = wrap("cor", method = "spearman"))) + theme_few()
ggsave("figures/corCH.png", plot = p, width = 10, height = 10)

########### Calculate Charleston metabolism ##############
##In UTC
posix.time.localtz= as.POSIXct(gamutCH_fill$time, "%m/%d/%y %H:%M",tz= "UTC")
posix.time.solar <- streamMetabolizer::convert_UTC_to_solartime(posix.time.localtz, longitude=-111.47, time.type = "mean solar") #Convert DateTime to Solar Time
datCH <- transmute(gamutCH_fill,
                     solar.time = posix.time.solar,
                     DO.obs = gamutCH_fill$ODOCH, 
                     DO.sat = gamutCH_fill$ODOCH/(gamutCH_fill$ODO_SatCH/100),
                     depth = gamutCH_fill$DepthCH,
                     temp.water = gamutCH_fill$EXWTempCH,
                     light = convert_PAR_to_SW(gamutCH_fill$PARIn_AvgCH, coef = 0.473) )

bayes_name <- mm_name(type='bayes', pool_K600='normal', err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_name

bayes_specs <- specs(bayes_name)
bayes_specs

mmCH <- metab(bayes_specs, data=datCH)
mmM <- metab(bayes_specs, data=datM)
mmCH


predict_DataCH <- predict_metab(mmCH)
write.csv(predict_DataCH, file = "PR_CH_predict_Data.csv" )


plot_metab_preds(mmCH)
ggsave("figures/CH_Metab_preds.png",width=8, height=4)


param_DataCH <- get_params(mmCH)
write.csv(param_DataCH, file = "PR_CH_param_Data.csv" )

predict_DO(mmCH) %>% head()


plot_DO_preds(mmCH)
ggsave("figures/CH_DO_preds.png",width=40, height=15)

mcmcCH <- get_mcmc(mmCH)
rstan::traceplot(mcmcCH, pars='K600_daily', nrow=10)

get_fit(mmCH)$overall %>%
  select(ends_with('Rhat'))
get_fit(mmCH) %>%
  lapply(names)


######### Calculate Mendon Rd metabolism ###########
### In UTC timezone
posix.time.localtz= as.POSIXct(gamutMendon_fill$time, tz= "UTC","%m/%d/%y %H:%M")
posix.time.solar <- streamMetabolizer::calc_solar_time(posix.time.localtz, longitude=-111.89) #Convert DateTime to Solar Time
datM <- transmute(gamutMendon_fill,
                    solar.time = posix.time.solar,
                    DO.obs = gamutMendon_fill$ODOMendon, 
                    DO.sat = gamutMendon_fill$ODOMendon/(gamutMendon_fill$ODO_SatMendon/100),
                    depth = gamutMendon_fill$DepthMendon,
                    temp.water = gamutMendon_fill$EXWTempMendon,
                    light = convert_PAR_to_SW(gamutMendon_fill$PARIn_AvgMendon, coef = 0.473) )



bayes_name <- mm_name(type='bayes', pool_K600='binned', err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_name

bayes_specs <- specs(bayes_name)
bayes_specs

mmM <- metab(bayes_specs, data=datM)
mmM


predict_DataM <- predict_metab(mmM)
write.csv(predict_DataM, file = "LR_M_predict_Data.csv" )


plot_metab_preds(mmM)
ggsave("figures/MR_Metab_preds.png",width=8, height=4)


param_DataM <- get_params(mmM)
write.csv(param_DataM, file = "LR_M_param_Data.csv" )

predict_DO(mmM) %>% head()


plot_DO_preds(mmM)
ggsave("figures/MR_DO_preds.png",width=40, height=15)

# mcmcM <- get_mcmc(mmM)
# rstan::traceplot(mcmcM, pars='K600_daily', nrow=3)
# 
# get_fit(mmM)$overall %>%
#   select(ends_with('Rhat'))
# get_fit(mmM) %>%
#   lapply(names)


############ Calculate Foothill Dr Metabolism #############

posix.time.localtz= as.POSIXct(gamutFD_fill$time, tz= "UTC","%m/%d/%y %H:%M")
posix.time.solar <- streamMetabolizer::calc_solar_time(posix.time.localtz, longitude=-111.83) #Convert DateTime to Solar Time
datFD <- transmute(gamutFD_fill,
                    solar.time = posix.time.solar,
                    DO.obs = gamutFD_fill$ODOFD, 
                    DO.sat = gamutFD_fill$ODOFD/(gamutFD_fill$ODO_SatFD/100),
                    depth = gamutFD_fill$DepthFD,
                    temp.water = gamutFD_fill$EXWTempFD,
                    light = convert_PAR_to_SW(gamutFD_fill$PARIn_AvgFD, coef = 0.473) )

Q_daily<-read.csv("Q_daily.csv", stringsAsFactors = FALSE)
Q_daily$date <- as.Date(Q_daily$date,format="%m/%d/%y")

bayes_name <- mm_name(type='bayes', pool_K600='binned', err_obs_iid=TRUE, err_proc_iid=TRUE)
bayes_name

bayes_specs <- specs(bayes_name,K600_lnQ_nodes_meanlog=log(15.5), K600_lnQ_nodes_centers=c(log(0.01),log(0.02),log(0.03),log(0.04),log(0.05),log(0.06),log(0.07),log(0.08),log(0.09),log(0.1),log(0.15),log(0.25)))
bayes_specs

mmFD <- metab(bayes_specs, data=datFD, data_daily=Q_daily)
mmFD


predict_DataFD <- predict_metab(mmFD)
predict_DataFD$GPP_ER = predict_DataFD$GPP/predict_DataFD$ER
write.csv(predict_DataFD, file = "RB_FD_predict_Data191115.csv" )


plot_metab_preds(mmFD)
ggsave("figures/FD_Metab_preds191115.png",width=8, height=4)


param_DataFD <- get_params(mmFD)
write.csv(param_DataFD, file = "RB_FD_param_Data191115.csv" )

predict_DO(mmFD) %>% head()


plot_DO_preds(mmFD)
ggsave("figures/FD_DO_preds191115.png",width=12, height=6)

mcmcFD <- get_mcmc(mmFD)
rstan::traceplot(mcmcFD, pars='K600_daily', nrow=6)

get_fit(mmFD)$overall %>%
  select(ends_with('Rhat'))
get_fit(mmFD) %>%
  lapply(names)


# 
# ###### Write output to avoid having to rerun again 
saveRDS(mmCH, file="mmCH190815.rda")
saveRDS(mmM, file="mmM190815.rda")
saveRDS(mmFD, file="mmFD191115.rda")
# 
# ##### Read in existing metabolism objects
 mmCH<-readRDS(file="mmCH190610.rda")
 mmM<-readRDS(file="mmM190610.rda")
 mmFD<-readRDS(file="mmFD190610.rda")

# plot_metab_preds(mmCH)
# plot_metab_preds(mmM)
