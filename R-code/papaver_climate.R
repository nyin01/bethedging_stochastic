clim <- read.csv("england_wintertemp.csv") #load csv of mean annual winter temps
clim <- subset(clim, Year>1659)
clim <- subset(clim, Year<1965)
colnames(clim) <- c("year", "mean_temp")

clim_day <- read.csv("winter_daytemp_england.csv") #load csv of daily winter temps
#re-format date column
clim_day$year <- substr(clim_day$Date, 1, 4) #separate year
clim_day$month <- substr(clim_day$Date, 6, 7) #separate month
clim_day <- clim_day %>% select(c("Date", "year", "month", "Value"))
clim_day$year <- as.numeric(clim_day$year)
clim_day$year[clim_day$month %in% c("11", "12")] <- clim_day$year[clim_day$month %in% c("11", "12")] + 1 #reassign year to "winter year" for november and december

clim_day_summary_0 <- clim_day %>% #calculate number of days below 0C
  group_by(year) %>%
  tally(Value < 0)
colnames(clim_day_summary_0) <- c("year", "below.0")

clim_day_summary_1 <- clim_day %>% #calculate number of days below -1.6C
  group_by(year) %>%
  tally(Value < -1.6)
colnames(clim_day_summary_1) <- c("year", "below.deepfreeze")

clim_day_summary <- cbind(clim_day_summary_0, below.deepfreeze = clim_day_summary_1$below.deepfreeze)
clim_day_summary$above_zero <- 89-clim_day_summary$below.0 #convert from count below 0 to count above 0
clim_day_summary$above_deepfreeze <- 89-clim_day_summary$below.deepfreeze #convert from count below -1.6 to count above -1.6

arthur_thresholds <- data.frame(year = c(1966, 1976, 1968), 
                                winter_type = c("mild", "mild", "harsh"),
                                mean_temp = c(4.4, 5.1, 3.6),
                                below.0 = c(7, 2, 11),
                                below.deepfreeze = c(3, 1, 6))

mean_temp_thresh = mean(arthur_thresholds$mean_temp[arthur_thresholds$year == 1968], arthur_thresholds$mean_temp[arthur_thresholds$year == 1966])
below_0_thresh = mean(arthur_thresholds$below.0[arthur_thresholds$year == 1968], arthur_thresholds$below.0[arthur_thresholds$year == 1966])
below_deepfreeze_thresh = mean(arthur_thresholds$below.deepfreeze[arthur_thresholds$year == 1968], arthur_thresholds$below.deepfreeze[arthur_thresholds$year == 1966])

prop_mild_meantemp = nrow(clim$mean_temp > mean_temp_thresh)/nrow(clim)
prop_mild_0 = nrow(clim_day_summary$below.0 < below_0_thresh)/nrow(clim_day_summary)
prop_mild_deepfreeze = nrow(clim_day_summary$below.deepfreeze < below_deepfreeze_thresh)/nrow(clim_day_summary)

ggplot(data=clim_day_summary) +
  geom_histogram(aes(x=below.deepfreeze)) +
  geom_vline(aes(xintercept=mean(below.deepfreeze)),
             color="black", linetype="solid", size=1) +
  geom_vline(aes(xintercept= arthur_thresholds$below.deepfreeze[arthur_thresholds$year == 1968]), color="lightblue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= arthur_thresholds$below.deepfreeze[arthur_thresholds$year == 1967]), color="firebrick2", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= arthur_thresholds$below.deepfreeze[arthur_thresholds$year == 1966]), color="firebrick", linetype="dashed", size=1) +
  theme_bw() +
  scale_x_reverse() +
  xlab("Number of days below -1.6C") +
  ylab("Frequency") +
  theme(text = element_text(size = 20)) +
  theme(text = element_text(family = "Times"))

ggplot(data=clim_day_summary) +
  geom_histogram(aes(x=below.0)) +
  geom_vline(aes(xintercept=mean(below.0)),
             color="black", linetype="solid", size=1) +
  geom_vline(aes(xintercept= arthur_thresholds$below.0[arthur_thresholds$year == 1968]), color="lightblue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= arthur_thresholds$below.0[arthur_thresholds$year == 1967]), color="firebrick2", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= arthur_thresholds$below.0[arthur_thresholds$year == 1966]), color="firebrick", linetype="dashed", size=1) +
  theme_bw() +
  scale_x_reverse() +
  xlab("Number of days below 0C") +
  ylab("Frequency") +
  theme(text = element_text(size = 20)) +
  theme(text = element_text(family = "Times"))

ggplot(data=clim) +
  geom_histogram(aes(x=mean_temp)) +
  geom_vline(aes(xintercept=mean(mean_temp)),
             color="black", linetype="solid", size=1) +
  geom_vline(aes(xintercept= arthur_thresholds$mean_temp[arthur_thresholds$year == 1968]), color="lightblue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= arthur_thresholds$mean_temp[arthur_thresholds$year == 1967]), color="firebrick2", linetype="dashed", size=1) +
  geom_vline(aes(xintercept= arthur_thresholds$mean_temp[arthur_thresholds$year == 1966]), color="firebrick", linetype="dashed", size=1) +
  theme_bw() +
  xlab("Mean Winter Temperature") +
  ylab("Frequency") +
  theme(text = element_text(size = 20)) +
  theme(text = element_text(family = "Times"))



