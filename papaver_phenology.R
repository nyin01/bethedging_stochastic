fall_wt_df <- read.csv("/Users/mayaweissman/Documents/GitHub/SignInversionNoise/julia/817papaver_fall_wt_e1.p1.csv")
fall_wt_df$bh.gmf <- 10.01^(fall_wt_df$pA)*0.735^(1-fall_wt_df$pA)
fall_wt_df$wt.gmf <- 18^(fall_wt_df$pA)*0.5^(1-fall_wt_df$pA)
fall_wt_df$dGMF <- fall_wt_df$bh.gmf - fall_wt_df$wt.gmf
fall_wt_df$exp_sign[fall_wt_df$dGMF < 0] <- "Deleterious"
fall_wt_df$exp_sign[fall_wt_df$dGMF > 0] <- "Beneficial"
fall_wt_df$NPfix[fall_wt_df$NPfix == 0] <- 10^(-6)

npf <- ggplot(data=fall_wt_df, aes(x=N, y=NPfix, color=as.factor(pA), shape = exp_sign)) +
  geom_rect(aes(xmin = 500, xmax = 5000, ymin = 10^(-4), ymax = 10^4), fill = "lightgrey", color = "lightgrey") +
  geom_line(size=2) +
  geom_point(size=5, alpha = 0.9) +
  geom_hline(yintercept=1, color="gray25", size=1.5, linetype="longdash") + 
  scale_color_viridis(discrete = TRUE, direction = -1, option="D", guide = guide_legend(reverse = FALSE), name = expression(P[Mild~Winter])) + 
  scale_x_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+     
  scale_y_continuous(trans = log10_trans(), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  labs(shape="Expected Sign") +
  xlab("Population Size (N)") + 
  ylab(expression(NP[fix])) +
  coord_cartesian(ylim= c(10^(-2), 10^(3))) +
  theme_bw()+ 
  theme(text = element_text(family = "Times")) +
  theme(text = element_text(size = 20)); npf

ggsave("912papavernpf.png", npf)

pa <- c(0,
        0.05,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1)
dgmf <- c(0.235,
          0.2394034155,
          0.238853011,
          0.2152932171,
          0.1438305279,
          -0.007444825872,
          -0.2875564522,
          -0.7710205301,
          -1.57013489,
          -2.852953018,
          -4.869509011,
          -7.99)
gmf_df <- cbind(pa,dgmf)
gmf_df <- as.data.frame(gmf_df)

dgmfp <- ggplot() +
  geom_rect(aes(xmin = 0.428, xmax = 0.593, ymin = -8, ymax = 1), fill = "lightgrey", color = "lightgrey", alpha = 0.5) +
  geom_line(data = gmf_df, aes(x = pa, y = dgmf), color = "blue", size = 2) +
  geom_hline(aes(yintercept = 0), color = "black", linetype ="dashed", size = 1) +
  geom_vline(aes(xintercept=0.396), color="black", linetype="dashed", size=1) +
  coord_cartesian(ylim= c(-2, 0.3)) +
  theme_bw()+ 
  xlab("Probability of a mild winter") +
  ylab("ΔGMF") +
  theme(text = element_text(family = "Times")) +
  theme(text = element_text(size = 20)); dgmfp

ggsave("912papaver_dgmf.png", dgmfp)

#salmonella
pa <- c(0,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.75,
        0.8,
        0.85,
        0.9,
        0.925,
        0.93,
        0.931,
        0.94,
        0.95,
        1)
dgmf <- c(0.00703,
          0.01101807636,
          0.01704151329,
          0.02590042318,
          0.03843060659,
          0.05508217748,
          0.07481893975,
          0.09252914549,
          0.09660950262,
          0.09327948176,
          0.07714584032,
          0.04002699912,
          0.009976540557,
          0.00280922871,
          0.001325178564,
          -0.01282378123,
          -0.03031746599,
          -0.152)
gmf_df <- cbind(pa,dgmf)
gmf_df <- as.data.frame(gmf_df)

dgmfp <- ggplot() +
  geom_line(data = gmf_df, aes(x = pa, y = dgmf), color = "blue", size = 2) +
  geom_hline(aes(yintercept = 0), color = "black", linetype ="dashed", size = 1) +
  geom_vline(aes(xintercept=0.9315), color="black", linetype="dashed", size=1) +
  #coord_cartesian(ylim= c(-2, 0.3)) +
  theme_bw()+ 
  xlab("Probability of no antibiotics") +
  ylab("ΔGMF") +
  theme(text = element_text(family = "Times")) +
  theme(text = element_text(size = 20)); dgmfp
ggsave("912salmonella_dgmf.png", dgmfp)


clim <- read.csv("/Users/mayaweissman/Documents/england_wintertemp.csv")
clim <- subset(clim, Year>1659)
clim <- subset(clim, Year<1965)

nrow(clim[clim$Win>4,])/nrow(clim)

clim_day <- read.csv("/Users/mayaweissman/Documents/winter_daytemp_england.csv")
clim_day$year <- substr(clim_day$Date, 1, 4)
clim_day$month <- substr(clim_day$Date, 6, 7)
clim_day <- clim_day %>% select(c("Date", "year", "month", "Value"))
clim_day$year <- as.numeric(clim_day$year)
clim_day$year[clim_day$month %in% c("11", "12")] <- clim_day$year[clim_day$month %in% c("11", "12")] + 1

clim_day_summary_0 <- clim_day %>%
  group_by(year) %>%
  tally(Value < 0)
colnames(clim_day_summary_0) <- c("year", "below.0")

clim_day_summary_1 <- clim_day %>%
  group_by(year) %>%
  tally(Value < -1.6)
colnames(clim_day_summary_1) <- c("year", "below.deepfreeze")

clim_day_summary <- cbind(clim_day_summary_0, below.deepfreeze = clim_day_summary_1$below.deepfreeze)
clim_day_summary$above_zero <- 89-clim_day_summary$below.0
clim_day_summary$above_deepfreeze <- 89-clim_day_summary$below.deepfreeze

ggplot(data=clim_day_summary) +
  geom_histogram(aes(x=below.deepfreeze)) +
  geom_vline(aes(xintercept=mean(below.deepfreeze)),
             color="black", linetype="solid", size=1) +
  geom_vline(aes(xintercept=6), color="lightblue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=1), color="firebrick", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=3), color="firebrick2", linetype="dashed", size=1) +
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
  geom_vline(aes(xintercept=11), color="lightblue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=2), color="firebrick", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=7), color="firebrick2", linetype="dashed", size=1) +
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
  geom_vline(aes(xintercept=3.6), color="lightblue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=5.1), color="firebrick", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=4.4), color="firebrick2", linetype="dashed", size=1) +
  theme_bw() +
  xlab("Mean Winter Temperature") +
  ylab("Frequency") +
  theme(text = element_text(size = 20)) +
  theme(text = element_text(family = "Times"))



