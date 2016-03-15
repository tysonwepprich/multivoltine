#miami blue phenology
library(lubridate)
library(mgcv)
library(ggplot2)

rawdata <- read.csv("data/KWNWR counts.csv", header = TRUE)
dat <- rawdata[complete.cases(rawdata), ] # remove blank lines

dat$date <- mdy(dat$date)
dat$day <- as.numeric(difftime(dat$date, min(dat$date), units = "days")) + 1
dat$week <- round(as.numeric(difftime(dat$date, min(dat$date), units = "weeks")) + 1)

# choose whether to add zero counts before and after to ground the GAM estimates
# add zero counts a month outside of monitoring season
zeroCounts <- expand.grid(unique(dat$island), c(-31,-30, max(dat$day) + 30, max(dat$day) + 31))
names(zeroCounts) <- c("island", "day")
zeroCounts$Count <- 0
dat_zeros <- rbind(dat[, c("island", "day", "Count")], zeroCounts)




# choose model
# use dat or dat_zeros
# also, adjust k (knots) to determine how flexible lines will be. 
# big difference between 10 or 20

# each island has same phenology, different population size
mod <- gam(Count ~ island + s(day, k = 20), 
           family = poisson(link = "log"), data = dat, 
           method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))

# each island has different phenology and population size
# k > 15 gave errors for too many parameters to estimate
mod <- gam(Count ~ island + s(day, k = 10, by = island), 
           family = poisson(link = "log"), data = dat, 
           method = "GCV.Cp", optimizer = "perf", gamma = 1.4, control = list(maxit = 500))


# get predicted values for each day and plot the abundance over time
start.surv <- min(dat$day)
end.surv <- max(dat$day)

newData <- expand.grid(unique(dat$island), c(start.surv:end.surv))
names(newData) <- c("island", "day")

# here's where daily estimated abundance for each island will be
newData$GAM.pred <- predict.gam(mod, newData, type = "response")

# plot estimated abundance over time from GAM curves
a <- ggplot(data = newData, aes(x = day, y = GAM.pred, group = island, color = island)) + geom_line()
a


# look at matrix of counts x weeks to see if Matechou model might be used
# doubtful with the # of NA's
# one option would be averaging counts at island over 2 weeks and use 2 weeks as time period, 
# but not sure if that'd work well
library(reshape)
library(dplyr)

# see if more than one count per week/island
# only happens once, take mean of 2 counts
counts_uniq <- dat %>% 
  group_by(island, week) %>%
  summarise(weekcount = mean(Count))

# weeks with no observations missing from this
count_matrix <- as.matrix(cast(counts_uniq, island ~ week, value = "weekcount"))
count_matrix <- round(count_matrix)
