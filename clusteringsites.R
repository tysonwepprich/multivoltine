
library(plyr)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)
library(mclust)


rawcounts <- read_csv("data/rawcounts.csv")
rawcovs <- read_csv("data/rawcovariates.csv")
rawcovs <- rawcovs %>%
  mutate(year = year(SiteDate))

# get all surveys at sites/years matching sites where species has ever been seen

species <- "Peck's Skipper"
counts <- rawcounts %>% filter(CommonName == species)
allsurvs <- rawcovs %>% 
  filter(SiteID %in% unique(counts$SiteID))
allcounts <- rawcounts %>%
  filter(SiteID %in% unique(counts$SiteID))

# losing empty years at this step
test <- allcounts %>%
  group_by(SiteID, year) %>%
  complete(SeqID, subtransect) %>%
  select(SiteID, SeqID, subtransect, year) %>%
  distinct()
surveys <- merge(test, allsurvs, all.x = TRUE)

test <- merge(surveys, counts, by = c("SeqID", "SiteID", "SiteDate", "Week", "year", "subtransect"), all.x = TRUE)

test$Total <- plyr::mapvalues(test[, "Total"], from = NA, to = 0)
test$count <- plyr::mapvalues(test[, "count"], from = NA, to = 0)
finaldat <- test %>%
  group_by(SiteID, year) %>%
  mutate(uniqsubtran = length(unique(subtransect)),
         present = ifelse(count > 0, 1, 0),
         CommonName = species)


library(mclust)

clustdat <- finaldat %>% 
  ungroup() %>%
  mutate(ordinal = yday(SiteDate))
clustdat <- clustdat %>%
  filter(Total > 0) %>%
  select(SiteID, year, ordinal, Total, lat) %>%
  distinct()
clustdat <- clustdat %>% 
  # filter(year == 2012) %>%
  filter(lat < 40) %>%
  select(ordinal)
clustdat <- clustdat$ordinal

mod <- densityMclust(clustdat)
plot(mod, what = "BIC")
plot(mod, what = "density")

geo <- finaldat %>%
  ungroup() %>%
  select(lat, lon) %>%
  distinct()

sites <- read.csv("data/OHsites_reconciled.csv")
names(sites)[1] <- "site"
siteGDD <- readRDS("data/growingDD_Daymet.RDS")
siteGDD <- siteGDD %>%
  group_by(site) %>% 
  filter(yday == 365) %>%
  summarise(meanGDD = mean(cumdegday))
sites <- merge(sites, siteGDD, by = "site")
sitemod <- densityMclust(scale(sites[,c(3:5)]), G = 1:15)
plot(sitemod, what = "density")

library(ggplot2)

sites$class <- as.character(sitemod$classification)
a <- ggplot(data = sites, aes(x = lon, y = lat, group = class, color = class)) + geom_point()
a

names(sites)[1] <- "SiteID"
dat <- merge(finaldat, sites, by = c("SiteID", "Description.x", "lat", "lon"))
dat <- dat %>%
  select(year, SiteDate, SiteID, Total, class) %>%
  distinct()

outlist <- list()
for (yr in unique(dat$year)){
  for (cls in unique(dat$class)){
    temp <- dat %>% 
      filter(year == yr, class == cls) %>%
      mutate(ordinal = yday(SiteDate)) %>%
      filter(Total > 0)
      
    if(nrow(temp) > 25){
      mod <- densityMclust(temp[,c("ordinal", "Total")], G = 1:3)
      plot(mod, what = "density", main = paste(yr, cls, sep = "_"))
    }
    
  }
}



