
phiout <- get_slurm_out(sjob8)
models <- phiout[13]
# models <- slurm_out[401]

counts <- count_array[,,10] #select species here, number corresponds to row in species_list

# select sites with enough individuals counted
siteRows <- which(rowSums(counts, na.rm = TRUE) >= 10)
counts <- counts[siteRows, ]
counts[is.na(counts)] <- -1


# Covariate (latitude)
lat <- cov_sites[siteRows, "lat"]
# lat <- cov_sites[siteRows, "YearGDD"]
cov.w <- cov.mu <- scale(lat)


temp <- models[[1]]

#---  for selected model ---#
FittedVal <- FittedVal.f(temp)



#Tyson's plots
S <- dim(counts)[1]
t <- dim(counts)[2]
df_all <- data.frame()
for (i in 1:S){
  df <- data.frame(Week = 1:t, Site = rep(i, t))
  df$Fit <- FittedVal[i,]
  df$Count <- counts[i,]
  df$Count[df$Count == -1] <- NA
  df$Phen <- temp$betta.est[i,]
  df$Latitude <- -cov.mu[i]
  
  df_all <- rbind(df_all, df)
}

library(dplyr)
df_all <- df_all %>% group_by(Site) %>% mutate(Tot = sum(Count, na.rm = TRUE)) %>% filter(Tot > 10)

library(ggplot2)

#plot fitted phenology by site, colored by latitude
c <- ggplot(data = df_all, aes(x = Week, y = Phen, group = Site)) + geom_line(aes(color = Latitude), size = .8) 
c + theme_bw() +
  scale_colour_gradient2(name = "Latitude", midpoint = mean(range(df_all$Latitude)), low = "red", mid = "yellow", high = "blue") 

#plot fitted values and actual counts for all sites
#Scaled latitude is multiplied by -1, so Southern sites at bottom of the plot
d <- ggplot(data = df_all, aes(x = Week, y = Count, group = Site)) + geom_point() + 
  facet_wrap(Latitude ~ Site, ncol = 4, scales = "free_y") + geom_line(aes(x = Week, y = Fit)) + theme_bw()
d


# weird estimates of 3rd generation outside of monitoring weeks, have small weights
# why even choose 3 generations over 2 if one of them will have no weight?
hist(temp$mu.est, breaks = 20)
t <- temp$mu.est[temp$w.est > 0.01]
hist(t, breaks = 20)
