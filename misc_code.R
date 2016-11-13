# for gamlist, half size of file by removing redundant gam model

fs <- list.files('gamOrdinal')
for (i in 1:length(fs)){
  gamlist <- readRDS(paste("gamOrdinal", fs[i], sep = "/"))
  if("gam" %in% class(gamlist$mod)){
    if("gam" %in% class(gamlist$modgdd)){
      gamlist$modgdd <- NULL
      saveRDS(gamlist, file = paste("gamOrdinal", fs[i], sep = "/"))
    }
  }
}

# looking at raw counts to see why some years suck at separating N and P


ReduceCounts <- function(list_index){
  counts <- dat[[list_index]]$counts
  raw_cutoff <- 10
  siteRows <- which(rowSums(counts, na.rm = TRUE) >= raw_cutoff)
  surv_present <- which(apply(counts, 1, function(x) length(which(x > 0))) >= 3)
  siteRows <- siteRows[which(siteRows %in% surv_present)]
  return(counts[siteRows, ])
}


ct <- ReduceCounts(3)

matplot(t(ct), type = "l")


which(paramIN$nRun == 35)



actual <- rowSums(ct, na.rm = TRUE)
fitted <- outList[[which(paramIN$nRun == 53)]][[1]]$N.est
fittedbad <- outList[[which(paramIN$nRun == 65)]][[1]]$N.est
fittedcom <- outList[[which(paramIN$nRun == 41)]][[1]]$N.est


plot(actual, fittedcom)



# List lengths vs abundance
dat[, listlength := length(CommonName), by = SeqID]
dat[, totalbfly := sum(Total), by = SeqID]




plot(fitted, fittedbad)

