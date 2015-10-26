# finding results when you clear the workspace and get_slurm_out no longer works

slurm_out <- list()
missing_files <- c()
tmpEnv <- new.env()
for (i in 0:3) {
  fname <- paste0("slr8743", "_", i, 
                  ".RData")
  if (fname %in% dir()) {
    load(fname, envir = tmpEnv)
    slurm_out <- c(slurm_out, get(".rslurm_result", 
                                  envir = tmpEnv))
  }
  else {
    missing_files <- c(missing_files, fname)
  }
}


outList <- slurm_out
outDF <- list()
for (i in 1:length(outList)){
  out <- outList[[i]]$pars
  out$bootstrap <- outList[[i]]$bootstrap
  out$ll.val <- outList[[i]]$ll.val
  out$time <- as.double(outList[[i]]$time, units = "mins")
  outDF[[i]] <- out
}

outDF <- do.call("rbind", outDF)
