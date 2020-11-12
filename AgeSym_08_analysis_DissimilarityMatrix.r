###########################################################
## Purpose: compute dissimilarity matrix for PAM clustering        
###########################################################

args <- commandArgs(TRUE)
cohort <- as.character(args[1])
base <- "/Users/jamesroe/Dropbox/OpenScienceFramework/AgeSym"
resdir <- file.path(base,"results/reproduceClustering",cohort)


# open libraries
#..............#
packages <- c("dplyr", "stringr", "magrittr")
lapply(packages, require, character.only = T)


# load trajectories within signficant clusters
#............................................#
fit_val = t(read.csv(file.path(resdir,"fitFitfsav5.csv"),stringsAsFactors = F, sep = ",", header = F))
end = dim(fit_val)[2]
mat.dist.fit = matrix(rep(0,end*end), nrow=end)


#compute dissimilaritiy matrix
# ...........................#
pb = txtProgressBar(min=1, max=end, style=3)
for(k in 1:end) {
  setTxtProgressBar(pb,k)
  for(j in 1:end){
  	mat.dist.fit[k,j] <- sum( (fit_val[,k] - fit_val[,j]) ^2 ) #dissimilarity matrix based on sum of squares between difference trajectories
  }
}

save("mat.dist.fit", "fit_val",
     file = file.path(resdir, "DS_summary.Rda"))
write.table(mat.dist.fit,file.path(resdir, "distancemat.csv"),quote=F,row.names = F, col.names = F, sep=" ")
quit()

