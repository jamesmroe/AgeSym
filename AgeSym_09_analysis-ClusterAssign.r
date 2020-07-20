###################################
## Purpose: PAM clustering protocol        
###################################

args <- commandArgs(TRUE)
outdir <- as.character(args[1])


# load libraries
# .............#
packages <- c("dplyr", "stringr", "magrittr", "cluster","ggplot2")
sapply(packages, require, character.only = T)


# load data
# ........#
vtx <- readLines(file.path(outdir,"nvtxfsav5.csv")) #clustering constrained within significant vertices on fsaverage5
load(file = file.path(outdir, "Opt.Rda")) #100 obs age model
load(file = file.path(outdir, "DS_summary.Rda")) #dissimilarity matrix
hemieffect = t(read.csv(file.path(outdir, "mapHCoef.fsav5.csv"), header = F)) #hemisphere effect on fsaverage5
age <- Opt$Age
if (!dim(mat.dist.fit)[2] == vtx) { print("matrix does not match with n of vertices"); quit() }


# compute silhouettes #
# ...................#
cl.sil = 2:7
sil.coef.dist = NULL
for (i in cl.sil) {
  print(i)
  
  # PAM clustering (2-7 solutions)
  cl.dist = pam(mat.dist.fit, k=i, diss=T) 

  # silhouette
  sil.dist = silhouette(cl.dist) 

  # extract silhouette coefficients 
  sil.coef.dist[i] = summary(sil.dist)$avg.width
}

sols = data.frame("cl"=c(2:max(cl.sil)),
                  "dist"=sil.coef.dist[-1])


#silhouette coefficients
ggplot(sols) +
  geom_vline(xintercept = 3, linetype=2,col="dark grey") +
  geom_point(aes(x=cl,y=dist),size=2,col="#00886e") +
  geom_line(aes(x=cl,y=dist),size=1,col="#00886e") +
  ylim(c(0,0.8)) +
  xlim(c(2,max(cl.sil))) +
  theme_classic() +
  labs(y="Silhouette coefficient",
       x="N clusters") +
  ggtitle("LCBC") +
  theme(text = element_text(size=12,family="Nimbus Sans Narrow"),
        axis.title.x = element_text(vjust=1),
        plot.title = element_text(hjust=0.5))



## plot cluster solutions ##
# .........................#
LM.dist = 3 #3 cluster solution
saveord=0
col.cl.main = c('lightblue1','orange1','gold1')
col.cl.mean = c('lightblue4','orange4','gold4')
ymin = min(fit_val[,])- 0.01
ymax = max(fit_val[,])+ 0.01


for (jj in 1:length(LM.dist)) {
  
  ncl <- LM.dist[jj]
  print(paste(ncl,"cluster solution"))
  cluster = pam(mat.dist.fit, k=ncl, diss=T)
  cl.ord = cluster$clustering

  for (cl in 1:ncl) {
  
    #vertex plot
    plot(0, type='n', xlab='Age', ylab='Thickness asymmetry (LH-RH)', 
         main=paste('N =', toString(length(which(cl.ord == cl)))), 
         xlim=c(min(age),max(age)), ylim=c(ymin,ymax))
    
    #plot vertex trajectories with added Hemisphere effect
    for (n in which(cl.ord == cl)) {
      lines(age,fit_val[,n] + hemieffect[n], 
            type = 'l', lty = 1, lwd = 2, col = col.cl.main[cl])
    }
    #mean trajectory with added Hemisphere effect
    lines(age,(rowMeans(fit_val[,which(cl.ord == cl)])) + mean(hemieffect[,which(cl.ord == cl)]), 
          type = 'l', lty = 1, lwd = 5, col = col.cl.mean[cl])
    #symmetry line
    lines(age, rep(0,length(age)), 
          type = "l", lty = 2, lwd = 1, col="black")
  
  }
    
  # write out cluster affiliations for surface plotting
  # if (saveord == 1) {
    # write(cl.ord, file = file.path(plotdir, paste('order_clusters', toString(ncl),'.txt', sep = "")), ncolumns = 1)
  # }
}



# Plot variation
# .............#
#cluster1
meancurve1 = rowMeans(fit_val[,which(cl.ord == 1)])+mean(hemieffect[,which(cl.ord == 1)])
sdd1 = apply(fit_val[,which(cl.ord == 1)],1,sd)
res1 = data.frame(meancurve1, sdd1, age)
#cluster2
meancurve2 = rowMeans(fit_val[,which(cl.ord == 2)])+mean(hemieffect[,which(cl.ord == 2)])
sdd2 = apply(fit_val[,which(cl.ord == 2)],1,sd)
res2 = data.frame(meancurve2, sdd2, age)
#cluster3
meancurve3 = rowMeans(fit_val[,which(cl.ord == 3)])+mean(hemieffect[,which(cl.ord == 3)])
sdd3 = apply(fit_val[,which(cl.ord == 3)],1,sd)
res3 = data.frame(meancurve3, sdd3, age)

  
#save curves for cross-cohort comparison
curves=list()
curves[[1]] = meancurve1
curves[[2]] = meancurve2
curves[[3]] = meancurve3
curves[[4]] = sdd1
curves[[5]] = sdd2
curves[[6]] = sdd3
curves[[7]] = Opt$Age
# save(curves, file = file.path(plotdir,"curves-LCBC.Rda"))


#plot means and SD's
res = data.frame(res1,res2,res3) %>% mutate(zer=0)
ggplot(res) +
  geom_line(aes(x = age, y = meancurve2), col = col.cl.main[2], size =2,alpha=0.4) +
  geom_line(aes(x = age, y = meancurve1), col = col.cl.main[1], size = 2) +
  geom_line(aes(x = age, y = meancurve3), col = col.cl.main[3], size = 2) +
  geom_ribbon(aes(x = age, ymin = meancurve2-sdd2, ymax = meancurve2+sdd2), alpha = 0.2, fill = col.cl.main[2]) +
  geom_ribbon(aes(x = age, ymin = meancurve1-sdd1, ymax = meancurve1+sdd1), alpha = 0.7, fill = col.cl.main[1]) +
  geom_ribbon(aes(x = age, ymin = meancurve3-sdd3, ymax = meancurve3+sdd3), alpha = 0.5, fill = col.cl.main[3]) +
  geom_line(aes(x = age, y = zer), col = "black", size = 1, linetype = 3) +
  xlab("Age") +
  ylab("Mean (SD)") +
  theme_classic() +
  theme(text = element_text(size=16),
        plot.title = element_text(hjust=0.5))
