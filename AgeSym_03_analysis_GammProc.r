########################################################
## Purpose: run vertex-wise factor-smooth GAMM models to 
##          compute smooth Age x Hemisphere interaction
########################################################


# Script inputs
#.............#
i=1
base="/Users/jamesroe/Dropbox/OpenScienceFramework/AgeSym/vtxGamm"
subset.size = 3                              #N vertices to compute in subset
outdir = base                                
database = file.path(base,"simdat.csv")      #simulated data for 3 vertices
nvtx = 155865		                             #N vertices in lh.cortex.labels
knots = 6		                                 #N knots used for GAMM fitting
print(paste(subset.size,outdir,database,nvtx,knots))


# Load packages
#.............#
packages = c("dplyr", "stringr", "numDeriv","gamm4","magrittr","ggplot2")
# sapply(packages, install.packages, character.only = T)
sapply(packages, require, character.only = T)


# Load VertexXSub mat
#..................#
vtxmat = read.csv(file.path(base, "simvtx.csv"), header = F) %>%
  t() %>% as.data.frame()


# Load database (ordered same as LabelXSub mat)
#.............................................#
# required variables:ID, Age, Hemisphere, Sex, Scanner
db = read.csv(database, sep = "\t", header = T, stringsAsFactors = F) %>%
  mutate(scanner_demean = ifelse(Site_Name == "ousAvanto", -0.5, 0.5),
         sex_demean=Sex-mean(Sex),
         scanner_demean = (scanner_demean-mean(scanner_demean)) /sd(scanner_demean)) %>%
  select(fsid_base,Age,hemi,sex_demean,scanner_demean)



Yorig = as.matrix(vtxmat)
N = ceiling(nvtx/subset.size)
print(N)
if (i == N) {
  end = dim(Yorig)[2]
} else {
  end = subset.size
}
print(paste("indexing", 1, ":", end))
Y = Yorig[,1:end]


# make 100 obs age model
#......................#
nn = 100
Opt = list()
AGE = db$Age #true age
Opt$Age = seq(min(AGE), max(AGE), length.out = nn) 
Opt$fake.frame = data.frame("Age" = Opt$Age,
                            "hemi" =rep(1,nn),
                            "sex_demean" = rep(0,nn),
                            "scanner_demean" = rep(0,nn))
pdat = rbind(Opt$fake.frame,
             Opt$fake.frame)
pdat$hemi[1:100] = 0 #prediction data


# filename=file.path(outdir,"Opt.Rda")
# if (!file.exists(filename)) {
#   save(Opt,file =filename)
# }



# vertex-wise GAMMs
#.................#
for (j in 1:end) {
  
  if (j == 1) {
    #set outputs
    gamm.trajectories = list()
    ogamm.trajectories = list()
    RR = list()
    pb = txtProgressBar(min=1, max=end, style=3)
    ph = list()
  }
  setTxtProgressBar(pb,j)
  
  
  #select vertex
  Y = data.frame((Yorig)[,j])
  names(Y) = "Y"
  dat = data.frame(Y,db)
  
  
  #factor-smooth GAMM models
  gamm.trajectories[[j]] = gamm4(Y ~ s(Age, by = as.factor(hemi), k = knots) + as.factor(hemi) + sex_demean + scanner_demean, 
                                 data = dat, random = ~ (1 |fsid_base))
  gamm.sum = summary(gamm.trajectories[[j]]$gam)
  g = gamm.trajectories[[j]]$gam
  
  
  #check worked - zero-centered hemispheric age trajectories
  #plot.gam(g,shade = T,  pages = 1, scale = 0, seWithMean = T)
  
  
  #simulate from posterior - Xp matrix
  #muliply matrix by model coefficients and sum rows to get predicted values
  Xp = predict(g, newdata = pdat, type = "lpmatrix") 
  #plot(Xp %*% coef(g))  # fitted trajectories
  
  
  #which cols of Xp relate to which hemi
  c1 = grepl("hemi\\)1", colnames(Xp)) #L
  c2 = grepl("hemi\\)0", colnames(Xp)) #R
  
  
  #which rows of Xp
  r1 = with(pdat, hemi == 1) #L
  r2 = with(pdat, hemi == 0) #R
  
  
  #differences between smooths
  #cols of differenced Xp matrix that aren't involved in comparison set to zero
  X = Xp[r1, ] - Xp[r2, ] #L-R
  X[, !(c1 | c2)] = 0
  X[, !grepl("s\\(", colnames(Xp))] = 0
  
  
  #get predicted vals of zero-centered hemispheric difference at zero of other covs
  #SE of the difference using variance-covariance matrix of estimated model coefficients
  dif = X %*% coef(g)
  se = sqrt(rowSums((X %*% vcov(g, unconditional =T)) * X))
  comp = data.frame("OptAge" = Opt$Age, dif, se)
  
  
  #ordered factor approach to get test statistics for GAMM interaction
  dat = mutate(dat,
               ohemi = ifelse(dat$hemi == 1, "left", "right"),
               ohemi = factor(ohemi, levels = c("left","right"),ordered = T))
  
  
  #estimate smooth for set reflevel and a smoothed difference between ref and other levels
  ogamm.trajectories[[j]] = gamm4(Y ~ as.factor(hemi) + s(Age) + s(Age, by = ohemi, k = knots) + sex_demean + scanner_demean, 
                                  data = dat, random = ~ (1 |fsid_base))
  ogamm.sum = summary(ogamm.trajectories[[j]]$gam)
  
  
  # #retrieve data used by plotgam
  plotData <- list()
  trace(mgcv:::plot.gam, at = list(c(27, 1)),
        quote({
          message("assigning into globalenv()'s plotData...")
          plotData <<- pd
        }))
  
  
  mgcv::plot.gam(gamm.trajectories[[j]]$gam, seWithMean = TRUE, pages = 1)
  plotdf=data.frame(Age=plotData[[1]]$x,
                    Lfit=plotData[[2]]$fit, #L
                    Rfit=plotData[[1]]$fit, #R
                    xl=plotData[[1]]$xlim,
                    Lse=plotData[[2]]$se,
                    Rse=plotData[[1]]$se
  )
  
  
  ph[[j]] = ggplot(plotdf) +
    geom_line(aes(Age,Lfit),col=" blue",alpha=0.7,size=1) +
    geom_line(aes(Age,Rfit),col="gold2",size=1) +
    # geom_line(aes(Age,Rfit-Lfit),col="black",size=1) +
    # geom_line(aes(Age,Lfit-Rfit),col="black",size=1) +
    geom_ribbon(aes(Age,ymin=Lfit-Lse,ymax = Lfit+Lse), alpha = 0.2,fill="blue") +
    geom_ribbon(aes(Age,ymin=Rfit-Rse,ymax = Rfit+Rse), alpha = 0.3,fill="gold2") +
    geom_line(data=comp,aes(Opt$Age,dif),col="#008571",size=1) +
    geom_ribbon(data=comp,aes(x=Opt$Age,ymin = dif-se, ymax = dif+se), alpha = 0.2,fill="#008571") +
    theme_classic() +
    theme(text=element_text(size=18)) +
    geom_text(y=-0.1, x=40, label = "s(LH-Age)", col="blue") +
    geom_text(y=-0.1-0.025, x=40, label = "s(RH-Age)", col="gold3") +
    geom_text(y=-0.1-0.05, x=40, label = "s(LH-Age)-s(RH-Age)", col="#008571") +
    ylab(NULL)
  
  
  #get stats
  RR$edf = RR$edf %>% cbind(., ogamm.sum$edf[2]) #edf
  RR$Fval = RR$Fval %>% cbind(., ogamm.sum$s.table[[2,3]]) #Fstat
  RR$p.spl = RR$p.spl %>% cbind(., -log10(ogamm.sum$s.pv[2])) #-log10 p-val for smooth difference
  RR$pp.spl = RR$pp.spl %>% cbind(., ogamm.sum$s.pv[2]) #p-val
  
  
  #get fitted trajectories
  RR$fit_valL = RR$fit_valL %>% cbind(.,
                                      predict.gam(gamm.trajectories[[j]]$gam, newdata = pdat[101:200,])) #L
  RR$fit_valR = RR$fit_valR %>% cbind(.,
                                      predict.gam(gamm.trajectories[[j]]$gam, newdata = pdat[1:100,])) #R
  RR$fit_valdiff = RR$fit_valdiff %>% cbind(., 
                                            dif) #zero-centered L-R difference
  RR$se = RR$se %>% cbind(., se)
  
  
  #get main effect of Hemisphere 
  RR$hT = RR$hT %>% cbind(., gamm.sum$p.t[[2]])
  RR$hCoef = RR$hCoef %>% cbind(., gamm.sum$p.coeff[[2]])
  RR$hP = RR$hP %>% cbind(., gamm.sum$p.pv[[2]])
  RR$hPlog = RR$hPlog %>% cbind(., -log10(gamm.sum$p.pv[[2]]))
  
}
ph[[3]]

# save(RR, file = file.path(split.dir, paste0("gamm.results.",str_pad(toString(i),3,pad = "0"),".Rda")))
# quit()              


#for further reading see Gavin Simpson's smooth term comparison procedure outlined here
#https://fromthebottomoftheheap.net/2017/10/10/difference-splines-i/
#https://fromthebottomoftheheap.net/2017/10/10/difference-splines-ii/
