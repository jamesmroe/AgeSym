#==============================================================#
# Purpose: cognitive change longitudinal analysis with LCBC data
#==============================================================#

args <- commandArgs(TRUE)
datadir <- as.character(args[1])


#load packages
#............#
tmp.packages = c("tidyverse","readr","devtools","itsadug","magrittr","gamm4","nlme")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
rm(list=ls(pattern="tmp*"))


#coded ROIs
labs=c("Exclfit3cl1.0001.label",
        "Exclfit3cl1.0002.label",
        "Exclfit3cl1.0003.label",
        "Exclfit3cl1.0008.label",
        "Exclfit3cl2.0001.label",
        "Exclfit3cl3.0001.label",
        "Exclfit3cl3.0002.label",
        "Exclfit3cl3.0003.label")


#load CVLT data
#.............#
longCVLT = read.csv(file.path(datadir,"Statfiles/longCVLT.csv"),stringsAsFactors = F,header=T, sep="\t")
longCVLT %<>%
  mutate(site = ifelse(Site_Name=="ousAvanto",0,1),
         sitescale = (site-mean(site))/sd(site),
         sex_demean = sex-mean(sex)
  )
#describe sample characteristics
length(unique(longCVLT$fsid_base))
nrow(longCVLT)
mean(longCVLT$nTimepoints)
mean(longCVLT$TimeBsl)
max(longCVLT$TimeBsl)
obsCVLT=longCVLT$Folder #observations


#main protocol - CVLT
#...................#
for (i in 1:length(labs)) {
  if (i ==1) {
    #set outputs
    CVLT.gam = list()
    CVLT.gamageout = list()
    CVLT.stable = list()
    CVLT.ptable = list()
    CVLT.stableageout = list()
    CVLT.ptableageout = list()
    hemidiffout1=list()
    meanCthout1=list()
  }
  whichlabel=labs[i]
  print(whichlabel)
  
  #subset
  DFF = DF %>% filter(Folder %in% obsCVLT)
  DFF$thickness = DFF[[whichlabel]]
  LH = DFF %>% filter(hemi == 1)
  RH = DFF %>% filter(hemi == 0)
  
  
  #compute asym and mean thickness
  longCVLT$hemidiff = LH$thickness-RH$thickness
  longCVLT$meanCth = (LH$thickness+RH$thickness) / 2
  
  
  #PCA on CVLT subtests
  PCA1 = prcomp(data.frame(longCVLT$CVLT_A_Total,
                           longCVLT$CVLT_5min_Free,
                           longCVLT$CVLT_30min_Free), center = T,scale. = T)
  summary(PCA1)
  Y = PCA1$x[,1]
  longCVLT$Y = Y*-1 #inverse PC1 to correlate negatively with age

  
  #GAMM models
  #age corrected
  mod1=gamm4(Y ~ s(hemidiff, k=6) + s(Age, k=6) + s(meanCth, k=6) + sex_demean + sitescale + VLT_Version, random = ~ (1|fsid_base), data = longCVLT)
  mod1sum=summary(mod1$gam)
  CVLT.gam[[i]]=mod1
  CVLT.stable[[i]] = mod1sum$s.table
  CVLT.ptable[[i]] = mod1sum$p.table
  
  #age uncorrected
  mod2=gamm4(Y ~ s(hemidiff, k=6) + s(meanCth, k=6) + sex_demean + sitescale + VLT_Version, random = ~ (1|fsid_base), data = longCVLT)
  mod2sum=summary(mod2$gam)
  CVLT.gamageout[[i]]=mod2
  CVLT.stableageout[[i]] = mod2sum$s.table
  CVLT.ptableageout[[i]] = mod2sum$p.table
  
  #get all asyms and means for F test
  hemidiffout1[[i]]=longCVLT$hemidiff
  meanCthout1[[i]]=longCVLT$meanCth
}
  

#load Matrix data
#...............#
longIQ = read.csv(file.path(datadir,"Statfiles/longIQ.csv"),stringsAsFactors = F,header=T, sep="\t")
longIQ %<>%
  mutate(site = ifelse(Site_Name=="ousAvanto",0,1),
         sitescale = (site-mean(site))/sd(site),
         sex_demean = sex-mean(sex)
  )
length(unique(longIQ$fsid_base))
nrow(longIQ)
mean(longIQ$nTimepoints)
mean(longIQ$TimeBsl)
max(longIQ$TimeBsl)
obsIQ=longIQ$Folder #observations


#main protocol - Matrix reasoning
#...............................#
for (i in 1:length(labs)) {
  if (i ==1) {
    #set outputs
    IQ.gam = list()
    IQ.gamageout = list()
    IQ.stable = list()
    IQ.ptable = list()
    IQ.stableageout = list()
    IQ.ptableageout = list()
    hemidiffout2 = list()
    meanCthout2 = list()
  }
  whichlabel=labs[i]; print(whichlabel)
  
  #subset
  DFF = DF %>% filter(Folder %in% obsIQ)
  DFF$thickness = DFF[[whichlabel]]
  LH = DFF %>% filter(hemi == 1)
  RH = DFF %>% filter(hemi == 0)
  
  
  #compute asym and mean thickness
  longIQ$hemidiff = LH$thickness-RH$thickness
  longIQ$meanCth = (LH$thickness+RH$thickness) / 2
  
  
  #GAMM models
  #age corrected
  mod1=gamm4(IQ_Matrix_Raw ~ s(hemidiff, k=6) + s(Age, k=6) + s(meanCth, k=6) + sex_demean + sitescale, random = ~ (1|fsid_base), data = longIQ)
  mod1sum=summary(mod1$gam)
  IQ.gam[[i]]=mod1
  IQ.stable[[i]] = mod1sum$s.table
  IQ.ptable[[i]] = mod1sum$p.table
  
  #age not corrected
  mod2=gamm4(IQ_Matrix_Raw ~ s(hemidiff, k=6) + s(meanCth, k=6) + sex_demean + sitescale, random = ~ (1|fsid_base), data = longIQ)
  mod2sum=summary(mod2$gam)
  IQ.gamageout[[i]]=mod2
  IQ.stableageout[[i]] = mod2sum$s.table
  IQ.ptableageout[[i]] = mod2sum$p.table

  #get all asyms and means for F test
  hemidiffout2[[i]]=longIQ$hemidiff
  meanCthout2[[i]]=longIQ$meanCth
}



#extract model outputs 
#....................#
#main models
for (i in 1:length(labs)) {
  if (i==1) { 
    #set outputs
    CVLTedf1=0; CVLTedf2=0; CVLTf1=0; CVLTf2=0; CVLTp1=0; CVLTp2=0; CVLTpsex=0; CVLTpsite= 0
    IQedf1=0; IQedf2=0; IQf1=0; IQf2=0; IQp1=0; IQp2=0; IQpsex=0; IQpsite= 0
  }
  CVLTedf1[i] = CVLT.stable[[i]][1,1]
  CVLTedf2[i] = CVLT.stable[[i]][3,1]
  CVLTf1[i] = CVLT.stable[[i]][1,3]
  CVLTf2[i] = CVLT.stable[[i]][3,3]
  CVLTp1[i] = CVLT.stable[[i]][1,4]
  CVLTp2[i] = CVLT.stable[[i]][3,4]
  CVLTpsex[i] = CVLT.ptable[[i]][2,4]
  CVLTpsite[i] = CVLT.ptable[[i]][3,4]
  
  
  IQedf1[i] = IQ.stable[[i]][1,1]
  IQedf2[i] = IQ.stable[[i]][3,1]
  IQf1[i] = IQ.stable[[i]][1,3]
  IQf2[i] = IQ.stable[[i]][3,3]
  IQp1[i] = IQ.stable[[i]] [1,4]
  IQp2[i] = IQ.stable[[i]] [3,4]
  IQpsex[i] = IQ.ptable[[i]][2,4]
  IQpsite[i] = IQ.ptable[[i]][3,4]
  
}

#age uncorrected
for (i in 1:length(labs)) {
  if (i==1) { 
    #set outputs
    CVLTageoutp1=0; CVLTageoutp2=0
    IQageoutp1=0; IQageoutp2=0
  }
  CVLTageoutp1[i] = CVLT.stableageout[[i]][1,4]
  CVLTageoutp2[i] = CVLT.stableageout[[i]][2,4]
  IQageoutp1[i] = IQ.stableageout[[i]] [1,4]
  IQageoutp2[i] = IQ.stableageout[[i]] [2,4]
}


#inspect [uncorrected] signficant effects 
#.......................................#
#(asym) 
# which(CVLTp1<0.05)
# which(IQp1<0.05)
# plot.gam(CVLT.gam[[3]]$gam, residuals=T)
# plot.gam(CVLT.gam[[5]]$gam, residuals=T)


#(mean thick)
# which(CVLTp2<0.05)
# which(IQp2<0.05)
# plot.gam(CVLT.gam[[4]]$gam, residuals=T)
# plot.gam(IQ.gam[[2]]$gam, residuals=T)
# plot.gam(IQ.gam[[8]]$gam, residuals=T)


#(asym age uncorrected)
# which(CVLTageoutp1<0.05)
# which(IQageoutp1<0.05)
# plot.gam(CVLT.gamageout[[6]]$gam, residuals=T,col="goldenrod1", shade=T, shade.col="goldenrod2")
# plot.gam(IQ.gamageout[[1]]$gam, residuals=T,col="darkolivegreen1", shade=T, shade.col="darkolivegreen3",cex=3)
# plot.gam(IQ.gamageout[[3]]$gam, residuals=T,col="darkolivegreen1", shade=T, shade.col="darkolivegreen3",cex=3)
# plot.gam(IQ.gamageout[[4]]$gam, residuals=T,col="darkolivegreen1", shade=T, shade.col="darkolivegreen3",cex=3)
# plot.gam(IQ.gamageout[[6]]$gam, residuals=T,col="darkolivegreen1", shade=T, shade.col="darkolivegreen3",cex=3)


#(mean thick age uncorrected)
# for (i in 1:length(labs)) {
#   plot.gam(CVLT.gamageout[[i]]$gam, select=2,residuals=T,col="cyan", shade=T, shade.col="cyan",cex=3)
# }
# for (i in 1:length(labs)) {
#   plot.gam(IQ.gamageout[[i]]$gam, select=2,residuals=T,col="pink4", shade=T, shade.col="pink4",cex=3)
# }


#FDR correction
pFDR=p.adjust(c(CVLTp1,IQp1),method="BH")


#build table output
#.................#
newlabs=c("rostral anterior cingulate",
          "superior frontal",
          "lateral orbitofrontal",
          "superior temporal",
          "frontal cortex",
          "insula",
          "inferior parietal",
          "caudal anterior cingulate"
)
solution = c(rep("Cluster1",4),rep("Cluster2",1),rep("Cluster3",3))

#CVLT table
out1=data.frame(labs,
                labnames = c(1:8),
                newlabs,
                solution,
                "edf" = format(round(CVLTedf1,2),nsmall=2),
                "F" = round(CVLTf1,2),
                "p"= round(CVLTp1,3),
                "pFDR"= round(pFDR[1:8],3),
                "meas" = rep("CVLT",8),
                "p.mThick"= ifelse(CVLTp2<0.001,"<0.001",round(CVLTp2,3)),
                "p.sex"= ifelse(CVLTpsex<0.001,"<0.001",round(CVLTpsex,3)),
                "p.site" = ifelse(CVLTpsite<0.001,"<0.001",round(CVLTpsite,3)),
                "p.ageout" = round(CVLTageoutp1,3),
                "p.ageoutmThick" = round(CVLTageoutp2,3)
)

# IQ table
out2=data.frame(labs,
                labnames = c(1:8),
                newlabs,
                solution,
                "edf" = format(round(IQedf1,2),nsmall=2),
                "F" = round(IQf1,2),
                "p"= round(IQp1,3),
                "pFDR"= round(pFDR[9:16],3),
                "meas" = rep("Matrix Reasoning",8),
                "p.mThick"= ifelse(IQp2<0.001,"<0.001",round(IQp2,3)),
                "p.sex"= ifelse(IQpsex<0.001,"<0.001",round(IQpsex,3)),
                "p.site" = ifelse(IQpsite<0.001,"<0.001",round(IQpsite,3)),
                "p.ageout" = round(IQageoutp1,3),
                "p.ageoutmThick" = round(IQageoutp2,3)
)



# Ftest of global model fits #
#............................#
# PCA on mean thickness - CVLT
pcadat=data.frame(l1=meanCthout1[[1]],
                  l2=meanCthout1[[2]],
                  l3=meanCthout1[[3]],
                  l4=meanCthout1[[4]],
                  l5=meanCthout1[[5]],
                  l6=meanCthout1[[6]],
                  l7=meanCthout1[[7]],
                  l8=meanCthout1[[8]]
)
pcaCth=prcomp(pcadat[,],center=T,scale.=T)
spcaCth=summary(pcaCth)
pcomp1meanCth=pcaCth$x[,1] #PC1
pcomp2meanCth=pcaCth$x[,2] #PC2

longCVLT=data.frame(longCVLT,
                    s1=hemidiffout1[[1]],
                    s2=hemidiffout1[[2]],
                    s3=hemidiffout1[[3]],
                    s4=hemidiffout1[[4]],
                    s5=hemidiffout1[[5]],
                    s6=hemidiffout1[[6]],
                    s7=hemidiffout1[[7]],
                    s8=hemidiffout1[[8]])

# CVLT model with all asymmetry predictors
modcom1 = gamm4(Y ~ s(s1, k=6) +
                  s(s2, k=6) +
                  s(s3, k=6) +
                  s(s4, k=6) +
                  s(s5, k=6) +
                  s(s6, k=6) +
                  s(s7, k=6) +
                  s(s8, k=6) +
                  s(Age, k=6) + s(pcomp1meanCth) + s(pcomp2meanCth) + sex_demean + sitescale, random = ~ (1|fsid_base), data=longCVLT)
modcom2 = gamm4(Y ~
                  s(Age, k=6) + s(pcomp1meanCth) + s(pcomp2meanCth) + sex_demean + sitescale, random = ~ (1|fsid_base), data=longCVLT)

#compare models CVLT
anova(modcom1$mer,modcom2$mer,test="F")



# PCA on mean thickness - Matrix
pcadat=data.frame(l1=meanCthout2[[1]],
                  l2=meanCthout2[[2]],
                  l3=meanCthout2[[3]],
                  l4=meanCthout2[[4]],
                  l5=meanCthout2[[5]],
                  l6=meanCthout2[[6]],
                  l7=meanCthout2[[7]],
                  l8=meanCthout2[[8]]
)
pcaCth=prcomp(pcadat[,],center=T,scale.=T)
spcaCth=summary(pcaCth)
pcomp1meanCth=pcaCth$x[,1] #PC1
pcomp2meanCth=pcaCth$x[,2] #PC2

longIQ=data.frame(longIQ,
                  s1=hemidiffout2[[1]],
                  s2=hemidiffout2[[2]],
                  s3=hemidiffout2[[3]],
                  s4=hemidiffout2[[4]],
                  s5=hemidiffout2[[5]],
                  s6=hemidiffout2[[6]],
                  s7=hemidiffout2[[7]],
                  s8=hemidiffout2[[8]])


modcom3 = gamm4(IQ_Matrix_Raw ~ s(s1, k=6) +
                  s(s2, k=6) +
                  s(s3, k=6) +
                  s(s4, k=6) +
                  s(s5, k=6) +
                  s(s6, k=6) +
                  s(s7, k=6) +
                  s(s8, k=6) +
                  s(Age, k=6) + s(pcomp1meanCth) + s(pcomp2meanCth) + sex_demean + sitescale, random = ~ (1|fsid_base), data=longIQ)
modcom4 = gamm4(IQ_Matrix_Raw ~ 
                  s(Age, k=6) + s(pcomp1meanCth) + s(pcomp2meanCth) + sex_demean + sitescale, random = ~ (1|fsid_base), data=longIQ)

#compare models Matrix
anova(modcom3$mer,modcom4$mer,test="F")
