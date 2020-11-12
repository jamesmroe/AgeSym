#========================================================#
# Purpose: use longitudinal AIBL data to assess 
#         group differences in asymmetry change over time 
#         between normal controls and Alzheimer's disease
#========================================================#

args <- commandArgs(TRUE)
datadir <- as.character(args[1])


#---load packages
tmp.packages = c("tidyverse","magrittr","lubridate")
tmpnew.packages = tmp.packages[!(tmp.packages %in% installed.packages()[,"Package"])]
if(length(tmpnew.packages)) {
  install.packages(tmpnew.packages)
}
sapply(tmp.packages, require, character.only = T)
rm(list=ls(pattern="tmp*"))


#---data (restricted) - AIBL data is open
#data is available at https://aibl.csiro.au/research/support/ pending application approval and compliance with the data usage agreement
long=read.csv(file.path(datadir,"Statfiles/labelbase_AIBLFULLclean.tsv"),stringsAsFactors = F, header=T, sep=" ")


#calc age_interval
long = long %>% arrange(RID) %>% 
  group_by(RID) %>% mutate(TimeBsl = AgeJitt - min(AgeJitt),
                           AgeBsl = min(AgeJitt),
                           nTimepoints = length(unique(ID)),
                           Site = ifelse(SITEID==1,0,1),
                           sex = ifelse(Sex=="f",0,1)) %>% 
  arrange(RID) %>% ungroup()


#set group
long %<>% mutate(group = includelong)
long$group = factor(long$group, c("NC","AD"))


#plot longitudinal trajectories
tmp1=long %>% filter(group=="AD")
tmp2=long %>% filter(group=="NC")
ggplot() +
  geom_line(data=long,aes(x=VISID,y=AIBLDIAGNOSIS,group=RID,col=group),position=position_jitter(w=0.08, h=0.08),alpha=0) +
  geom_line(data=tmp1,aes(VISID,AIBLDIAGNOSIS,group=RID,col=group),position=position_jitter(w=0.08, h=0.08),alpha=0.6) +
  geom_line(data=tmp2,aes(VISID,AIBLDIAGNOSIS,group=RID),position=position_jitter(w=0.1, h=0.1),col="#d8b707",alpha=0.1) +
  scale_color_manual(values=rev(c("#01401c","#d8b707"))) +
  labs(x = "Timepoint",  y="Diagnosis") +
  ggtitle("AIBL (Longitudinal sample)") +
  theme(
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black", size = 18, family = "Helvetica Neue Light"),
    axis.ticks = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
    axis.title.x = element_text(color = "black", size = 22, vjust = 1, margin = margin(0,20,0,0)),
    axis.text = element_text(color = "black", size = 18))


#coded ROI's
labs = c("Exclfit3cl1.0001.label",
         "Exclfit3cl1.0002.label",
         "Exclfit3cl1.0003.label",
         "Exclfit3cl1.0008.label",
         "Exclfit3cl2.0001.label",
         "Exclfit3cl3.0001.label",
         "Exclfit3cl3.0002.label",
         "Exclfit3cl3.0003.label")


# FUNCTION to get SE for LME prediction
predict.SE.lme = function(model.lme, pdat) {
  # predict fit
  pdat$pred = predict(model.lme, pdat, level = 0)
  # predict se.fit
  designmat = model.matrix(eval(eval(model.lme$call$fixed)[-2]),pdat)
  predvar = diag(designmat %*% mod1$varFix %*% t(designmat))
  pdat$SE = sqrt(predvar)
  return(list(pdat,
              designmat,
              predvar))
}


#---model per ROI
coefout=list()
pp=list()
long = long[!is.na(long$group),]
for (i in 1:length(labs)) {
  whichlabel=labs[[i]]
  long$thickasym = long[[whichlabel]]
  
  Y = long$thickasym
  long$group = as.factor(long$group)
  long$group = relevel(long$group, "NC")
  
  #mean center covariates
  long %<>% mutate(csex = sex - mean(sex),
                     site = Site - mean(Site)/sd(Site)) 
  
  #LME model
  mod1=nlme::lme(Y ~ group*TimeBsl + AgeBsl + csex + site,
                 random = ~ 1|factor(RID),
                 data= long)
  modsum1=summary(mod1)
  
  
  #inspect coefficients
  coefout[[i]] = coef(modsum1)
  

  #predict LME trajectories
  pdat = long %>% select(RID,group,TimeBsl,AgeBsl,csex,site) %>% 
    mutate(csex=0,
           AgeBsl=mean(AgeBsl),
           site=0)
  pred1=predict(mod1,newdata=pdat,level=0,se=T)
  pdat$pred1=pred1
  newdat = predict.SE.lme(mod1, pdat)[[1]]
  newdat$zer = 0
  newdat$CI = 1.96*newdat$SE
  
  
  #---roi plots
  pp[[i]]=ggplot() +
    geom_line(data=newdat,aes(x=TimeBsl,y=zer),color="black",size=0.75,linetype = 2,alpha =1) +
    geom_line(data=newdat, aes(x=TimeBsl,y=pred,col=group),alpha=1, size = 1.5) +
    geom_line(data=newdat, aes(x=TimeBsl,y=pred-CI,col=group),alpha=0.7, size = 0.5, linetype=1) +
    geom_line(data=newdat, aes(x=TimeBsl,y=pred+CI,col=group),alpha=0.7, size = 0.5, linetype=1) +
    geom_ribbon(data=newdat, aes(x=TimeBsl,ymin=pred-CI,ymax=pred+CI,fill=group), alpha=0.2) +
    scale_color_manual(values=rev(c("#01401c","#d8b707"))) +
    scale_fill_manual(values=rev(c("#01401c","#d8b707"))) +
    ggtitle(whichlabel) +
    labs(x = "Years", y = "LH-RH thickness") +
    theme(
      plot.background = element_rect(fill = "white"),
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      text = element_text(color = "black", size = 18, family = "Helvetica Neue Light"),
      plot.title = element_text(hjust = 0.5),
      axis.ticks = element_blank(),
      axis.line = element_line(color = "black", size=1),
      axis.title.y = element_text(color = "black", size = 22, vjust =-1, margin = margin(0,20,0,0)),
      axis.title.x = element_text(color = "black", size = 22, vjust = 0, margin = margin(0,20,0,0)),
      axis.text = element_text(color = "black", size = 22))
}


#---build table
newlabs=c("rostral anterior cingulate",
          "superior frontal",
          "lateral orbitofrontal",
          "anterior temporal",
          "frontal cortex",
          "insula",
          "inferior parietal",
          "caudal anterior cingulate")

solution = c(rep("1",4),rep("2",1),rep("3",3))
df1=t1=p1=c=0
for (i in 1:length(labs)) {
  c=c+1
  #extract model stats for group*time
  df1[c] = coefout[[i]][7,3]
  p1[c] = coefout[[i]][7,5]
  t1[c] = coefout[[i]][7,4]
}

#FDR correction
pfdr1=p.adjust(p1,"BH")
out1=data.frame("ROI" = c(1:8),
                "Clustering.solution" = solution,
                "df" = df1,
                "ADvControls.T" = round(t1,2),
                "ADvControls.p" = ifelse(p1<0.001,"<0.001",round(p1,3)),
                "ADvControls.pFDR" = ifelse(pfdr1<0.001,"<0.001",round(pfdr1,2)),
                "Sig" = ifelse(pfdr1 <=.05, "*","")
                )
out1
