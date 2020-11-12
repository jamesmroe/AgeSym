#########################################################
## Purpose: merge vertex-wise GAMM results and write maps
#########################################################

args = commandArgs(TRUE)
outdir = as.character(args[1])
nvtx = 155865 #N vertices in lh.cortex.label
savefiles = 0


# Load packages
#.............#
packages = c("dplyr", "stringr", "magrittr")
# sapply(packages, install.packages, character.only = T)
sapply(packages, require, character.only = T)


# Check inputs
#............#
split.dir = file.path(outdir, "split.data")
setwd(split.dir)
MMs = list.files(pattern = "label.matrix")
RRs = list.files(pattern = "gamm.results")
if (!length(RRs) == length(MMs)) {
  print("the number of LabelxSub matrices and GAMM results files are not equal. Quitting"); quit()
}


# Merge results
#.............#
for (i in 1:length(RRs)) {
  j = RRs[i]
  
  if (i == 1) {
    #set outputs
    edf = Fval = p.spl = pp.spl = NULL 
    fit_valL = fit_valR = fit_valdiff = NULL
    hT = hCoef = hP = hPlog = NULL 
    pb = txtProgressBar(min=1, max=length(RRs), style=3)
  }
  setTxtProgressBar(pb,i)
  load(j)
  edf %<>% cbind(.,RR$edf)
  Fval %<>% cbind(.,RR$Fval)
  p.spl %<>% cbind(.,RR$p.spl)
  pp.spl %<>% cbind(.,RR$pp.spl)
  fit_valL %<>% cbind(.,RR$fit_valL)
  fit_valR %<>% cbind(.,RR$fit_valR)
  fit_valdiff %<>% cbind(.,RR$fit_valdiff)
  hT %<>% cbind(.,RR$hT)
  hCoef %<>% cbind(.,RR$hCoef)
  hP %<>% cbind(.,RR$hP)
  hPlog %<>% cbind(.,RR$hPlog)
}
vars = c("edf","Fval","p.spl","pp.spl",
         "fit_valL","fit_valR","fit_valdiff",
         "hT","hCoef","hP","hPlog")
for (v in vars) {
  if (!dim(get(v))[2] == nvtx) { print(paste(v,"does not match with n of vertices")); quit() }
}


#FDR correction
#Age x Hemi effects
FDRcorr = p.adjust(pp.spl,method = "BY")


#fix significance ceiling and get direction of hemi effects
hPlog[hP== 0.000000e+00] = ceiling(max(hPlog[is.finite(hPlog)]))+1
hPlog[which(hT < 0)] = hPlog[which(hT < 0)]*-1


#FDR correction
#Hemi effects
tmpHp = p.adjust(hP,method = "BY")
logHpcor = -log10(tmpHp)
logHpcor[hP== 0.000000e+00] = ceiling(max(logHpcor[is.finite(logHpcor)]))+1
logHpcor[which(hT < 0)] = logHpcor[which(hT < 0)]*-1


if (savefiles==1) {
  print("saving csv files of maps")
  
  save("edf", "Fval", "p.spl", "fit_valL", "fit_valR", "fit_valdiff", "fit_derivdiff",
              file = file.path(outdir, "gamm_summary.Rda"))
  
  write.table(-log10(FDRcorr),
              file = file.path(outdir, "mapSigFDR.csv"), sep=",", col.names = F, row.names = F, quote = F)
  #uncorrected signficance
  write.table(p.spl, 
              file = file.path(outdir, "mapSig.csv"), sep=",", col.names = F, row.names = F, quote = F)
  #F stats
  write.table(Fval, 
              file = file.path(outdir, "mapFval.csv"), sep=",", col.names = F, row.names = F, quote = F)
  #edf
  write.table(edf, 
              file = file.path(outdir, "mapEdf.csv"), sep=",", col.names = F, row.names = F, quote = F)
  #Left trajectory
  write.table(t(fit_valL), 
              file = file.path(outdir, "fitValL.csv"), sep=",", col.names = F, row.names = F, quote = F)
  #Right trajectory
  write.table(t(fit_valR), 
                file = file.path(outdir, "fitValR.csv"), sep=",", col.names = F, row.names = F, quote = F)
  #diff trajectory
  write.table(t(fit_valdiff), 
                file = file.path(outdir, "fitValDiff.csv"), sep=",", col.names = F, row.names = F, quote = F)
  
  
  #hemi effects
  write.table(hT, 
              file = file.path(outdir, "mapHT.csv"), sep=",", col.names = F, row.names = F, quote = F)
  write.table(hCoef, 
              file = file.path(outdir, "mapHCoef.csv"), sep=",", col.names = F, row.names = F, quote = F)
  write.table(hPlog, 
              file = file.path(outdir, "mapHP.csv"), sep=",", col.names = F, row.names = F, quote = F) 
  write.table(logHpcor,
              file = file.path(outdir, "mapHPFDR.csv"), sep=",", col.names = F, row.names = F, quote = F)
}

quit()
