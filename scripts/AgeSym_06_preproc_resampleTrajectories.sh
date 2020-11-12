#!/bin/bash

########################################################################
## Purpose: resample significant asymmetry trajectories and hemi effects      
########################################################################

base=/Users/jamesroe/Dropbox/OpenScienceFramework/AgeSym

#--------
#c=${1} #LCBC / CamCan / BaseII / Betula
c=CamCan
#--------

cd $indir

#output
resdir="${base}/results/reproduceClustering/$c"
if [ -d $resdir ]; then
	echo "> removing results dir"
	rm -r $resdir
fi
echo "> copying to results"
mkdir -p $resdir
cp $base/reproduceClustering/$c/* $resdir
	

cd $resdir


echo "SOURCING FREESURFER"
export FREESURFER_HOME=/Applications/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=/Users/jamesroe/Dropbox/PHD_project/Paper3/Materials/subjectsdir
echo "SOURCING FSL"
FSLDIR=/Applications/fsl
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH


template=40tvs8ref_sym_20


if [ $c == "LCBC" ]; then
	#take LCBC FDR-corrected significance map as input
	#mask by p(FDR)=.001 and remove clusters <300mm2
	#LCBC hardcoded here
	echo "> masking map"
	mkdir surfclust
	mri_binarize --i mapSigFDR.LCBC.mgh --min 3.0 --o mSigFDR30.LCBC.mgh
	mri_surfcluster \
	--in mapSigFDR.LCBC.mgh \
	--subject $template \
	--hemi lh \
	--surf inflated \
	--mask mSigFDR30.LCBC.mgh \
	--olab surfclust/cluster \
	--minarea 300 \
	--o SigFDR30.LCBC.mgh \
	--thmin 3.0
	mri_binarize --i SigFDR30.LCBC.mgh --min 0.000000001 --o maskSigFDR30.LCBC.mgh 
fi


#mask hemi effect by p(FDR)=.001
mri_binarize --i mapHPFDR.mgh --abs --min 3.0 --o mHPFDR30.mgh
mris_calc --output mapHCoef30.mgh mapHCoef.mgh mul mHPFDR30.mgh


#make 4D mask to constrain clustering within LCBC results
for i in `seq 100`; do 
	echo "--i maskSigFDR30.LCBC.mgh \\"
done > tmp1
echo "#!/bin/bash" > tmphead
echo "concatdir=${1}" >> tmphead
echo "echo running mri_concat" >> tmphead
echo "mri_concat \\" >> tmphead
echo "--o m1004d.mgh" > tmptail
cat tmphead tmp1 tmptail > 1004d.sh
sh 1004d.sh
rm tmp*


#mask cohort-specific trajectories
mris_calc --output m1004dfit.mgh 1004dfit.mgh mul m1004d.mgh	


echo "> resampling to fsaverage"
from=40tvs8ref_sym_20; to=fsaverage
mris_apply_reg \
--src m1004dfit.mgh \
--trg m1004dfit.fsav.mgh \
--streg $SUBJECTS_DIR/$from/surf/lh.sphere.reg $SUBJECTS_DIR/$to/surf/lh.$from.sphere.reg
mris_apply_reg \
--src mapHCoef.mgh \
--trg mapHCoef.fsav.mgh \
--streg $SUBJECTS_DIR/$from/surf/lh.sphere.reg $SUBJECTS_DIR/$to/surf/lh.$from.sphere.reg
mris_apply_reg \
--src maskSigFDR30.LCBC.mgh \
--trg maskSigFDR30.LCBC.fsav.mgh \
--streg $SUBJECTS_DIR/$from/surf/lh.sphere.reg $SUBJECTS_DIR/$to/surf/lh.$from.sphere.reg


echo "> resampling to fsaverage5"
from=fsaverage; to=fsaverage5
mris_apply_reg \
--src m1004dfit.fsav.mgh \
--trg m1004dfit.fsav5.mgh \
--streg $SUBJECTS_DIR/$from/surf/lh.sphere.reg $SUBJECTS_DIR/$to/surf/lh.sphere.reg
mris_apply_reg \
--src mapHCoef.fsav.mgh \
--trg mapHCoef.fsav5.mgh \
--streg $SUBJECTS_DIR/$from/surf/lh.sphere.reg $SUBJECTS_DIR/$to/surf/lh.sphere.reg
mris_apply_reg \
--src maskSigFDR30.LCBC.fsav.mgh \
--trg maskSigFDR30.LCBC.fsav5.mgh \
--streg $SUBJECTS_DIR/$from/surf/lh.sphere.reg $SUBJECTS_DIR/$to/surf/lh.sphere.reg


#get vertices on fsaverage5
mri_binarize --i maskSigFDR30.LCBC.fsav5.mgh --o finalSigFDR30.LCBC.fsav5.mgh --min 0.1
mkdir surfclust.fsav5
mri_surfcluster \
--in finalSigFDR30.LCBC.fsav5.mgh  \
--subject fsaverage5 \
--hemi lh \
--surf inflated \
--olab surfclust.fsav5/cluster \
--thmin 1
lab="fsav5-ALL.label"
mri_mergelabels \
-d surfclust.fsav5 \
-o surfclust.fsav5/$lab
echo "> writing vertex index"
lc=`cat surfclust.fsav5/$lab | wc -l`
tail -n`echo $((lc-2))` surfclust.fsav5/$lab | awk '{print $1}' > vtxfsav5.csv
echo $((lc-2)) > nvtxfsav5.csv

rm *fsav.*