#!/bin/bash

########################################################################
## Purpose: resample significant asymmetry trajectories and hemi effects      
########################################################################

outdir=${1}
cd $outdir

echo "SOURCING FREESURFER"
export FREESURFER_HOME=/cluster/projects/p23/tools/mri/freesurfer/current
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=$FREESURFER_HOME/subjects
echo "SOURCING FSL"
FSLDIR=/cluster/projects/p23/tools/mri/fsl/current
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH


template=40tvs8ref_sym_20


#mask by p(FDR)=.001 and remove clusters <300mm2
echo "> masking map"
mkdir surfcluster
mri_binarize --i mapSigFDR.mgh --min 3.0 --o mSigFDR30.mgh
mri_surfcluster \
--in mapSigFDR.mgh \
--subject $template \
--hemi lh \
--surf inflated \
--mask mSigFDR30.mgh \
--olab surfclust/cluster \
--minarea 300 \
--o SigFDR30.mgh \
--thmin 3.0
mri_binarize --i SigFDR30.mgh --min 0.000000001 --o maskSigFDR30.mgh 
	

#mask hemi effect by p(FDR)=.001
mri_binarize --i mapHPFDR.mgh --abs --min 3.0 --o mHPFDR30.mgh
mris_calc --output mapHCoef30.mgh mapHCoef.mgh mul mHPFDR30.mgh


#make 4D mask
for i in `seq 100`; do 
	echo "--i maskSigFDR30.mgh \\"
done > tmp1
echo -e '#!/bin/bash\n\nconcatdir=${1}\n\necho "running mri_concat"\n\nmri_concat \' > tmphead
echo -e "--o m1004d.mgh" > tmptail
cat tmphead tmp1 tmptail > 1004d.sh
sh 1004d.sh
rm tmp*


#mask significant trajectories
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
--src maskSigFDR30.mgh \
--trg maskSigFDR30.fsav.mgh \
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
--src maskSigFDR30.fsav.mgh \
--trg maskSigFDR30.fsav5.mgh \
--streg $SUBJECTS_DIR/$from/surf/lh.sphere.reg $SUBJECTS_DIR/$to/surf/lh.sphere.reg


#get vertices on fsaverage5
mri_binarize --i maskSigFDR30.fsav5.mgh --o finalSigFDR30.fsav5.mgh --min 0.1
mkdir surfclust.fsav5
mri_surfcluster \
--in finalSigFDR30.fsav5.mgh  \
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
tail -n`echo $((lc-2))` surfclust.fsav5/$lab | awk '{print $1}' > $outdir/vtxfsav5.csv
echo $((lc-2)) > $outdir/nvtxfsav5.csv

