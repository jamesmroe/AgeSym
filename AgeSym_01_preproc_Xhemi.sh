#!/bin/bash

#######################################################################
## Purpose: register individual left and right thickness maps to LH_Sym
#######################################################################

echo "SOURCING FREESURFER"
export FREESURFER_HOME=${1}
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=${2}

sub=${3}
subfolder_recon=$SUBJECTS_DIR/$sub
file=thickness
template=40tvs8ref_sym_20


#Xhemi (implies opposite hemisphere)
if [ ! -e "$subfolder_recon/surf/lh.$template.sphere.reg" ]; then
	surfreg --s $sub --t $template --lh
fi

if [ ! -e "$subfolder_recon/xhemi/surf/lh.$template.sphere.reg" ]; then			
	surfreg --s $sub --t $template --lh --xhemi
fi


#Register data to sym
#LEFT to LEFT
mris_apply_reg \
--src $SUBJECTS_DIR/$sub/surf/lh.$file \
--trg $datadir/lh.lh.${file}.${sub}.${template}.mgh \
--streg $SUBJECTS_DIR/$sub/surf/lh.$template.sphere.reg $SUBJECTS_DIR/$template/surf/lh.sphere.reg
						
#RIGHT to LEFT
mris_apply_reg \
--src $SUBJECTS_DIR/$sub/surf/rh.$file \
--trg $datadir/rh.lh.${file}.${sub}.${template}.mgh \
--streg $SUBJECTS_DIR/$sub/xhemi/surf/lh.$template.sphere.reg $SUBJECTS_DIR/$template/surf/lh.sphere.reg


sm=8
for hemi in lh rh; do
echo "> smoothing data with FWHM ${sm}"
	mri_surf2surf \
	--srcsurfval $datadir/$hemi.lh.${file}.${sub}.${template}.mgh \
	--fwhm-trg $sm \
	--srcsubject $template \
	--trgsubject $template \
	--hemi lh \
	--surfreg sphere.reg \
	--trgsurfval $datadir/$hemi.lh.${file}.${sub}.${template}.sm${sm}.mgh \
	--cortex \
	--noreshape
	
	if [ -e $datadir/$hemi.lh.${file}.${sub}.${template}.sm${sm}.mgh ]; then
		rm $datadir/$hemi.lh.${file}.${sub}.${template}.mgh
	fi
done


#make new cortex label excluding corpus callosum
for hemi in lh; do
	if [ ! -e "$SUBJECTS_DIR/$template/label/${hemi}.CORTEX.label" ]; then
		echo "> making new ${hemi}.cortex.label"
		mv $SUBJECTS_DIR/$template/label/${hemi}.cortex.label $SUBJECTS_DIR/$template/label/${hemi}.CORTEX.label
		
		mri_annotation2label \
		--annotation aparc.DKTatlas40 \
		--subject $template \
		--hemi ${hemi} \
		--surf inflated \
		--sd $SUBJECTS_DIR \
		--outdir $SUBJECTS_DIR/$template/label/aparc.DKTatlas40/${hemi}

		mri_mergelabels \
		-d $SUBJECTS_DIR/$template/label/aparc.DKTatlas40/${hemi} \
		-o $SUBJECTS_DIR/$template/label/${hemi}.cortex.label
		
	fi
done
