#!/bin/bash

#set -x

maindir=
outdir=

HPF=100
TR=2.3
del_vols=2
nr_vols=68 #vols - del_vols

###################################
# obtain regressors of WM and CSF #
###################################
cd $maindir

for i in `ls -d sub-???`; do
 #for ses in ses-BL; do #ses-TT ses-PT; do
 # for run in run-1 run-2; do
anatdir=$maindir/$i/anat
funcdir=$maindir/$i/func

outdir_sub=$outdir/"$i"
mkdir $outdir_sub

# combine WM and CSF regressor (tail + 5 because we have to remove 3 volume later)
awk -F"\t" '{print $1 " " $5}' $funcdir/"$i"_task-faces_run-1_desc-confounds_regressors.tsv | tail -n +4 > $outdir_sub/nuisance_timeseries

#remove first 2 volumes from BIDS output --> mean signal will be calculated for these all -2 volumes
fslroi $funcdir/"$i"_task-faces_run-1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold.nii.gz $outdir_sub/"$i"_task-faces_run-1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold_preproc_volsdel.nii.gz ${del_vols} ${nr_vols}

# create tempMean
fslmaths $outdir_sub/"$i"_task-faces_run-1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold_preproc_volsdel.nii.gz -Tmean $outdir_sub/"$i"_task-faces_run-1_tempMean.nii.gz

# Regress WM and CSF from your main signal
fsl_glm -i $outdir_sub/"$i"_task-faces_run-1_space-MNI152NLin6Asym_desc-smoothAROMAnonaggr_bold_preproc_volsdel.nii.gz -d $outdir_sub/nuisance_timeseries --demean --out_res=$outdir_sub/"$i"_task-faces_run-1_residual.nii.gz


############################
### BAND PASS FILTERING ####
############################

# high pass only
# bptf = high-pass filter / TR / 2
BPTF=`echo ${HPF} / ${TR} / 2 | bc -l`

fslmaths $outdir_sub/"$i"_task-faces_run-1_residual.nii.gz -bptf ${BPTF} -1 -add $outdir_sub/"$i"_task-faces_run-1_tempMean.nii.gz $outdir_sub/filtered_func_data.nii.gz

done
