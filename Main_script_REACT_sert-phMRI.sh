
##################################################################################################
#### READ ME ####

# This script contains all steps for the REACT analysis for sert-phMRI. If parts of the analysis have to be run again, please copy the parts in question from this script into a new slurm bash script.

# The first part contains all steps followed for the rs-fMRI REACT analysis 
# The second part contains all steps followed for the task-fMRI PPI REACT analysis

##################
##################################################################################################


#######################
#### rs-fMRI REACT ####
#######################

# Steps in the script:
# 1. Post-fMRIprep processing
#	1.1 Regressing out CSF and WM
#	1.2 Band pass filtering
# 2. REACT
# 3. Randomise
# 	3.1 Main effects
#		3.1.1 Main effects all subs
#		3.1.2 Main effects per group
# 	3.2 Dose-dependent effects
#		3.2.1 T-tests
#		3.2.2 
# 	3.3 Correlation SPECT maps
# 4. Extracting cluster location and individual values
#	4.1 Main effects
#	4.2 Main effects per group
#	4.3 Dose-dependent effects
#
##################################################################################################

USER=
type=fMRI

####################################
### Step 1: Post-fMRI processing ###
####################################

./$scriptdir/Postprocessing_rsfmri.sh # added as a separate script

#####################
### Step 2: REACT ###
#####################

### Define variables
outdir=
atlasdir=
atlas=

cd $outdir

### Create masks ###

react_masks subject_list.txt "$atlasdir"/PETatlas_"$atlas".nii.gz "$atlasdir"/react-fmri-main/data/gm_mask.nii.gz out_masks_"$atlas"


### Run react ###

for i in sub-6??; do

react "$outdir"/"$i"/filtered_func_data.nii.gz "$outdir"/out_masks_"$atlas"/mask_stage1.nii.gz "$outdir"/out_masks_"$atlas"/mask_stage2.nii.gz $atlasdir/PETatlas_"$atlas".nii.gz "$outdir"/REACT_"$atlas"/"$i"

done

#########################
### Step 3: Randomise ###
#########################

### Define variables
atlas=5HT1A_SERT
outdir=
designdir=
niftidir=
atlasdir=
reactdir=

### Prep randomise

cd $reactdir

## All conditions
# Create lists
find $PWD -name '*IC0.nii.gz' | sort > $designdir/5HT1A.txt
find $PWD -name '*IC1.nii.gz' | sort > $designdir/SERT.txt


## Subset conditions
for condition in plac low high plac_low plac_high low_high; do
	for i in `cat $designdir/"$condition"_numonly.txt`; do echo "$reactdir"/sub-"$i"_react_stage2_IC0.nii.gz; done > $designdir/5HT1A_"$condition".txt
	for i in `cat $designdir/"$condition"_numonly.txt`; do echo "$reactdir"/sub-"$i"_react_stage2_IC1.nii.gz; done > $designdir/SERT_"$condition".txt
done

## Create 4D files

cd $designdir

for map in 5HT1A SERT; do
fslmerge -t $niftidir/4D_"$map" `cat "$map".txt`
fslmaths $niftidir/4D_"$map" -mul -1 $niftidir/4D_"$map"_neg

	for condition in plac low high ; do
	fslmerge -t $niftidir/4D_"$map"_"$condition" `cat "$map"_"$condition".txt`
	fslmaths $niftidir/4D_"$map"_"$condition" -mul -1 $niftidir/4D_"$map"_"$condition"_neg
	for condition in plac_low plac_high low_high; do
	fslmerge -t $niftidir/4D_"$map"_"$condition" `cat "$map"_"$condition".txt`
	done
done

### 3.1 Main effects ###

## Define mask
mask=

## 3.1.1 All subs
for map in 5HT1A SERT; do
randomise -i $niftidir/4D_"$map" -o $outdir/OneSamp_"$map" -m $mask -1 -v 5 -T -c 2.3
randomise -i $niftidir/4D_"$map"_neg -o $outdir/OneSamp_"$map"_neg -m $mask -1 -v 5 -T -c 2.3

## relation with SERT post2 values
randomise -i $niftidir/4D_"$map".nii.gz -o $outdir/corr_SPECT2thal_"$map" -d $designdir/SPECT_BPpost_thal_no618.mat -t $designdir/correlation.con -m $mask -v 5 -T -c 2.3

	## 3.1.2 Main effect per group
	for condition in plac low high; do 
	randomise -i $niftidir/4D_"$map"_"$condition" -o $outdir/OneSamp_"$map"_"$condition" -m $mask -1 -v 5 -T -c 2.3
	randomise -i $niftidir/4D_"$map"_"$condition"_neg -o $outdir/OneSamp_"$map"_"$condition"_neg -m $mask -1 -v 5 -T -c 2.3
	done

	### 3.2 Dose-dependent effect ###
	for condition in plac_low plac_high low_high; do
	randomise -i $niftidir/4D_"$map"_"$condition".nii.gz -o $outdir/TwoSamp_"$map"_"$condition" -d $designdir/design_groups_"$condition".mat -t $designdir/design_two_groups.con -m $mask -v 5 -T
	randomise -i $niftidir/4D_"$map"_"$condition"_neg.nii.gz -o $outdir/TwoSamp_"$map"_"$condition"_neg -d $designdir/design_groups_"$condition".mat -t $designdir/design_two_groups.con -m $mask -v 5 -T
	done

done

#########################################################
### Step 4: Extracting location and individual values ###
#########################################################

### Define variables
atlas=5HT1A_SERT
outdir=
designdir=
niftidir=
atlasdir=
reactdir=

#	4.2 Main effects per group
#	4.3 Dose-dependent effects
### 4.1 Main effects ###
# Convert to MNI
	for map in 5HT1A SERT; do
	flirt -in OneSamp_"$map"_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out OneSamp_"$map"_tfce_corrp_tstat1_flirt.nii.gz
	flirt -in OneSamp_"$map"_neg_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out OneSamp_"$map"_neg_tfce_corrp_tstat1_flirt.nii.gz
	# Extract areas
autoaq -i OneSamp_"$map"_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o  Autoaq_OneSamp_"$map"_tfce_corrp_tstat1_output_cortical.txt
autoaq -i OneSamp_"$map"_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o OneSamp_"$map"_tfce_corrp_tstat1_output_subcortical.txt

autoaq -i OneSamp_"$map"_neg_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o Autoaq_OneSamp_"$map"_neg_tfce_corrp_tstat1_output_cortical.txt
autoaq -i OneSamp_"$map"_neg_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o Autoaq_OneSamp_"$map"_neg_tfce_corrp_tstat1_output_subcortical.txt

	# Convert to MNI
	flirt -in corr_SPECT2thal_"$map"_tfce_corrp_tstat2.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out corr_SPECT2thal_"$map"_tfce_corrp_tstat2_flirt.nii.gz
	# Extract areas
	autoaq -i corr_SPECT2thal_"$map"_tfce_corrp_tstat2_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o TwoSamp_"$map"_"$condition"_tfce_corrp_tstat2_output_cortical.txt
	autoaq -i corr_SPECT2thal_"$map"_tfce_corrp_tstat2_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o TwoSamp_"$map"_"$condition"_tfce_corrp_tstat2_output_subcortical.txt

### 4.2 Main effects per group ###
for condition in plac low high; do
	# Convert to MNI
	flirt -in OneSamp_"$map"_"$condition"_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out OneSamp_"$map"_"$condition"_tfce_corrp_tstat1_flirt.nii.gz
	flirt -in OneSamp_"$map"_"$condition"_neg_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out OneSamp_"$map"_"$condition"_neg_tfce_corrp_tstat1_flirt.nii.gz
	# Extract areas
	autoaq -i OneSamp_"$map"_"$condition"_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o  OneSamp_"$map"_"$condition"_tfce_corrp_tstat1_output_cortical.txt
	autoaq -i $OneSamp_"$map"_"$condition"_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o OneSamp_"$map"_"$condition"_tfce_corrp_tstat1_output_subcortical.txt

	autoaq -i OneSamp_"$map"_"$condition"_neg_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o OneSamp_"$map"_"$condition"_neg_tfce_corrp_tstat1_output_cortical.txt
	autoaq -i OneSamp_"$map"_"$condition"_neg_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o OneSamp_"$map"_"$condition"_neg_tfce_corrp_tstat1_output_subcortical.txt
	done

### 4.3 Dose-dependent effects ###
for condition in low_high; do # Only significant difference
	# Convert to MNI
	flirt -in TwoSamp_"$map"_"$condition"_tfce_corrp_tstat2.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out TwoSamp_"$map"_"$reg"_tfce_corrp_tstat2_flirt.nii.gz
	# Extract areas
	autoaq -i TwoSamp_"$map"_"$condition"_tfce_corrp_tstat2_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.9875 -o TwoSamp_"$map"_"$condition"_tfce_corrp_tstat2_output_cortical.txt
	autoaq -i TwoSamp_"$map"_"$condition"_tfce_corrp_tstat2_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.9875 -o TwoSamp_"$map"_"$condition"_tfce_corrp_tstat2_output_subcortical.txt
	done
done


##################################################################################################
##################################################################################################


#############################
#### task-fMRI REACT PPI ####
#############################

# Steps in the script:
# 1. Post-fMRIprep processing
#	1.1 Regressing out CSF and WM
#	1.2 Band pass filtering
# 2. PPI in FSL Feat
#	2.1 Extract timeseries from SERT and 5HT1A density maps
#	2.2 Run PPI in Feat
#	2.3 FSL fake registration
# 3. Randomise
# 	3.1 Main effects
#		3.1.1 Main effects all subs
#		3.1.2 Main effects per group
#		3.1.3 Relation with thalamic binding
# 	3.2 Dose-dependent effects with small-volume correction
#		3.2.1 Dose-dependent effects with volume determined on placebo group
# 4. Extracting cluster location and individual values
#	4.1 Main effects
#	4.2 Main effects per group
#	4.3 Dose-dependent small volume correction
#	4.4 Dose-dependent small volume correction placebo
#
#
##################################################################################################

USER=
type=PPI

####################################
### Step 1: Post-fMRI processing ###
####################################

#########################
### Step 2: REACT PPI ###
#########################

maindir=
scriptdir=

cd $maindir
for sub in `ls -d sub-???`; do

	# find and check if filename is present
	feat_dir=$maindir/$sub/
	if [ -d $feat_dir ]; then
	
	### 2.1 Extract timeseries from 5HT1A and SERT density maps ###
	fsl_glm -i $maindir/"$sub"/filtered_func_data.nii.gz -d $scriptdir/masks/masks_original/PETatlas_5HT1A_SERT.nii.gz -o $maindir/"$sub"/"$sub"_task-faces_run-1_tc_5HT1A_SERT.txt --demean -m $scriptdir/masks/out_masks_5HT1A_SERT/mask_stage1.nii.gz 

	cat $maindir/"$sub"/"$sub"_task-faces_run-1_tc_5HT1A_SERT.txt > $maindir/"$sub"/"$sub"_task-faces_run-1_tc_5HT1A_SERT.csv
	tc_all=` cat $maindir/"$sub"/"$sub"_task-faces_run-1_tc_5HT1A_SERT.csv`
	echo $tc_all | awk -F '  ' '{print $1}' >> $maindir/"$sub"/"$sub"_task-faces_run-1_tc_all_5HT1A.txt
	echo $tc_all | awk -F '  ' '{print $3}' >> $maindir/"$sub"/"$sub"_task-faces_run-1_tc_all_SERT.txt
	
	### 2.2 Run PPI in Feat ###
	cp $scriptdir/Run_PPI_REACT/PPI_REACT_5HT1A_SERT.fsf $maindir/"$sub"/PPI_REACT_5HT1A_SERT_"$sub"_task-faces_run-1.fsf
	sed -i -e 's/CHANGEPP/'$sub'/g' $maindir/"$sub"/PPI_REACT_5HT1A_SERT_"$sub"_task-faces_run-1.fsf
	sed -i -e 's/CHANGEVALUE/'0'/' $maindir/"$sub"/PPI_REACT_5HT1A_SERT_"$sub"_task-faces_run-1.fsf
	feat $maindir/"$sub"/PPI_REACT_5HT1A_SERT_"$sub"_task-faces_run-1.fsf
	
	fi
done


### 2.3 FSL fake registration ###

# See ReadMe file

#########################
### Step 3: Randomise ###
#########################

### Define variables
atlas=5HT1A_SERT
outdir=
designdir=
niftidir=
atlasdir=
ppidir=

### Prep randomise
## All conditions
# Create lists
find $PWD -name 'cope1.nii.gz' -path '*/sub*_faces_run-1_PPI_5HT1A_SERT.feat/*'| sort > $designdir/_task.txt
find $PWD -name 'cope2.nii.gz' -path '*/sub*_faces_run-1_PPI_5HT1A_SERT.feat/*'| sort > $designdir/_5HT1A_phys.txt
find $PWD -name 'cope3.nii.gz' -path '*/sub*_faces_run-1_PPI_5HT1A_SERT.feat/*'| sort > $designdir/_SERT_phys.txt
find $PWD -name 'cope4.nii.gz' -path '*/sub*_faces_run-1_PPI_5HT1A_SERT.feat/*'| sort > $designdir/_5HT1A_PPI.txt
find $PWD -name 'cope5.nii.gz' -path '*/sub*_faces_run-1_PPI_5HT1A_SERT.feat/*'| sort > $designdir/_SERT_PPI.txt

## Subset conditions
for condition in plac low high plac_low plac_high low_high; do
	for i in `cat $designdir/"$condition"_numonly.txt`; do echo "$ppidir"/sub-"$i"/sub-"$i"_faces_run-1_PPI_5HT1A_SERT.feat/stats/cope1.nii.gz; done > $designdir/task_"$condition".txt
	for i in `cat $designdir/"$condition"_numonly.txt`; do echo "$ppidir"/sub-"$i"/sub-"$i"_faces_run-1_PPI_5HT1A_SERT.feat/stats/cope2.nii.gz; done > $designdir/5HT1A_phys_"$condition".txt
	for i in `cat $designdir/"$condition"_numonly.txt`; do echo "$ppidir"/sub-"$i"/sub-"$i"_faces_run-1_PPI_5HT1A_SERT.feat/stats/cope3.nii.gz; done > $designdir/SERT_phys_"$condition".txt
	for i in `cat $designdir/"$condition"_numonly.txt`; do echo "$ppidir"/sub-"$i"/sub-"$i"_faces_run-1_PPI_5HT1A_SERT.feat/stats/cope4.nii.gz; done > $designdir/5HT1A_PPI_"$condition".txt
	for i in `cat $designdir/"$condition"_numonly.txt`; do echo "$ppidir"/sub-"$i"/sub-"$i"_faces_run-1_PPI_5HT1A_SERT.feat/stats/cope5.nii.gz; done > $designdir/SERT_PPI_"$condition".txt
done

## Create files
for map in task 5HT1A_phys SERT_phys 5HT1A_PPI SERT_PPI; do
	fslmerge -t 4D_"$type" `cat $designdir/_"$type".txt`
	fslmaths 4D_"$type" -mul -1 4D_"$type"_neg
	for condition in plac low high; do
		fslmerge -t 4D_"$map"_"$condition"
		fslmaths 4D_"$map"_"$condition" -mul -1 4D_"$map"_"$condition"_neg
	done
		for condition in plac_low plac_high low_high; do
		fslmerge -t 4D_"$map"_"$condition"
	done
done

## Define mask
mask=

### 3.1 Main effects ###
## 3.1.1 All subs
for map in task 5HT1A_phys SERT_phys 5HT1A_PPI SERT_PPI; do
	randomise -i $niftidir/4D_"$map" -o $outdir/OneSamp_"$map" -m $mask -1 -v 5 -T -c 2.3
	randomise -i $niftidir/4D_"$map"_neg -o $outdir/OneSamp_"$map"_neg -m $mask -1 -v 5 -T -c 2.3

	## relation with SERT post2 values
	randomise -i $niftidir/4D_"$map".nii.gz -o $outdir/corr_SPECT2thal_"$map" -d $designdir/SPECT_BPpost_thal_no618.mat -t $designdir/correlation.con -m $mask -v 5 -T -c 2.3
	## 3.1.2 Main effect per group
	for condition in plac low high; do 
	randomise -i $niftidir/4D_"$map"_"$condition" -o $outdir/OneSamp_"$map"_"$condition" -m $mask -1 -v 5 -T -c 2.3
	randomise -i $niftidir/4D_"$map"_"$condition"_neg -o $outdir/OneSamp_"$map"_"$condition"_neg -m $mask -1 -v 5 -T -c 2.3
	done

	### 3.2 Dose-dependent effect ###
	for condition in plac_low plac_high low_high; do
	randomise -i $niftidir/4D_"$map"_"$condition".nii.gz -o $outdir/TwoSamp_"$map"_"$condition" -d $designdir/design_groups_"$condition".mat -t $designdir/design_two_groups.con -m $mask -v 5 -T
	done

done

#########################################################
### Step 4: Extracting location and individual values ###
#########################################################

### Set variables
maindir=
dir_4D=

cd $outdir


### 4.1 Main effects ###
cd $maineffdir

## Extract cluster location
for map in task 5HT1A_phys SERT_phys 5HT1A_PPI SERT_PPI; do
	# Convert to MNI
	flirt -in $maineffdir/OneSamp_"$map"_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out $maineffdir/OneSamp_"$map"_"$reg"_tfce_corrp_tstat1_flirt.nii.gz
	flirt -in $maineffdir/OneSamp_"$map"_neg_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out $maineffdir/OneSamp_"$map"_"$reg"_neg_tfce_corrp_tstat1_flirt.nii.gz
	
	# Extract areas
	autoaq -i $maineffdir/OneSamp_"$map"_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o  $maineffdir/Autoaq_OneSamp_"$map"_"$reg"_tfce_corrp_tstat1_output_cortical.txt
	autoaq -i $maineffdir/OneSamp_"$map"_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o $maineffdir/Autoaq_OneSamp_"$map"_"$reg"_tfce_corrp_tstat1_output_subcortical.txt

	autoaq -i $maineffdir/OneSamp_"$map"_neg_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o $maineffdir/Autoaq_OneSamp_"$map"_"$reg"_neg_tfce_corrp_tstat1_output_cortical.txt
	autoaq -i $maineffdir/OneSamp_"$map"_neg_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o $maineffdir/Autoaq_OneSamp_"$map"_"$reg"_neg_tfce_corrp_tstat1_output_subcortical.txt
done

## Correlation SPECT values
# Convert to MNI
flirt -in corr_SPECT2thal_5HT1A_PPI_tfce_corrp_tstat2.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out corr_SPECT2thal_5HT1A_PPI_tfce_corrp_tstat2_flirt.nii.gz

# Extract areas
autoaq -i corr_SPECT2thal_5HT1A_PPI_tfce_corrp_tstat2_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o  corr_SPECT2thal_5HT1A_PPI_tfce_corrp_tstat2_output_cortical.txt

autoaq -i corr_SPECT2thal_5HT1A_PPI_tfce_corrp_tstat2_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o corr_SPECT2thal_5HT1A_PPI_tfce_corrp_tstat2_output_subcortical.txt

# Extract individual values
# Create thresholded and binarised masks
fslmaths corr_SPECT2thal_5HT1A_PPI_tfce_corrp_tstat2.nii.gz -thr 0.95 -bin corr_SPECT2thal_5HT1A_PPI_tfce_corrp_tstat2_sign_bin

# Extract values from 4D file
fslstats -t $dir_4D/4D_5HT1A_PPI.nii.gz -k $corr_SPECT2thal_5HT1A_PPI_sign_bin.nii.gz -M > corr_SPECT2thal_5HT1A_PPI_individualvals.txt

### Main effects per group ###
cd $maingroupdir

## 4.2.1 Extract cluster location
for map in task 5HT1A_phys SERT_phys 5HT1A_PPI SERT_PPI; do
	for condition in plac low high; do
	# Convert to MNI
	flirt -in $maingroupdir/OneSamp_"$map"_"$condition"_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out $maingroupdir/OneSamp_"$map"_"$reg"_"$condition"_tfce_corrp_tstat1_flirt.nii.gz
	flirt -in $maingroupdir/OneSamp_"$map"_"$condition"_neg_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out $maingroupdir/OneSamp_"$map"_"$reg"_"$condition"_neg_tfce_corrp_tstat1_flirt.nii.gz
	# Extract areas
	autoaq -i $maingroupdir/OneSamp_"$map"_"$reg"_"$condition"_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o  $maingroupdir/Autoaq_OneSamp_"$map"_"$reg"_"$condition"_tfce_corrp_tstat1_output_cortical.txt
		autoaq -i $maingroupdir/OneSamp_"$map"_"$reg"_"$condition"_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o $maingroupdir/Autoaq_OneSamp_"$map"_"$reg"_"$condition"_tfce_corrp_tstat1_output_subcortical.txt

	autoaq -i $maingroupdir/OneSamp_"$map"_"$reg"_"$condition"_neg_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o $maingroupdir/Autoaq_OneSamp_"$map"_"$reg"_"$condition"_neg_tfce_corrp_tstat1_output_cortical.txt
	autoaq -i $maingroupdir/OneSamp_"$map"_"$reg"_"$condition"_neg_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o $maingroupdir/Autoaq_OneSamp_"$map"_"$reg"_"$condition"_neg_tfce_corrp_tstat1_output_subcortical.txt
	done
done


### 4.3 Dose-dependent effects small-volume ###
### Significant effects: 5HT1A PPI low vs high

cd $smallvoldir

## Extract cluster location
# Convert to MNI
flirt -in $smallvoldir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out $smallvoldir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_flirt.nii.gz

# Extract areas
autoaq -i $smallvoldir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.9875 -o $smallvoldir/Autoaq_group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_output_cortical.txt

autoaq -i $maingroupdir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.9875 -o $smallvoldir/Autoaq_group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_output_subcortical.txt

## Extract individual vaues
# Create thresholded and binarised masks
fslmaths $smallvoldir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1.nii.gz -thr 0.9875 -bin $smallvoldir/group_PPI_plac_low_5HT1A_tfce_corrp_sign_bin

# Extract values from 4D file
fslstats -t $dir_4D/4D_5HT1A_PPI_plac_low.nii.gz -k $smallvoldir/OneSamp_5HT1A_interaction_tfce_corrp_tstat1_sign_bin.nii.gz -M > 5HT1A_PPI_plac_low_sign_individualvals.txt


### 4.4 Small-volume main placebo effect ###

cd $smallvolplacdir

## 4.4.1 Extract cluster location
# Convert to MNI
flirt -in $smallvolplacdir/group_PPI_plac_high_5HT1A_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out $smallvolplacdir/group_PPI_plac_high_5HT1A_tfce_corrp_tstat1_flirt.nii.gz

flirt -in $smallvolplacdir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1.nii.gz -ref MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out $smallvolplacdir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_flirt.nii.gz

flirt -in $smallvolplacdir/group_diff_5HT1A_groups_tfce_corrp_tstat7.nii.gz -ref /MNI152_T1_2mm.nii.gz -applyxfm -usesqform -out $smallvolplacdir/group_diff_5HT1A_groups_tfce_corrp_tstat7_flirt.nii.gz

# Extract areas
autoaq -i $smallvolplacdir/group_PPI_plac_high_5HT1A_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o $smallvolplacdir/Autoaq_group_PPI_plac_high_5HT1A_tfce_corrp_tstat1_output_cortical.txt
autoaq -i $smallvolplacdir/group_PPI_plac_high_5HT1A_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o $smallvolplacdir/Autoaq_group_PPI_plac_high_5HT1A_tfce_corrp_tstat1_output_subcortical.txt

autoaq -i $smallvolplacdir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o $smallvolplacdir/Autoaq_group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_output_cortical.txt
autoaq -i $smallvolplacdir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o $smallvolplacdir/Autoaq_group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_output_subcortical.txt

autoaq -i $smallvolplacdir/group_diff_5HT1A_groups_tfce_corrp_tstat7_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.95 -o $smallvolplacdir/Autoaq_group_diff_5HT1A_groups_tfce_corrp_tstat7_output_cortical.txt
autoaq -i $smallvolplacdir/group_diff_5HT1A_groups_tfce_corrp_tstat7_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.95 -o $smallvolplacdir/Autoaq_group_diff_5HT1A_groups_tfce_corrp_tstat7_output_subcortical.txt

## 4.4.2 Extract individual values
# Create thresholded and binarised masks
fslmaths $smallvolplacdir/group_PPI_plac_high_5HT1A_tfce_corrp_tstat1.nii.gz -thr 0.95 -bin $smallvolplacdir/group_PPI_plac_high_5HT1A_tfce_corrp_tstat1_sign_bin

fslmaths $smallvolplacdir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1.nii.gz -thr 0.95 -bin $smallvolplacdir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_sign_bin

fslmaths $smallvolplacdir/group_diff_5HT1A_groups_tfce_corrp_tstat7.nii.gz -thr 0.95 -bin $smallvolplacdir/group_diff_5HT1A_groups_tfce_corrp_tstat7_sign_bin

# Extract values from 4D file
fslstats -t $dir_4D/4D_5HT1A_PPI_plac_high.nii.gz -k $smallvolplacdir/group_PPI_plac_high_5HT1A_tfce_corrp_tstat1_sign_bin.nii.gz -M > $smallvolplacdir/5HT1A_PPI_plac_high_sign_individualvals.txt

fslstats -t $dir_4D/4D_5HT1A_PPI_plac_low.nii.gz -k $smallvolplacdir/group_PPI_plac_low_5HT1A_tfce_corrp_tstat1_sign_bin.nii.gz -M > $smallvolplacdir/5HT1A_PPI_plac_low_sign_individualvals.txt

fslstats -t $dir_4D/4D_5HT1A_PPI.nii.gz -k $smallvolplacdir/group_diff_5HT1A_groups_tfce_corrp_tstat7_sign_bin.nii.gz -M > $smallvolplacdir/5HT1A_PPI_plac_low_high_sign_individualvals.txt
