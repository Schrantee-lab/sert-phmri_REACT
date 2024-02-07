#!/bin/bash
#
# Script to generate required 'registered' versions of image data for
# FSL/FEAT so it will perform group analysis. By default this package
# performs a registration to the standard space AFTER the first level
# analysis. However, we have already done this (and it would not work
# well with rCBV data). Hence, we fool it into thinking it has done
# so itself by generating the necessary files.
# Run this script following first level FEAT analysis.
#
# Arguments: list of .feat directories
#
# Calling sequence: fsl_gpAnal_prep.sh feat_dir1 [feat_dir2 [...]]
#
# Example use:
#  fsl_fake_reg.sh *.feat
#
# Modification history:
#
# 2013-12-10  pg   Removed dependency on unit.mat by creating the file HERE.
#                  Included test to skip existing reg-directory.
#                  Included early-abort flag to stop script after first error.
#                  Removed requirement to run from data directory
#                  Replaced hard-coded FSL folder with predefined env. variable.
#                  Replaced passing environment variable name with proper dynamic argument list.
#                  Made it whitespace-proof
# 2011-03-08  ag   Update to FSL v4 and MRI-1 workstation; IIT@NEST
#                  For the sciprt to work, make sure the default FSL output type is NIFTI
#                  (Default is NIFTI_GZ). This can be done changing fls.sh script in usr/share/fsl/4.1/etc/
# 2005-08-11  ajs  Update to template v.3
# 2005-05-19  ajs  Original
# -

# stop after any error:
set -e

if [ $# -eq 0 ];
then
	MYNAME=`basename "$0"`
	echo "Usage:"
	echo " ${MYNAME} feat_dir1 [feat_dir2 [...]]"
	echo
	echo "Example:"
	echo " ${MYNAME} *_2.feat"
	echo
fi

for DD in "$@";
do
	if [ ! -e "${DD}/example_func.nii.gz" ]; then
		echo "Missing ${DD}/example_func.nii.gz; skipping this directory"
	elif [ -d "${DD}/reg" ]; then
		#echo "Directory ${DD}/reg exists; skipping this directory"
		cd "${DD}/reg"

	else
		echo "Creating unity registration in $DD..."
		mkdir "${DD}/reg"
		OLDDIR=$PWD
	  	cd "${DD}/reg"

		# create the 3 required transformation matrices (using unit transf.)
		for MM in standard2example_func.mat highres2standard.mat example_func2standard.mat;	do
			cat > "${MM}" <<EOF
1  0  0  0
0  1  0  0
0  0  1  0
0  0  0  1
EOF
		done

	#	ln -s ../example_func.nii.gz example_func.nii.gz
	#	ln -s "/scratch/opt/fsl-6.0.0/data/standard/MNI152_T1_2mm_brain.nii.gz" standard.nii.gz
	#	ln -s "/scratch/opt/fsl-6.0.0/data/standard/MNI152_T1_2mm_brain.nii.gz" highres.nii.gz
		cp  ../example_func.nii.gz example_func.nii.gz
		cp  "/opt/amc/fsl-6.0.0/data/standard/MNI152_T1_2mm_brain.nii.gz" standard.nii.gz
		cp  "/opt/amc/fsl-6.0.0/data/standard/MNI152_T1_2mm_brain.nii.gz" highres.nii.gz
		cd "$OLDDIR"
	fi
done

echo "done"
exit 0
