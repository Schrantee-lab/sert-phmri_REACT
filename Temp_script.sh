#!/bin/bash
#SBATCH --job-name=slurm_randomise_mc
#SBATCH --mem=8G               # max memory per node
#SBATCH --cpus-per-task=2      # max CPU cores per MPI process
#SBATCH --time=0-06:58       # time limit (DD-HH:MM)
#SBATCH --partition=rng-short  # rng-short is default, but use rng-long if time exceeds 7h
# ---mail to yourself if job has ended or failed (may end up in SPAM)
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=d.e.boucherie@amsterdamumc.nl


################################################################################

outdir=/home/deboucherie/lood_storage/divi/Projects/sert-phmri/NIFTI/BIDS/output/REACT-PPI/randomise_allmasks/Small_volume_correction/significant
niftidir=/home/deboucherie/lood_storage/divi/Projects/sert-phmri/NIFTI/BIDS/output/REACT-PPI/randomise_allmasks/4D_allmasks/
designdir=/home/deboucherie/lood_storage/divi/Projects/sert-phmri/NIFTI/BIDS/output/REACT-PPI/designs/group_files/daphne/


autoaq -i $outdir/group_allmasks_PPI_plac_low_5HT1A_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Cortical Structural Atlas" -t 0.9843 -o $outdir/Autoaq_ggroup_allmasks_PPI_plac_low_5HT1A_output_cortical.txt
autoaq -i $outdir/group_allmasks_PPI_plac_low_5HT1A_tfce_corrp_tstat1_flirt.nii.gz -a "Harvard-Oxford Subcortical Structural Atlas" -t 0.9843 -o $outdir/Autoaq_group_allmasks_PPI_plac_low_5HT1A_output_subcortical.txt

