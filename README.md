# sert-phmri_REACT
This folder contains all code associated with the manuscript of Boucherie et al. (2023), which was published in Journal of Psychopharmacology. 
You can find this manuscript here: https://journals.sagepub.com/doi/epub/10.1177/02698811231211154
It was previously published as a preprint, which can be found here: https://www.biorxiv.org/content/10.1101/2022.10.20.512503v3

# Main batch script and file descriptions
The file Main_script_REACT_sert-phMRI.sh contains all relevant steps for the analysis. This script is divided into two parts:
- Part 1: resting-state fMRI REACT analysis
- Part 2: task-based fMRI PPI analysis

In the script, several other files are called. 
For the rs-fMRI analysis, this is the script used for postprocessing (after fMRIprep): Postprocessing_rsfmri.sh

For the tb-fMRI analysis, these are the scripts used for:
- postprocessing (Postprocessing_tbfmri.sh);
- First-level PPI analysis design in FSL Feat (PPI_REACT_5HT1A_SERT.fsf);
- fake registration after the first-level PPI analysis in Feat (fsl_fake_reg.sh) (https://www.youtube.com/watch?v=U3tG7JMEf7M)
