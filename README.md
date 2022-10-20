# Lifespan_MEG

Analysis of age effects on MEG-derived markers and cortical thickness in an adult lifespan sample 

A pre-print describing the analysis and the results can be accessed through XX.

The provided repository contains all relevant scripts and main results but requires some external tools that need to be either placed in the Matlab path, included in the scripts/utilities directory, or in the R path.

External tools

- Matlab (required, tested with version MATLAB 9.5.0.1298439 (R2018b))
- Fieldtrip (required, tested with version 20191127) https://www.fieldtriptoolbox.org/
- SPM12 (required, tested with version 7487, https://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
- Freesurfer (required, tested with version 6.0.0), https://surfer.nmr.mgh.harvard.edu/
- PALM FSL (required), https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/PALM
- R (required, tested with version 3.6.3 (2020-02-29)) plus R packages as indicated in the /R folder

Data

- Raw data can be obtained via: https://camcan-archive.mrc-cbu.cam.ac.uk/dataaccess/ upon request

Analysis steps

1. Pre-processing of single subject MEG data using the script 'cs_process_subject_neuromag.m' 
2. Export processed MEG and cortical thickness measures for statistical analyses (vertex-level and globally) using the scripts 'cs_export_metrics_all_subj.m' and 'cs_export_thickness.m'
3. Save demographics of the individuals included using 'cs_save_demographics.m'
4. For statistical analyses as described in the pre-print and plotting the figures use 'cs_statistics_and_visualization.m' 
