# QCMD_analyze
This repository includes a Matlab script, a README file, and a copy of a working data file. The Matlab script was created with the intent to calculate Sauerbrey mass and thin film compliance from QCM-D data measurements. The script is still undergoing changes to improve the amount of work required from the user. This script was developed by Katelyn Lunny using MATLAB R2023b. This script should work for older versions of Matlab as well, but is not a guarantee. 

# Using QCM_analyze.m
## Setup
- To use this Matlab script, please download the attached m file and run it first with the given csv file to ensure the script is able to run
-  The script has two main outputs:
  - SauerbreyMass.csv 
  - ThinFilmCompliance.csv
- The values the script should produce when using the given csv (07092024_Fn at 50 ug per ml at pH 4 + dSF pH4 then pH 7 PBS wash on func. Au at 25C_n=4_KL.csv) can be found in the SauerbreyMassTest.csv and ThinFilmComplianceTest.csv attached to this repository
- If the numbers do not match, please check the following
  - [ ] the issue and liscensing of the Matlab installed on the computer. This script was developed using MATLAB R2023b, but should work for older verisons of Matlab, but is not a guarantee. Please try updating Matlab to resolve the issue.

## Use
- Once you have confirmed that the output from running the m file with the provided csv file matches the expected (found in the SauerbreyMassTest.csv and ThinFilmComplianceTest.csv), you may then proceed with using QCM-D_analyze with your own data.
- Checklist for use:
  - [ ] File of interest is in csv format
  - [ ] File of interest includes all data from QCM-D ('Add all' from QCM-D export)
  - [ ] Time stamps for each step in the experiment

- In its current state, the script is unable to automatically detect step changes in the experiment.
  -To adjust the ranges for each experimental step, please adjust the ranges in lines 72-79 (USER INPUT section)
    - This ranges correspond to the row numbers from the data array for the experimental step.
    - The original baseline must include the entire baseline (i.e. 0min to 10min)
    - All other steps must be the last 2-3 minutes of the step
    - To confirm you have selected the correct data points, please refer to the 'Frequency Against Time with Overlay of Time Stamps for Steps' graph, which should plot the selected ranges in a different color/line style overlayed with the frequency vs. time
  
- For any issues experienced in running QCM_analyze.m, please check for:
  - [ ] File of interest is in csv format
  - [ ] File of interest includes all data from QCM-D ('Add all' from QCM-D export)
  - [ ] Time stamps for each step in the experiment
  - [ ] any missing data (specifically in the dissipation changes). At its current state, any missing data from dissipation changes can result in erroneous Jf prime values. If you discover any missing dissipation changes, please remove the corresponding del_bandcondition or ndel_condition from the Jf prime calculation.


