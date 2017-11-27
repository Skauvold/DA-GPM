Matlab scripts with supporting functions and input data used to carry out computations, analyse results and create figures in connection with the article "Data Assimilation for a Geological Process Model Using the Ensemble Kalman Filter" by J. Skauvold and J. Eidsvik. Original manuscript submitted to Basin Research April 2017. Revised manuscript submitted October 2017.
https://arxiv.org/abs/1711.07763

List of files included in this repository:

Scripts for performing computations and/or analysing results (see Note 2 below):
- ANSCaseEnKF_bw.m (Script for running EnKF on North Slope, Alaska case)
- enkftrial.m (Script for running a synthetic data trial using the EnKF)
- enstrial.m (Script for running a synthetic data trial using the ensemble smoother)
- trialsummary.m (Script for computing statistics from EnKF and EnS trial results on synthetic case)
- trialsummaryMDA.m (Script for computing statistics from MDA trial results on synthetic case)

Scripts for making figures (see also Note 3 below):
- plot2Dsection.m (Script for creating Fig. 1)
- makeSynthCaseFigure.m (Script for creating Fig. 2)
- makeSLSSfigure.m (Script for creating Fig. 4)
- plotresults_blockwise.m (Script for creating Fig. 7, Fig. 8, and some other plots)
- makeANSparamFigure.m (Script for creating Fig. 9)

Supporting functions:
- mkreffun.m (Function for generating a reference realisation. Called in enkftrial.m)
- crps.m (Function for computing CRPS, called in trialsummary.m and trialsummaryMDA.m)
- eta2p.m (Function for converting Gaussian variables into proportions)
- gr_obs_model.m (Observation operator for North Slope, Alaska case)

Data files:
- myref.mat (State variables of the reference realisation shown in Fig. 2)
- LUTSSLL_{EnKF,EnS,MDA}.mat (Data plotted in Fig. 4)
- dzArrayPrior.mat (Ensemble of unconditional realisations of thickness for North Slope case)
- run1_arrays.mat (Posterior GR and thickness realisations for North Slope case)
- thetaPriorRealisations.mat (Samples from prior distribution of parameters)
- Tunalik1data.mat (Low-res depth and gamma ray data from Tunalik 1 well log)
- priorResults.mat (Unconditional thickness and gamma ray realisations for North Slope case)
- Tunalik1grlog.mat (High-res depth and gamma ray data from Tunalik 1 well log)
- grArray.mat (Unconditional gamma ray realisations for North Slope case)
- thetapost3km (Posterior realisations of parameters for North Slope case with very uncertain thickness observations)
- thetapost30m (Posterior realisations of parameters for North Slope case with moderately uncertain thickness observations)

Note 1: The file ANSCaseEnKF_output.mat, which contains the analysis ensemble from an EnKF run on the North Slope, Alaska case is required for running the plot2Dsection.m and plotresults_blockwise.m scripts. Unfortunately, it is too large (~1.27 GB) to include in the repository. The file is available at http://folk.ntnu.no/skauvold/GPM-EnKF/.

Note 2: These scripts are not intended to be run as they are. Some of them depend on the simulator (GPM), and all of them depend on supporting functions not included in this repository. For confidentiality reasons, some lines in these scripts have been replaced by comments.

Note 3: These scripts can be run as they are to produce all figures in the paper, except figures 3, 5 and 10. Fig. 3 was made using the plot editor GUI in Matlab. The map in Fig. 5 was adapted from the original figure using the Inkscape vector graphics editor. Fig. 10 was created using the TikZ package in LaTeX. Two of the scripts require a file which is not included in this repository, but can be obtained elsewhere (see Note 1).

Note 4: The code (although not necessarily the files included in this repository) makes use of inpaint_nans and jsonlab from the Matlab Central file exchange.


Jacob Skauvold
14 Oct 2017

For more information, contact jacob dot skauvold at ntnu dot no.
