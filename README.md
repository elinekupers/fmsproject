# fmsproject

Welcome to the MATLAB based code repository of the MEG Forward Models and neural Synchrony (FMS) project

## Goal
Our goal in this project is to get a better understanding of the underlying neural synchrony 
of the brain responses that produce a particular MEG scalp topography. 

Inferring the underlying synchrony at the scalp level is challenging, because the responses
that your measure (for example with MEG), are not only affected by synchrony, but also by 
changes in neuronal amplitude and cancellation or summing of cortical sources over space.
Disentangling these factors is not easy.

In this project, we make predictions about the projection of the sources to the sensors.
We simulate neural signals in V1-V3 on the cortical surface, with either the same phase, 
randomized phases, or something in between.

We combine these simulated neural responses with an MEG forward model, to predict the MEG
responses on the sensors (axial gradiometers), and compare them to measured MEG data from
individuals viewing a full-field high contrast-reversing checkerboard stimulus.


## Repository overview

- external 			: functions from other repositories
- figurescripts		: contains master script s_makeAllFMSFigures.m, calls separate makeFigureX functions, where X is figure number in manuscript
	- subfunctions 	: folder with additional functions needed to run figure scripts
		- makeFigure2 - Visualize visually evoked steady state and broadband MEG sensor responses
		- makeFigure3 - Visualize MEG sensor data on topographic maps for SSVEF and broadband.
		- makeFigure5 - Visualize model predictions simulating the MEG sensor responses for coherent and incoherent neural sources.
		- makeFigure7 - Visualize model predictions that simulating the MEG sensor responses that CANNOT CANCEL for coherent and incoherent neural sources.
		- makeSupplementaryFigure1 - Visualize all subject's stimulus-locked and broadband response topographies
		- makeSupplementaryFigure2 - Visualize all subject's V1-V3 model predictions
		- makeSupplementaryFigure3 - Visualize V1 only model predictions for group, example subject and all individual subjects 
		- makeSupplementaryFigure4 - Visualize model predictions for BEM headmodel vs OS headmodelfor group, example subject 
		- makeSupplementaryFigure5 - Visualize model predictions for different visual areas for group, and example subject 
		- makeSupplementaryFigure6 - Visualize model predictions for different stimulus eccentricitie	
		- makeSupplementaryFigure7 - Visualize model predictions for different neural synchrony levels
		- Supplementary Figure for footnote in introduction:
			Visualize SL data component as difference (stim-blank) in power (amplitude squared), instead of plotting amplitudes (not squared).					
	
- external			: code from other toolboxes 
	- mne_alignment : Python code to do additional check MRI/MEG data alignment
	- CircStat2012a : Function to make von Mises distribution by Philip Berens (https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox-directional-statistics)
	- knk_utils 	: Functions from Kendrick Kay's utils toolbox (https://github.com/kendrickkay/knkutils)
	- nppDenoise 	: Functions from NoisePool-PCA denoising toolbox (https://github.com/elinekupers/noisepoolPCADenoise)

- fmsRootPath.m 		: function to call local path of this repository folder

## Data
- headmodel from brainstorm:
	currently on Winawerlab server (_/Projects/MEG/brainstorm_db/SSMEG/anat/wlsubjXXX_)
- templates from Freesurfer:
	currently on Winawerlab server (_/Freesurfer_subjects/wlsubjXXX/surf/ r/lh.benson14_XXX_)
- templates matched to Brainstorm surfaces:
	currently on Winawerlab server (_/Projects/MEG/brainstorm_db/SSMEG/anat/wlsubjXXX/ r/lh.benson14_XXX_)
- MEG data: 
	currently on Winawerlab server (_/Projects/MEG/SSMEG/XX_SSMEG_XXX_wlsubjXXX_ and
								_Projects/MEG/SSMEG/fullOnly/XX_SSMEG_XXX_wlsubjXXX_)

## Workflow and examples
1. How to get downsampled V1-V3 anatomical templates that can be used to make predictions with Brainstorm gain matrix?
	1. Compute FreeSurfer's recon-all segmentation of subject's T1w image (https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all)
	2. Use Benson Docker on Freesurfer subject to create V1-V3 anatomical template(https://github.com/noahbenson/neuropythy/)
	3. Combine and downsample FS hemispheres to BS with this repository script 's_visualAreasFS2BS'

2. How to get V1 anatomical template prediction for a particular subject and visualize it?
	1. Have access to the SSMEG folder in Brainstorm database
	2. Run figure script: `makeFigure1('exampleSubject', subjectnr)`

3. How to get SSMEG data for particular subject and visualize?
	1. Copy following matlab files from _/Projects/MEG/SSMEG_ on server to _data/wlsubjXXX_ folder in repository:
		- _sXX_conditions.m_
		- _sXX_denoisedData_bb.mat_
		- _sXX_denoiseData_sl.mat_
		- _sXX_denoisedts.mat_
		- _sXX_sensorData.mat_
	2. Run figure script: `makeFigure2('exampleSubject', subjectnr)`

4. How to plot data against model prediction for all subjects?
	1. Run figure script: `makeFigure4()`

5. How to plot model prediction with alternative forward model (no cancellation allowed)
	1. Run figure script: `makeFigure6()`

## Dependencies
- Brainstorm toolbox   (to align MEG and MRI and compute gain matrix, https://neuroimage.usc.edu/brainstorm/)
- FieldTrip toolboxes  (for visualizing meshes, http://www.fieldtriptoolbox.org/)
- meg_utils repository (for analysis and visualizing, https://github.com/WinawerLab/meg_utils)

This code was created with MATLAB 9.1 and the following build in toolboxes:
- Bioinformatics Toolbox 	 				4.7 
- Computer Vision System Toolbox 			7.2
- Control System Toolbox 	 				10.1 
- Curve Fitting Toolbox 	 				3.5.4 
- DSP System Toolbox 	 					9.3 
- Database Toolbox 	 					 	7.0 
- Global Optimization Toolbox 			 	3.4.1 
- Image Acquisition Toolbox 	 			5.1 
- Image Processing Toolbox 	 			 	9.5 
- Instrument Control Toolbox 	 			3.10 
- Mapping Toolbox 	 					 	4.4 
- Neural Network Toolbox 	 				9.1 
- Optimization Toolbox 	 				 	.5 
- Parallel Computing Toolbox 	 			6.9 
- Partial Differential Equation Toolbox 	2.3 
- RICOH MEG Reader toolbox for MATLAB 	 	1.0.3 
- Signal Processing Toolbox 	 			7.3 
- Simulink 	 							 	8.8 
- Simulink Control Design 	 			  	4.4 
- Statistics and Machine Learning Toolbox  	11.0 
- Symbolic Math Toolbox 	 				7.1 
- System Identification Toolbox 			9.5 
- Wavelet Toolbox 	 					 	4.17 
- Yokogawa MEG Reader toolbox for MATLAB 	1.5.2 


