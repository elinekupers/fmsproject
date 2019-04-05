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

Content:
- external 			: functions from other repositories
- figurescripts		: contains separate functions to make manuscript figures 1-4:
						makeFigure1 - compute and plot the forward model prediction for 1 or average across all subjects
						makeFigure1Mixtures - same as makeFigure1, but to compute many in between model predictions with different values for von Mises kappa (width of the distribution) 
						makeFigure2 - load and plot measured MEG data for 1 or average across all subjects
						makeFigure3 - load and plot measured MEG data and overlay with contourlines of model predictions for all subjects
						makeFigure4 - compute adjusted forward model, that does not allow for responses to cancel					
	- subfunctions 	: functions needed to run makeFigureXX
- mne_alignment		: folder with Python script to check MRI/MEG data alignment

fmsMakeAllFigures.m : master script that calls all makeFigureX functions
fmsRootPath.m 		: function to call local path of this repository folder

## Data
- headmodel from brainstorm:
	currently on Winawerlab server (/Projects/MEG/brainstorm_db/SSMEG/anat/wlsubjXXX)
- templates from Freesurfer:
	currently on Winawerlab server (/Freesurfer_subjects/wlsubjXXX/surf/ r/lh.benson14_XXX)
- templates matched to Brainstorm surfaces:
	currently on Winawerlab server (/Projects/MEG/brainstorm_db/SSMEG/anat/wlsubjXXX/ r/lh.benson14_XXX)
- MEG data: 
	currently on Winawerlab server (/Projects/MEG/SSMEG/XX_SSMEG_XXX_wlsubjXXX and
								Projects/MEG/SSMEG/fullOnly/XX_SSMEG_XXX_wlsubjXXX)

## Workflow and examples
1. How to get downsampled V1-V3 Templates?
	1.1 Compute FreeSurfer's recon-all of subject's T1w image (https://surfer.nmr.mgh.harvard.edu/fswiki/recon-all)
	1.2 Use Benson Docker on Freesurfer subject to create V1-V3 (https://github.com/noahbenson/neuropythy/)
	1.3 Combine and downsample FS hemispheres to BS with this repository script 's_visualAreasFS2BS'

2. How to get V1 template prediction for a particular subject and visualize it?
	2.1 Have access to the SSMEG folder in Brainstorm database
	2.2 Run figure script: makeFigure1('exampleSubject', subjectnr)

3. How to get SSMEG data for particular subject and visualize?
	3.1 Copy following matlab files from /Projects/MEG/SSMEG on server to data/wlsubjXXX folder in repository:
		sXX_conditions.m
		sXX_denoisedData_bb.mat
		sXX_denoiseData_sl.mat
		sXX_denoisedts.mat
		sXX_sensorData.mat
	3.2 Run figure script: makeFigure2('exampleSubject', subjectnr)

4. How to plot data against model prediction for all subjects?
	4.1 Run figure script: makeFigure3()

5. How to plot model prediction with alternative forward model (no cancellation allowed)
	5.1 Run figure script: makeFigure4()


