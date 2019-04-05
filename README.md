# fmsproject

Welcome to the code repository of the MEG Forward Models and neural Synchrony (FMS) project

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