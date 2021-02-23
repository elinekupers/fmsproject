%% s_makeAllFMSFigures
%
% Script to make all data figures from the manuscript: 
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
%
% By Eline Kupers, NYU (2019)

%% General params
saveFig = false;

%% FIGURE 2
% Visualize spectrum of single MEG sensor response
subjectsToPlot = 12;
makeFigure2('subjectsToPlot', subjectsToPlot)

%% FIGURE 3
% Visualize MEG sensor data on topographic maps for stimulus-locked and broadband.
makeFigure3('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig)

%% FIGURE 5
% Visualize model predictions simulating the MEG sensor responses for 
% coherent and incoherent neural sources.
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig)

%% FIGURE 8
% Visualize model predictions that simulating the MEG sensor responses 
% that CANNOT CANCEL for coherent and incoherent neural sources.
makeFigure8('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig)


%% Supplementary Figure 1
% Visualize all subject's stimulus-locked and broadband response topographies
subjectsToPlot = 1:12;
makeFigure3('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', false, 'saveFig', saveFig)


%% Supplementary Figure 2 
% Visualize broadband responses in smaller (10 Hz) frequency bins
subjectsToPlot = 12;
makeSupplementaryFigure2('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig)

%% Supplementary Figure 3 
% Visualize stimulus-locked response topographies computed as geometric
% mean from 12 Hz + all harmonics up to 150 Hz (except 60, 120 Hz, due to line noise).
subjectsToPlot = 12;
makeFigure3('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'amplitudeType','amplitudesHigherHarmonics','saveFig', saveFig)

%% Supplementary Figure 4 
% Visualize all subject's V1-V3 model predictions using simple coherent/incoherent sine wave simulations
subjectsToPlot = 1:12;
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', false, 'saveFig', saveFig)

%% Supplementary Figure 5 
% Visualize V1 only model predictions for group, example subject and all individual subjects 
subjectsToPlot = 12;
area = 'V1';
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig, 'area', area)

area = 'V23';
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig, 'area', area)

%% Supplementary Figure 6
% Visualize model predictions using a more complete, biologically plausible
% simulation of ECOG signals in V1-V3 cortical locations
makeSupplementaryFigure6('subjectsToPlot', 1:12, 'plotMeanSubject', false, 'saveFig', saveFig)

%% Supplementary Figure 7
% Visualize model predictions for different visual areas for group, and
% example subject: ROIs from Benson18 or Wang15 vs headmodels BEM or OS

% Panel A & B: Overlapping spheres head model
subjectsToPlot = 12;
headmodelType = 'OS';

makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig, 'area', 'benson18atlas', 'headmodelType', headmodelType)
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig, 'area', 'wang15atlas', 'headmodelType', headmodelType)

% Panel C & D: Boundary element head model
headmodelType = 'BEM';

makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig, 'area', 'benson18atlas', 'headmodelType', headmodelType)
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig, 'area', 'wang15atlas', 'headmodelType', headmodelType)

%% Supplementary Figure 8
% Visualize model predictions for different stimulus eccentricities
subjectsToPlot = 12;
makeSupplementaryFigure8('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig)

%% Supplementary Figure 9
% Visualize model predictions for different neural synchrony levels
subjectsToPlot = 12;
makeSupplementaryFigure9('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig)

%% Supplementary Figure for footnote in introduction
% Visualize SL data component as difference (stim-blank) in power (amplitude squared),
% instead of plotting amplitudes (not squared).
subjectsToPlot = 12;
makeFigure3('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', saveFig, 'useSLPower', true, 'useBBPower', true)

%% Figure 1, panel 3 showing separate time points for simulated V1-V3 time series
makeFigure1('subjectsToPlot', 1, 'saveFig', saveFig)



