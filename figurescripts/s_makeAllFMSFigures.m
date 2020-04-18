%% s_makeAllFMSFigures
%
% Script to make all data figures from the manuscript: 
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
%
% By Eline Kupers, NYU (2019)

%% FIGURE 2
% Visualize visually evoked steady state and broadband MEG sensor responses
subjectsToPlot = 12;
makeFigure2('subjectsToPlot', subjectsToPlot)

%% FIGURE 3
% Visualize model predictions simulating the MEG sensor responses for 
% coherent and incoherent neural sources.
makeFigure3('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)

%% FIGURE 5
% Visualize model predictions simulating the MEG sensor responses for 
% coherent and incoherent neural sources.
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)

%% FIGURE 7
% Visualize model predictions that simulating the MEG sensor responses 
% that CANNOT CANCEL for coherent and incoherent neural sources.
makeFigure7('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)


%% Supplementary Figure 1
% Visualize all subject's stimulus-locked and broadband response topographies
subjectsToPlot = 1:12;
makeFigure3('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', false, 'saveFig', false)

%% Supplementary Figure 2 
% Visualize all subject's V1-V3 model predictions
subjectsToPlot = 1:12;
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', false, 'saveFig', false)

%% Supplementary Figure 3 
% Visualize V1 only model predictions for group, example subject and all individual subjects 
subjectsToPlot = 1:12;
area = 'V1';
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false, 'area', area)

%% Supplementary Figure 4
% Visualize model predictions for BEM headmodel vs OS headmodelfor group, example subject 
subjectsToPlot = 12;
headmodelType  = 'BEM';
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false, 'headmodelType', headmodelType)

%% Supplementary Figure 5
% Visualize model predictions for BEM headmodel vs OS headmodelfor group, example subject 
subjectsToPlot = 12;
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false, 'area', 'V123')
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false, 'area', 'benson18atlas')
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false, 'area', 'wang15atlas')


%% Supplementary Figure 6
% Visualize model predictions for different stimulus eccentricities
subjectsToPlot = 12;
makeSupplementaryFigure6('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)

%% Supplementary Figure 7
% Visualize model predictions for different neural synchrony levels
subjectsToPlot = 12;
makeSupplementaryFigure7('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)

%% Supplementary Figure for footnote in introduction
% Visualize SL data component as difference (stim-blank) in power (amplitude squared),
% instead of plotting amplitudes (not squared).
subjectsToPlot = 12;
makeFigure3('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false, 'useSLPower', true)



