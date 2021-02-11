% s_reviewFigures
%
% Script to make all data figures for reviewers of manuscript: 
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
%
% By Eline Kupers, NYU (2020)


%% Supplementary Figure V2/V3 
% Visualize model predictions for different visual areas for group, and example subject 
subjectsToPlot = 1:12;
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true, 'area', 'V23')


%% Supplementary Figure Benson18 Wang15 vs headmodel BEM/OS for review
% Visualize model predictions for different visual areas for group, and example subject 
subjectsToPlot = 1:12;
headmodelType = 'BEM';

makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true, 'area', 'benson18atlas', 'headmodelType', headmodelType)
makeFigure5('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true, 'area', 'wang15atlas', 'headmodelType', headmodelType)


%% New Figure 3 (SL and BB POWER)
% Visualize visually evoked steady state and broadband MEG sensor responses
subjectsToPlot = 12;
makeFigure3('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true, 'useSLPower', true)


%% New Figure 3 (SL coherent spectrum and BB POWER)
% Visualize visually evoked steady state and broadband MEG sensor responses
subjectsToPlot = 1:12;
makeFigure3('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true, 'useSLPower', true, 'useSLIncohSpectrum', false)
