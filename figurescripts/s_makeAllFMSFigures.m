%% s_makeAllFMSFigures
%
% This is a script to make figures from the manuscript: 
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.

%% FIGURE 2
% Visualize visually evoked steady state and broadband MEG sensor responses
subjectsToPlot = 12;
makeFigure1(subjectsToPlot)

%% FIGURE 3
% Visualize model predictions simulating the MEG sensor responses for 
% coherent and incoherent neural sources.
makeFigure2('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)

%% FIGURE 5
% Visualize model predictions simulating the MEG sensor responses for 
% coherent and incoherent neural sources.
makeFigure4('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)

%% FIGURE 7
% Visualize model predictions that simulating the MEG sensor responses 
% that CANNOT CANCEL for coherent and incoherent neural sources.
makeFigure6('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)


%% Supplementary Figure 1
% Visualize all subject's stimulus-locked and broadband response topographies
subjectsToPlot = 1:12;
makeFigure2('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', false, 'saveFig', false)

%% Supplementary Figure 2 
% Visualize all subject's V1-V3 model predictions
subjectsToPlot = 1:12;
makeFigure4('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', false, 'saveFig', false)

%% Supplementary Figure 3 
% Visualize V1 only model predictions for group, example subject and all individual subjects 
subjectsToPlot = 1:12;
area = 'V1';
makeFigure4('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true, 'area', area)

%% Supplementary Figure 4
% Visualize model predictions for BEM headmodel vs OS headmodelfor group, example subject 
subjectsToPlot = 12;
headmodelType  = 'BEM';
makeFigure4('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true, 'headmodelType', headmodelType)

%% Supplementary Figure 5
% Visualize model predictions for BEM headmodel vs OS headmodelfor group, example subject 
subjectsToPlot = 12;
areas  = {'benson18atlas', 'wang15atlas'};
for a = 1:length(areas)
    makeFigure4('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true, 'area', areas{a})
end

%% Supplementary Figure 6
% Visualize model predictions for different stimulus eccentricities
subjectsToPlot = 12;
makeSupplementaryFigure6('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true)

%% Supplementary Figure 7
% Visualize model predictions for different neural synchrony levels
subjectsToPlot = 12;
makeSupplementaryFigure7('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true)

%% Supplementary Figure for footnote in introduction
% Visualize SL data component as difference (stim-blank) in power (amplitude squared),
% instead of plotting amplitudes (not squared).
subjectsToPlot = 12;
makeFigure2('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', true, 'useSLPower', true)



