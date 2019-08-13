%% s_makeAllFMSFigures


% This is a script to make figures from the manuscript: XXX by XXX

%% FIGURE 3
% Visualize visually evoked steady state and broadband MEG sensor responses
subjectsToPlot = 12;
makeFigure3(subjectsToPlot);

%% FIGURE 4
% Visualize model predictions simulating the MEG sensor responses for 
% coherent and incoherent neural sources.
makeFigure4('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)

%% FIGURE 4
% Visualize model predictions simulating the MEG sensor responses for 
% coherent and incoherent neural sources.
makeFigure4('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)

%% FIGURE 6
% Visualize model predictions that simulating the MEG sensor responses 
% that CANNOT CANCEL for coherent and incoherent neural sources.
makeFigure4('subjectsToPlot', subjectsToPlot, 'plotMeanSubject', true, 'saveFig', false)