% s_reviewFigures
%
% Script to make all data figures for reviewers of manuscript: 
%   A visual encoding model links magnetoencephalography
%   signals to neural synchrony in human cortex.
%       by Kupers, Benson, Winawer (YEAR) JOURNAL.
%
%
% By Eline Kupers, NYU (2020)

%% General params
saveFig = false;

%% Make movie of observed and simulated time series
makeFigure3Movie('subjectToPlot', 1)
makeFigure5Movie('subjectToPlot', 1)

%% Plot topomaps of 11, 12, 13 Hz amplitudes using coherent spectrum
makeSupplementaryFigureXX_coherentSpectrum('subjectsToPlot',12, 'plotMeanSubject',true,'singleColorbarFlag', true, 'saveFig', saveFig)

%% Plot coherent spectrum of single and average across subject
makeFigure3('subjectsToPlot', 12, 'plotMeanSubject', true, 'amplitudeType', 'amplitudesCoherentSpectrum', 'saveFig', saveFig)

makeFigure3('subjectsToPlot', 12, 'plotMeanSubject', true, 'withoutS10', true, 'saveFig', saveFig)

