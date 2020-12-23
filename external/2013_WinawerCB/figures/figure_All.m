%% figure_All
% Script to generate all figures

% Make all figures
figure1c_onOffTS
figure1d_onOffSpectra
figure2cd_AB_and_SL_timeSeries
figure3bc_pRF_Centers_Exponents
figure4ab_PRFcrossvalidation
figure5_Data
figure5_Simulation
figureS5_TravellingWave

% Clean up
close all; clear all

% See what was made
ls(fullfile(ecogPRFrootPath, 'Scratch')) 
