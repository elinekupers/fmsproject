%% Script to plot Figure 4ab
%
% Winawer, Kay, Foster, Parvizi, Wandell.
% *Asynchronous broadband signals are the principal source of the BOLD
% response in human visual cortex*
% _Current Biology, 2013_
%
% This figure compares the prf model accuracy for a linear fit and a CSS
% fit, for both the asynchronous broadband (panel a) and stimulus locked
% stimulus-locked (pabel b) time series.
%
% Copyright Jonathan Winawer, 2013

%%
% A note on reproducibility: The plots produced from this script differ
% very slightly from the plots in Figure 4 in the publication. This is due
% to a slight difference in the stimulus descriptions used for solving the
% pRF models. In this script, the stimuli used as inputs for the pRF models
% are binary masks. For the publication the stimuli were floats
% approximating a binary mask, with occasional pixel values differing
% slightly from 0 or 1 due to an imperfection in the algorithm that
% converted the image indices used for the experiments into into binary
% contrast masks. 

%% Set up paths

% Path to save the eps figures
savepth   = fullfile(ecogPRFrootPath, 'scratch');


%% Load the PRF model solutions
%
% The pRF models are pre-solved. If you would like to resolve them, run the
% following script:
%%
%  s_ecogSolvePRFs

% Load the precomputed AB and SL pRF model solutions. This includes the CSS
% model ('exp') and the linear model ('noexp') for each component of the
% ECoG signal.
prf.abExp   = load(fullfile(ecogPRFrootPath, 'data', 'PRF_xval_exp_ab'));
prf.abNoExp = load(fullfile(ecogPRFrootPath, 'data', 'PRF_xval_noexp_ab'));
prf.slExp   = load(fullfile(ecogPRFrootPath, 'data', 'PRF_xval_exp_sl'));
prf.slNoExp = load(fullfile(ecogPRFrootPath, 'data', 'PRF_xval_noexp_sl'));

%% Select appropriate channels
% For summarizing population data, we impose the following selection
% criterion:
%%
% # Variance explained from either the exp or noexp model exceeds 20%
% # This is determined separately for the AB and SL data sets
%
% For further details, see supplementary table 1 and supplementary methods
% section 'Channel Selection'.

abok = prf.abExp.params.r > 20 | prf.abNoExp.params.r > 20;
slok = prf.slExp.params.r > 20 | prf.slNoExp.params.r > 20;


%% Figure 4a: Scatterplot of AB variance explained, CSS versus Linear models
fH = figure;set(fH, 'Color', 'w');
set(gca, 'FontSize', 16, 'XLim' ,[0 100], 'YLim', [0 100]); axis square
hold on;

title('Cross-validated variance explained, Asynchronous Broadband')
xlabel('Linear model')
ylabel('CSS model')

colors = 'rgbw';

for subj = 1:4
    
    % Plot V1V2V3 channels as filled circles. Rectify (ie, negative r becomes
    % 0, otherwise we lose points on the plots)
    toPlot = abok & prf.abExp.params.subj == subj & prf.abExp.params.isV1V2V3;
    x = max(prf.abNoExp.params.r(toPlot), 0);
    y = max(prf.abExp.params.r(toPlot), 0);
    plot(x,y, 'ko', 'MarkerFaceColor', colors(subj), 'MarkerSize', 12);

    
    % Plot channels outside V!V2V3 as filled diamonds
    toPlot = abok & prf.abExp.params.subj == subj & ~prf.abExp.params.isV1V2V3;
    x = max(prf.abNoExp.params.r(toPlot), 0);
    y = max(prf.abExp.params.r(toPlot), 0);
    plot(x,y, 'kd', 'MarkerFaceColor', colors(subj), 'MarkerSize', 12);

end

% Plot the identity line
plot([0 100], [0 100], 'k-')

%% Figure 4b: Scatterplot of SL variance explained, CSS versus Linear models
fH(2) = figure;set(fH, 'Color', 'w');
set(gca, 'FontSize', 16, 'XLim' ,[0 100], 'YLim', [0 100]); axis square
hold on;

title('Cross-validated variance explained, Stimulus-locked')
xlabel('Linear model')
ylabel('CSS model')

colors = 'rgbw';

for subj = 1:4
    
    % Plot V1V2V3 channels as filled circles. Rectify (ie, negative r becomes
    % 0, otherwise we lose points on the plots)
    toPlot = slok & prf.slExp.params.subj == subj & prf.slExp.params.isV1V2V3;
    x = max(prf.slNoExp.params.r(toPlot), 0);
    y = max(prf.slExp.params.r(toPlot), 0);
    plot(x,y, 'ko', 'MarkerFaceColor', colors(subj), 'MarkerSize', 12);

    
    % Plot channels outside V!V2V3 as filled diamonds
    toPlot = slok & prf.slExp.params.subj == subj & ~prf.slExp.params.isV1V2V3;
    x = max(prf.slNoExp.params.r(toPlot), 0);
    y = max(prf.slExp.params.r(toPlot), 0);
    plot(x,y, 'kd', 'MarkerFaceColor', colors(subj), 'MarkerSize', 12);

end

% Plot the identity line
plot([0 100], [0 100], 'k-')


%% Save 'em

hgexport(fH(1), fullfile(savepth, 'Figure4a_AB_pRF_Accuracy.eps'));
hgexport(fH(2), fullfile(savepth, 'Figure4b_SL_pRF_Accuracy.eps'));

return