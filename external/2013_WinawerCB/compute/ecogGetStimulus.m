function stimulus = ecogGetStimulus(stimTypes)
% Stimulus is a cell array of stimulus apertures. Each stimulus
% description is pixels x time, where pixels is a vector corresponding to
% the concatenated colums of the binarized image aperture.

% define stimulus
stimulus = cell(1, length(stimTypes));

for ii = 1:length(stimTypes)
    aperture =  getAperture(stimTypes{ii});
    aperture = reshape(aperture, size(aperture, 1)*size(aperture, 2), [])';
    stimulus{ii} = double(aperture);
end


function aperture = getAperture(stim)
% Return the 3-d array of stimulus apertures x time for a given experiment,
% specified by the input stim.

switch stim
    case {'bar' 'barLC' 'barMC'}
        fname = 'bar';
    case {'barThin' 'barThinLC' 'barThinMC'}
        fname = 'barThin';
    case {'barThick' 'barThickLC' 'barThickMC'}
        fname = 'barThick';
    case {'onoff' 'onoffLC'}
        fname = 'OnOff';
    case {'barRandom' 'barRandomThin' 'barRandomBars' 'barRandomWide'}
        fname = 'barRandom';
    case 'barRandomCRF'
        fname = 'contrasts';
    otherwise
        error('Stimulus ''%s'' not recognized', stim)
end

tmp = load(fullfile(ecogPRFrootPath, 'stims', fname));
aperture = tmp.images;

% For experiments with event-related sequence (subject 4 only), stimuli of
% different types were interleaved in each run (for example, thin, normal,
% and thick bars) in random order). the name 'barRandomThin' refers to only
% 1/3 of the stimuli, which, after sorting, corresponds to trials 1:48.
% barRandomBars (medium size bars) corresponds to trials 49:96 (after time
% series has been sorted), and so on.
switch stim
    case 'barRandomThin'
        aperture = aperture(:,:,1:48);
    case 'barRandomBars'
        aperture = aperture(:,:,49:96);
    case 'barRandomWide'
        aperture = aperture(:,:,97:144);
    otherwise
        % do nothing
end

