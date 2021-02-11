function [x, y, s] = ecogPix2Deg(subj, x0, y0, s0)
%Convert pixels to degrees for pRF solutions
%
% [x, y, s] = ecogFitPRFPix2Deg(subj, x0, y0, s0)
%
% We compute PRF models in pixels, but would like to report them in
% degrees.
%
% Inputs:   
%   subj: subject number (1-4)
%   x0:   x center of pRF in pixels
%   y0:   y center of pRF in pixels
%   s0:   pRF size in pixels (pRF size = sigma/sqrt(n))
%
% Outputs
%   x:  x center of pRF in degrees of visual angle (0 is fixation)
%   y:  y center of pRF in degrees of visual angle (0 is fixation)
%   s:  pRF size in degrees

% S1: 61 cm
% S2: 54 cm
% S3: 57.5 cm
% S4: 47 cm

% Make sure all input arguments are row vectors

nPoints = max([length(x0) length(y0) length(s0)]);
if isempty(x0), x0 = zeros(1,nPoints); end
if isempty(y0), y0 = zeros(1,nPoints); end
if isempty(s0), s0 = zeros(1,nPoints); end

if length(subj) == 1, subj = repmat(subj, nPoints, 1); end
x0 = x0(:)'; y0 = y0(:)'; s0 = s0(:)'; subj = subj(:)';

% get viewing distance in cm (different for each subject)
vd = zeros(size(subj));
for ii = 1:length(subj)
    switch subj(ii)
        case 1, vd(ii) = 61; %cm
        case 2, vd(ii) = 54;
        case 3, vd(ii) = 57.5;
        case 4, vd(ii) = 47;
        otherwise
            error('Unknown viewing distance for subjec %d', subj);
    end
end

% height from center of screen to top of 15 inch apple macbook (gunjou or
% blue moon) used for ecog experiments
rcm = 20.7/2; % radius in cm
rdeg = rad2deg(atan(rcm./vd)); % radius in degrees

numpix  = 101; % number of pixels in stimulus description (one side) as input to fitprf
pix2deg = rdeg*2/numpix;
xC      = (numpix+1)/2; % center location in pixels
yC      = xC;            

% convert pixels to degrees
x = pix2deg.*(x0-xC);
y = pix2deg.*(y0-yC);
s = pix2deg.*s0;

end