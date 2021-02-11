function rootPath = ecogPRFrootPath
% Return the path to the root ECoG PRF directory
%
% This function must reside in the directory at the base of the ECoG PRF
% directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(ecogPRFrootPath,'data')

rootPath=which('ecogPRFrootPath');

rootPath = fileparts(rootPath);

return