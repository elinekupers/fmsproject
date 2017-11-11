function rootPath = fmsRootPath()
% Return the path to the root forward model synchrony directory
%
% This function must reside in the directory at the base of the
% forward modeling code directory structure.  It is used to determine the location of various
% sub-directories.
% 
% Example:
%   fullfile(DFDrootpath,'code')

rootPath=which('fmsRootPath');

rootPath=fileparts(rootPath);

return
