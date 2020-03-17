function interp_retinotopy(bsdir, fsdir, sub_fs, sub_bs, proj_bs, highResFlag)

% Function to interpret retinotopy from this subject's freesurfer directory
% as outputted by the hub.docker.com/r/nben/occipital_atlas Docker (in
% volume files) to this subject's brainstorm data.
%
% interp_retinotopy(bsdir, fsdir, sub_fs, sub_bs, proj_bs)
%
% INPUTS: 
%   bsdir:      Path to brainstorm database 
%   fsdir:      Path to freesurfer subjects
%   sub_fs:     name of subject in freesurfer folder
%   sub_bs:     name of subject in brainstorm folder
%   proj_bs:    name of project in brainstorm folder
%
% example:
%   interp_retinotopy([], [], 'wlsubj010', 'wlsubj010', 'SSMEG')
%   interp_retinotopy([], [], 'wlsubj011', 'wlsubj011', 'SSMEG')
%   interp_retinotopy([], [], 'wlsubj011', 'wlsubj011', 'SSMEG')
%   interp_retinotopy([], [], 'wl_subj014', 'wlsubj014', 'GAMMA')
%
% Notes:
%  - the occipital_atlas Docker usually names its files *.template_*.mgz;
%  - You must have already started the brainstorm GUI for this script to work
%    correctly.
%  - You need to have the following toolboxes on your path: 
%       addpath(genpath('/Applications/freesurfer/matlab'));
%       addpath(genpath('/Applications/freesurfer/fsfast/toolbox'));
%       addpath(genpath('~/matlab/git/toolboxes/brainstorm3'));

%% Setup Paths and Configuration

if ~exist('bsdir','var')||isempty(bsdir)
    bsdir = '/Volumes/server/Projects/MEG/brainstorm_db';
end

if ~exist('fsdir','var')||isempty(fsdir)
    fsdir = '/Volumes/server/Freesurfer_subjects';
end

if ~exist('sub_fs','var')||isempty(sub_fs)
    error('Need to specify freesurfer subject to run interp_retinotopy function')
end

if ~exist('sub_bs','var')||isempty(sub_bs)
    error('Need to specify brainstorm subject to run interp_retinotopy function')
end

if ~exist('proj_bs','var')||isempty(proj_bs)
    error('Need to specify brainstorm project to run interp_retinotopy function')
end

if ~exist('highResFlag','var')||isempty(highResFlag)
    highResFlag = 0;
end

sub_fsdir = sprintf('%s/%s', fsdir, sub_fs);
sub_bsdir = sprintf('%s/%s/anat/%s', bsdir, proj_bs, sub_bs);
setenv('SUBJECTS_DIR', fsdir);


%% Load Retinotopy Data

% d = dir(fullfile(sprintf('%s/surf/*benson14_varea*', sub_fsdir)));
% 
% if regexp(d(1).name, 'mgz', 'ONCE')
%     retino_fname = @(hem, type)(sprintf('%s/surf/%s.benson14_%s.mgz', sub_fsdir, hem, type));
%     surfdat      = @(lh, rh, rc)(setfield(setfield([], 'lh', lh.vol(:)), 'rh', rc*rh.vol(:)));
% 
%     sub_fs_angle = surfdat(MRIread(retino_fname('lh', 'angle')), MRIread(retino_fname('rh', 'angle')), -1);
%     sub_fs_eccen = surfdat(MRIread(retino_fname('lh', 'eccen')), MRIread(retino_fname('rh', 'eccen')),  1);
%     sub_fs_areas = surfdat(MRIread(retino_fname('lh', 'varea')), MRIread(retino_fname('rh', 'varea')),  1);
% 
%     
% else
%     retino_fname = @(hem, type)(sprintf('%s/surf/%s.benson14_%s', sub_fsdir, hem, type));
%     surfdat      = @(lh, rh, rc)(setfield(setfield([], 'lh', lh(:)), 'rh', rc*rh(:)));
% 
%     sub_fs_angle = surfdat(read_curv(retino_fname('lh', 'angle')), read_curv(retino_fname('rh', 'angle')), -1);
%     sub_fs_eccen = surfdat(read_curv(retino_fname('lh', 'eccen')), read_curv(retino_fname('rh', 'eccen')),  1);
%     sub_fs_areas = surfdat(read_curv(retino_fname('lh', 'varea')), read_curv(retino_fname('rh', 'varea')),  1);
% end


% Get Wang Atlas ROIs
wang_fname   = @(hem, type)(sprintf('%s/surf/%s.wang2015_atlas.mgz', sub_fsdir, hem));
surfdat      = @(lh, rh, rc)(setfield(setfield([], 'lh', lh.vol(:)), 'rh', rc*rh.vol(:)));
sub_fs_wang  = surfdat(MRIread(wang_fname('lh', 'wang2015_atlas')), MRIread(wang_fname('rh', 'wang2015_atlas')),  1);




%% Apply the interpolation


if highResFlag
    subFolder = 'highres';
    sub_fs_interp = @(dat)(tess_fs2bst(sub_bs, sub_fs, dat.lh, dat.rh, 'tess_cortex_pial_high'));
else
    subFolder = 'lowres';
    sub_fs_interp = @(dat)(tess_fs2bst(sub_bs, sub_fs, dat.lh, dat.rh, 'tess_cortex_pial_low'));
   
end

sub_bsdir = fullfile(sub_bsdir, subFolder);
if ~exist(sub_bsdir, 'dir'); mkdir(sub_bsdir); end
    

% sub_bs_angle = sub_fs_interp(sub_fs_angle);
% sub_bs_eccen = sub_fs_interp(sub_fs_eccen);
% sub_bs_areas = sub_fs_interp(sub_fs_areas);
sub_bs_wang = sub_fs_interp(sub_fs_wang);

sub_bs_angle_mgzfile = sprintf('%s/benson14angle_overlay.mgz', sub_bsdir);
sub_bs_eccen_mgzfile = sprintf('%s/benson14eccen_overlay.mgz', sub_bsdir);
sub_bs_areas_mgzfile = sprintf('%s/benson14areas_overlay.mgz', sub_bsdir);
sub_bs_wang_mgzfile = sprintf('%s/wang15areas_overlay.mgz', sub_bsdir);


sub_bs_angle_matfile = sprintf('%s/benson14angle_overlay.mat', sub_bsdir);
sub_bs_eccen_matfile = sprintf('%s/benson14eccen_overlay.mat', sub_bsdir);
sub_bs_areas_matfile = sprintf('%s/benson14areas_overlay.mat', sub_bsdir);
sub_bs_wang_matfile = sprintf('%s/wang15areas_overlay.mat', sub_bsdir);


% MRIwrite(struct('vol', sub_bs_angle), sub_bs_angle_mgzfile);
% MRIwrite(struct('vol', sub_bs_eccen), sub_bs_eccen_mgzfile);
% MRIwrite(struct('vol', sub_bs_areas), sub_bs_areas_mgzfile);
MRIwrite(struct('vol', sub_bs_wang), sub_bs_wang_mgzfile);

% save(sub_bs_angle_matfile, 'sub_bs_angle');
% save(sub_bs_eccen_matfile, 'sub_bs_eccen');
% save(sub_bs_areas_matfile, 'sub_bs_areas');

save(sub_bs_wang_matfile, 'sub_bs_wang');





return
