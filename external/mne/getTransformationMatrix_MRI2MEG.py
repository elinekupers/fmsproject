import mne
from mne import io
from mne.datasets import sample


import numpy as np
from mayavi import mlab
import os

dir_name = '/Volumes/server/Projects/MEG/SSMEG/fullOnly'

## Session 1
# subj_pth1 = '01_SSMEG_R1372_03.11.2018_wlsubj048/raw/'
# input_fname = os.path.join(dir_name, subj_pth1, 'R1372_SSMEG_Block1_3.12.18.sqd')
# mrk = os.path.join(dir_name, subj_pth1, 'R1372_Marker1_3.12.18.sqd')
# elp = os.path.join(dir_name, subj_pth1, 'R1372_3.12.18.elp')
# hsp = os.path.join(dir_name, subj_pth1, 'R1372_3.12.18.hsp')

## Session 2
# subj_pth2 = '02_SSMEG_R1373_03.13.2018_wlsubj046/raw/'
# input_fname = os.path.join(dir_name, subj_pth2, 'R1373_SSMEG_3.13.18.sqd')
# mrk = os.path.join(dir_name, subj_pth2, 'R1373_Marker1_3.13.18.sqd')
# elp = os.path.join(dir_name, subj_pth2, 'R1373_3.13.18.elp')
# hsp = os.path.join(dir_name, subj_pth2, 'R1373_3.13.18.hsp')
# output = os.path.join(dir_name, subj_pth2, 'wlsubj046_SSMEG.fif')

## Session 3
# subj_pth3 = '03_SSMEG_R1374_03.14.2018_wlsubj039/raw/'
# input_fname = os.path.join(dir_name, subj_pth3, 'R1374_SSMEG_Block1_03.14.2018.sqd')
# mrk = os.path.join(dir_name, subj_pth3, 'R1374_Marker1_03.14.2018.sqd')
# elp = os.path.join(dir_name, subj_pth3, 'R1374_03.14.18.elp')
# hsp = os.path.join(dir_name, subj_pth3, 'R1374_03.14.18.hsp')
# output = os.path.join(dir_name, subj_pth3, 'wlsubj039_SSMEG.fif')

## Session 4
# subj_pth4 = '04_SSMEG_R1021_08.30.2018_wlsubj059/raw/'
# input_fname = os.path.join(dir_name, subj_pth4, 'R1021_SSMEG2_run1-8.sqd')
# mrk = os.path.join(dir_name, subj_pth4, 'R1021_Marker1_8.30.18.sqd')
# elp = os.path.join(dir_name, subj_pth4, 'R1021_8.30.18.elp')
# hsp = os.path.join(dir_name, subj_pth4, 'R1021_8.30.18.hsp')
# output = os.path.join(dir_name, subj_pth4, 'wlsubj059_SSMEG.fif')

## Session 5
subj_pth5 = '05_SSMEG_R1109_11.12.2018_wlsubj067/raw/'
input_fname = os.path.join(dir_name, subj_pth5, 'R1109_SSMEG_fullOnly_11.12.18.sqd')
mrk = os.path.join(dir_name, subj_pth5, 'R1109_Marker1_11.12.18.sqd')
elp = os.path.join(dir_name, subj_pth5, 'R1109_11.12.18.elp')
hsp = os.path.join(dir_name, subj_pth5, 'R1109_11.12.18.hsp')
output = os.path.join(dir_name, subj_pth5, 'wlsubj067_SSMEG.fif')

stim = list(range(161,168))

raw = mne.io.read_raw_kit(input_fname, mrk=mrk, elp=elp, hsp=hsp, stim=stim, allow_unknown_format=True)
raw.save(output)
raw.close()
    

# Create fif file via bash (not sure how it works without additional allow_unknown_format=true)
# mne kit2fiff --input='/Volumes/server/Projects/MEG/SSMEG/fullOnly/02_SSMEG_R1373_03.13.2018_wlsubj046/raw/R1373_SSMEG_3.13.18.sqd' --mrk='/Volumes/server/Projects/MEG/SSMEG/fullOnly/02_SSMEG_R1373_03.13.2018_wlsubj046/raw/R1373_Marker1_3.13.18.sqd' --elp='/Volumes/server/Projects/MEG/SSMEG/fullOnly/02_SSMEG_R1373_03.13.2018_wlsubj046/raw/R1373_3.13.18.elp' --hsp='/Volumes/server/Projects/MEG/SSMEG/fullOnly/02_SSMEG_R1373_03.13.2018_wlsubj046/raw/R1373_3.13.18.hsp' --output='/Volumes/server/Projects/MEG/SSMEG/fullOnly/02_SSMEG_R1373_03.13.2018_wlsubj046/raw/wlsubj046_SSMEG.fif' --debug


# Compute Boundary Element Model (BEM) to create 3 head surfaces
mne make_scalp_surfaces --subject='wlsubj067'  --force --overwrite --verbose --subjects-dir='/Volumes/server/Freesurfer_subjects'

# Open GUI with mn coreg
mne coreg --subjects-dir='/Volumes/server/Freesurfer_subjects' --subject='wlsubj048' --fiff='/Volumes/server/Projects/MEG/SSMEG/fullOnly/01_SSMEG_R1372_03.11.2018/raw/R1372_SSMEG.fif'

# Read the trans.fif file
trans = mne.read_trans('/Volumes/server/Projects/MEG/SSMEG/fullOnly/01_SSMEG_R1372_03.11.2018_wlsubj048/raw/wlsubj048-trans.fif')
trans = mne.read_trans('/Volumes/server/Projects/MEG/SSMEG/fullOnly/02_SSMEG_R1373_03.13.2018_wlsubj046/raw/wlsubj046-trans.fif')
trans = mne.read_trans('/Volumes/server/Projects/MEG/SSMEG/fullOnly/03_SSMEG_R1374_03.14.2018_wlsubj039/raw/wlsubj039-trans.fif')
trans = mne.read_trans('/Volumes/server/Projects/MEG/SSMEG/fullOnly/04_SSMEG_R1021_08.30.2018_wlsubj059/raw/wlsubj059-trans.fif')
trans = mne.read_trans('/Volumes/server/Projects/MEG/SSMEG/fullOnly/05_SSMEG_R1109_11.12.2018_wlsubj067/raw/wlsubj067-trans.fif')

# [[ 0.9993735   0.03445512 -0.00808774  0.00507086]
#  [-0.03386335  0.86448795 -0.50151157  0.02088116]
#  [-0.01028789  0.50147128  0.86511314 -0.09503751]
#  [ 0.          0.          0.          1.        ]]



# Visualize alignment
subjects_dir = os.path.join(dir_name, subj_pth)
raw_fname = os.path.join(subjects_dir, 'R1372_SSMEG.fif')
trans_fname = os.path.join(subjects_dir, 'wlsubj048-trans.fif')

raw = mne.io.read_raw_fif(raw_fname)
trans = mne.read_trans(trans_fname)

mne.viz.plot_alignment(raw.info, trans=trans, subject='wlsubj048',
                       subjects_dir='/Volumes/server/Freesurfer_subjects', surfaces='head-dense',
                       show_axes=True, dig=True, eeg=[], meg='sensors',
                       coord_frame='meg')

