import mne
from mne import io
from mne.datasets import sample


import numpy as np
from mayavi import mlab
import os

dir_name = '/Volumes/server/Projects/MEG/SSMEG/fullOnly'
subj_pth = '01_SSMEG_R1372_03.11.2018_wlsubj048/raw/'
input_fname = os.path.join(dir_name, subj_pth, 'R1372_SSMEG_Block1_3.12.18.sqd')
mrk = os.path.join(dir_name, subj_pth, 'R1372_Marker1_3.12.18.sqd')
elp = os.path.join(dir_name, subj_pth, 'R1372_3.12.18.elp')
hsp = os.path.join(dir_name, subj_pth, 'R1372_3.12.18.hsp')
stim = list(range(161,168))

raw = mne.io.read_raw_kit(input_fname, mrk=mrk, elp=elp, hsp=hsp, stim=stim, allow_unknown_format=True)

# Compute Boundary Element Model (BEM) to create 3 head surfaces
mne make_scalp_surfaces --subject='wlsubj048'  --force --overwrite --verbose --subjects-dir='/Volumes/server/Freesurfer_subjects'

mne coreg --subjects-dir='/Volumes/server/Freesurfer_subjects' --subject='wlsubj048' --fiff='/Volumes/server/Projects/MEG/SSMEG/fullOnly/01_SSMEG_R1372_03.11.2018/raw/R1372_SSMEG.fif'


# trans = mne.read_trans('/Volumes/server/Projects/MEG/SSMEG/fullOnly/01_SSMEG_R1372_03.11.2018_wlsubj048/raw/wlsubj048-trans.fif')

# [[ 0.9993735   0.03445512 -0.00808774  0.00507086]
#  [-0.03386335  0.86448795 -0.50151157  0.02088116]
#  [-0.01028789  0.50147128  0.86511314 -0.09503751]
#  [ 0.          0.          0.          1.        ]]




subjects_dir = os.path.join(dir_name, subj_pth)
raw_fname = os.path.join(subjects_dir, 'R1372_SSMEG.fif')
trans_fname = os.path.join(subjects_dir, 'wlsubj048-trans.fif')

raw = mne.io.read_raw_fif(raw_fname)
trans = mne.read_trans(trans_fname)

mne.viz.plot_alignment(raw.info, trans=trans, subject='wlsubj048',
                       subjects_dir='/Volumes/server/Freesurfer_subjects', surfaces='head-dense',
                       show_axes=True, dig=True, eeg=[], meg='sensors',
                       coord_frame='meg')


mne.viz.plot_bem(subject=subject, subjects_dir=subjects_dir,
                 brain_surfaces='white', orientation='coronal')