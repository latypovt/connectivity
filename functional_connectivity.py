import argparse
import numpy as np
import nibabel as nib
import pandas as pd
import os
import pandas as pd
from nilearn.input_data import NiftiLabelsMasker
from nilearn.connectome import ConnectivityMeasure
from nilearn.interfaces.fmriprep import load_confounds, load_confounds_strategy
from tqdm import tqdm
from nilearn import datasets
import shutil
import json

def get_labels(parc='DK'):
    if parc == 'DK':
        labels = pd.read_csv('/d/gmi/1/timurlatypov/apps/connectivity/fs_default.csv', index_col=0, header=None)
    elif parc == 'Destrieux':
        labels = pd.read_csv('/d/gmi/1/timurlatypov/apps/connectivity/fs_a2009s.csv', index_col=0, header=None)
    labels = labels[1].values.tolist()
    return labels

def extract_time_series(fmri_file, parc_file, confounds_file=None, list_of_confounds=None):
    #if confounds_file and list_of_confounds:
    confounds_list = [fmri_file, confounds_file, list_of_confounds]
    fmri_img = nib.load(fmri_file)
    parc_img =  nib.load(parc_file)
    confounds, sample_mask = load_confounds(fmri_file, strategy=['high_pass', 'motion', 'wm_csf', 'global_signal'], motion='full', wm_csf='full', global_signal='full')
    labels = get_labels(parc='DK')

    masker = NiftiLabelsMasker(labels_img=parc_img, labels=labels, smoothing_fwhm=6, keep_masked_labels=True, t_r=1.5)
    time_series = masker.fit_transform(fmri_img, confounds=confounds, sample_mask=sample_mask)
    returned_labels = masker.region_names_
    return time_series, returned_labels

def compute_functional_connectivity(time_series):
    correlation_measure = ConnectivityMeasure(kind='correlation')
    correlation_matrix = correlation_measure.fit_transform([time_series])[0]
    return correlation_matrix

def list_sessions(subject_dir):
    return [d for d in os.listdir(subject_dir) if os.path.isdir(os.path.join(subject_dir, d)) and d.startswith('ses-')]

def list_runs(session_dir):
    return [f for f in os.listdir(session_dir) if f.endswith('preproc_bold.nii.gz')]

def main():
    parser = argparse.ArgumentParser(description="Compute functional connectivity matrix from fmriprep output.")
    parser.add_argument('subject_id', type=str, help='Subject ID')
    parser.add_argument('output_dir', type=str, help='Directory to save the connectivity matrices')
    args = parser.parse_args()

    subject_dir =  args.subject_id 
    sessions = list_sessions(subject_dir)

    # Load the aparc+aseg atlas


    for ses_id, session in enumerate(sessions):
        session_dir = os.path.join(subject_dir, session, 'func')
        runs = list_runs(session_dir)
        for index, run in enumerate(tqdm(runs, desc=f'Functional {subject_dir, session}')):
            fmri_file = os.path.join(session_dir, run)
            fmri_mask_file = fmri_file.replace('_space-T1w_desc-preproc_bold.nii.gz', '_desc-brain_mask.nii.gz')
            parc_file = fmri_file.replace('task-rest_space-T1w_desc-preproc_bold.nii.gz', 'DK_T1space.nii.gz') 
            parc_file = parc_file.replace('func', 'parc')
            confound_file = fmri_file.replace('space-T1w_desc-preproc_bold.nii.gz', 'desc-confounds_timeseries.tsv')
            confound_json = confound_file.replace('.tsv', '.json')
            time_series, out_labels = extract_time_series(fmri_file, parc_file, confounds_file=confound_file, list_of_confounds=confound_json)
            correlation_matrix = compute_functional_connectivity(time_series)

            output_file = os.path.join(args.output_dir, f'{args.subject_id}_ses-{ses_id+1}_run-{index+1}_fc_matrix.npy')
            np.save(output_file, correlation_matrix)
            json.dump(out_labels, open(output_file.replace('.npy', '.json'), 'w'))
if __name__ == "__main__":
    main()