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
    confounds, sample_mask = load_confounds(fmri_file, strategy=['high_pass', 'motion', 'wm_csf', 'global_signal'], motion='basic', wm_csf='basic', global_signal='basic')
    labels = get_labels(parc='DK')

    masker = NiftiLabelsMasker(labels_img=parc_img, smoothing_fwhm=6, keep_masked_labels=True, t_r=1.5)
    time_series = masker.fit_transform(fmri_img, confounds=confounds, sample_mask=sample_mask)
    returned_labels = []#masker.region_names_
    return time_series, returned_labels

def compute_functional_connectivity(time_series):
    correlation_measure = ConnectivityMeasure(kind='correlation')
    correlation_matrix = correlation_measure.fit_transform([time_series])[0]
    return correlation_matrix

def list_sessions(subject_dir):
    list = [d for d in os.listdir(subject_dir) if os.path.isdir(os.path.join(subject_dir, d)) and d.startswith('ses-')]
    return list

def list_runs(session_dir):
    return [f for f in os.listdir(session_dir) if f.endswith('preproc_bold.nii.gz')]

def main():
    parser = argparse.ArgumentParser(description="Compute functional connectivity matrix from fmriprep output.")
    parser.add_argument('subject_id', type=str, help='Subject ID')
    parser.add_argument('output_dir', type=str, help='Directory to save the connectivity matrices')
    parser.add_argument('--no_ses', action='store_true', help='No session directory')
    args = parser.parse_args()

    subject_dir =  args.subject_id 
    if args.no_ses:
        sessions = [subject_dir]
    else:
        sessions = list_sessions(subject_dir)

    for ses_id, session in enumerate(sessions):
        if args.no_ses:
            session_dir = os.path.join(subject_dir, 'func')
        else:
            session_dir = os.path.join(subject_dir, session, 'func')
        runs = list_runs(session_dir)

        session_pbar = tqdm(enumerate(runs), desc=f'Functional {subject_dir, session}', total=len(runs), leave=False)
        for index, run in session_pbar:
            fmri_file = os.path.join(session_dir, run)
            parc_file = [d for d in os.listdir(os.path.join(subject_dir,'parc')) if d.endswith('DK_T1space.nii.gz')]
            parc_file = os.path.join(subject_dir, 'parc', parc_file[0])
            confound_file = fmri_file.replace('space-T1w_desc-preproc_bold.nii.gz', 'desc-confounds_timeseries.tsv')
            confound_json = confound_file.replace('.tsv', '.json')
            time_series, out_labels = extract_time_series(fmri_file, parc_file, confounds_file=confound_file, list_of_confounds=confound_json)
            correlation_matrix = compute_functional_connectivity(time_series)

            if args.no_ses:
                os.makedirs(os.path.join(args.output_dir, 'tvb'), exist_ok=True)
                output_file = os.path.join(args.output_dir, f'{args.subject_id}_run-{index+1}_fc_matrix.npy')
                np.save(output_file, correlation_matrix)
                json.dump(out_labels, open(output_file.replace('.npy', '.json'), 'w'))
                np.savetxt(os.path.join(args.output_dir, 'tvb', f'{args.subject_id}_run-{index+1}_FC_ts.txt'), time_series, delimiter=' ', fmt='%f')
                np.savetxt(os.path.join(args.output_dir, 'tvb', f'{args.subject_id}_run-{index+1}_FC.csv'), correlation_matrix, delimiter=',', fmt='%f')
            else:
                os.makedirs(os.path.join(args.output_dir, session, 'tvb'), exist_ok=True)
                output_file = os.path.join(args.output_dir, session, f'{args.subject_id}_{session}_run-{index+1}_fc_matrix.npy')
                np.save(output_file, correlation_matrix)
                json.dump(out_labels, open(output_file.replace('.npy', '.json'), 'w'))
                np.savetxt(os.path.join(args.output_dir, session, 'tvb', f'{args.subject_id}_{session}_run-{index+1}_FC_ts.txt'), time_series, delimiter=' ', fmt='%f')
                np.savetxt(os.path.join(args.output_dir, session, 'tvb', f'{args.subject_id}_{session}_run-{index+1}_FC.csv'), correlation_matrix, delimiter=',', fmt='%f')

        session_pbar.close()
        print(f'Functional connectivity matrices saved to {args.output_dir}')
        
if __name__ == "__main__":
    main()