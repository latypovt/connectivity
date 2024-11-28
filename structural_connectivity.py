import argparse
import nibabel as nib
import numpy as np
from dipy.tracking.utils import connectivity_matrix
from dipy.tracking.streamline import transform_streamlines
from dipy.tracking.streamline import Streamlines
from dipy.io.image import load_nifti
from tqdm import tqdm
import os
from nilearn.input_data import NiftiLabelsMasker

def load_tractography(trk_file):
    tractogram = nib.streamlines.load(trk_file)
    streamlines = tractogram.streamlines
    return streamlines

def load_parcellation(parc_file):
    parc_data, affine = load_nifti(parc_file)
    return parc_data, affine

def list_sessions(subject_dir):
    return [d for d in os.listdir(subject_dir) if os.path.isdir(os.path.join(subject_dir, d)) and d.startswith('ses-')]

def list_runs(session_dir):
    return [f for f in os.listdir(session_dir) if f.endswith('tracking_prob_wm_seed_0.trk')]

def compute_connectivity_matrix(streamlines, parc_data, affine):
    # Check if the affine matrix is square
    #print unqiue values in parc_data
    parc_data = parc_data.astype(int)

    conn_matrix, group = connectivity_matrix(streamlines, affine, parc_data, return_mapping=True, mapping_as_streamlines=False)
    conn_matrix = conn_matrix[1:,1:]

    # drop first 3 rows and last 3 columns from the matrix array

    return conn_matrix, group

def main():
    parser = argparse.ArgumentParser(description="Compute structural connectivity matrix from tractography and parcellation.")
    parser.add_argument('subject_id', type=str, help='Subject ID')
    parser.add_argument('output_dir', type=str, help='Directory to save the connectivity matrices')
    args = parser.parse_args()

    subject_dir =  args.subject_id 
    sessions = list_sessions(subject_dir)

    for ses_id, session in enumerate(sessions):
        session_dir = os.path.join(subject_dir, session, 'dwi')
        runs = list_runs(session_dir)
        for index, run in enumerate(tqdm(runs, desc=f'Structural  {subject_dir, session}')):
            streamline_file = os.path.join(session_dir, run)
                # Load tractography data
            streamlines = load_tractography(streamline_file)
            parc_file = streamline_file.replace('_run-1__pft_tracking_prob_wm_seed_0.trk', '_DK_DiffusionSpace.nii.gz')
            parc_file = parc_file.replace('dwi', 'parc')
            parc_data, affine = load_parcellation(parc_file)
            conn_matrix, group = compute_connectivity_matrix(streamlines, parc_data, affine)

           
            output_file = os.path.join(args.output_dir, f'{args.subject_id}_ses-{ses_id+1}_run-{index+1}_sc_matrix.npy')
            np.save(output_file, conn_matrix)




if __name__ == "__main__":
    main()