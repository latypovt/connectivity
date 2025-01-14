import argparse
import nibabel as nib
import numpy as np
from dipy.tracking.utils import connectivity_matrix, path_length, density_map, length
from dipy.io.image import load_nifti
from tqdm import tqdm
import os
from nilearn.input_data import NiftiLabelsMasker

def load_tractography(trk_file):
    """
    Load tractography data from a TRK file.

    Args:
    - trk_file (str): Path to the TRK file

    Returns:
    - streamlines (list): A list of streamlines (each streamline is a 2D numpy array of shape (N, 3))
    """
    tractogram = nib.streamlines.load(trk_file)
    streamlines = tractogram.streamlines
    return streamlines

def load_parcellation(parc_file):
    parc_data, affine = load_nifti(parc_file)
    return parc_data, affine

def list_sessions(subject_dir):
    list = [d for d in os.listdir(subject_dir) if os.path.isdir(os.path.join(subject_dir, d)) and d.startswith('ses-')]

    return list

def list_runs(session_dir):
    return [f for f in os.listdir(session_dir) if f.endswith('tracking_prob_wm_seed_0.trk')]

def compute_seeds_per_region(parc_data, seeds_per_voxel=10):
    """
    Computes the total number of seeds per region in a parcellation map.
    
    Args:
    - parcellation_map (str): Path to the parcellation map (e.g., 'aparc+aseg.nii.gz')
    - seeds_per_voxel (int): Number of seeds per voxel (default is 10)
    
    Returns:
    - seeds_per_region (list): A list containing the total number of seeds per region.
    """
    seg_data = nib.load(parc_data).get_fdata()
    labels = np.unique(seg_data)
    # Initialize a list to store the number of seeds per region
    seeds_per_region = []

    # Iterate through each region label and compute the number of seeds for that region
    for label in labels:
        if label == 0:
            continue  # Skip background (label 0)
        
        # Count the number of voxels belonging to this region
        num_voxels = np.sum(seg_data == label)
        
        # Compute the total number of seeds for this region
        total_seeds = num_voxels * seeds_per_voxel
        
        # Append the result to the list
        seeds_per_region.append(total_seeds)
    
    return seeds_per_region

# Example usage:
parcellation_map_path = 'aparc+aseg.nii.gz'  # Path to your parcellation map
seeds_per_region = compute_seeds_per_region(parcellation_map_path)

print("Seeds per region:", seeds_per_region)



def compute_connectivity_matrix(streamlines, parc_data, affine):
    """
    Computes the structural connectivity matrix from a set of streamlines and a parcellation map.

    Args:
    - streamlines (list): A list of streamlines (each streamline is a 2D numpy array of shape (N, 3))
    - parc_data (str): Path to the parcellation map (e.g., 'aparc+aseg.nii.gz')
    - affine (numpy.ndarray): Affine matrix of the parcellation map

    Returns:
    - conn_matrix (numpy.ndarray): A 2D numpy array representing the structural connectivity matrix
    - group (dict): A dictionary containing the mapping of each streamline to a pair of regions
    """

    # Check if the affine matrix is square
    #print unqiue values in parc_data
    parc_data = parc_data.astype(int)

    conn_matrix, group = connectivity_matrix(streamlines, affine, parc_data, return_mapping=True, mapping_as_streamlines=True)
    
    # Initialize matrix to store mean streamline lengths
    mean_length_matrix = np.zeros_like(conn_matrix, dtype=float)
    # Calculate mean streamline length for each pair of regions
    bundle_pbar = tqdm(zip(group.keys(), group.values()), desc='Streamline length', total=len(group.keys()), leave=False)
    for (i, j), bundle in bundle_pbar:
        mean_length_matrix[i, j] = np.mean(list(length(bundle)))
        mean_length_matrix[j, i] = mean_length_matrix[i, j]
    bundle_pbar.close()

    conn_matrix = conn_matrix[1:,1:]
    mean_length_matrix = mean_length_matrix[1:,1:]
    # drop first 3 rows and last 3 columns from the matrix array

    return conn_matrix, group, mean_length_matrix




def main():
    parser = argparse.ArgumentParser(description="Compute structural connectivity matrix from tractography and parcellation.")
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
            session_dir = os.path.join(subject_dir, 'dwi')
        else:
            session_dir = os.path.join(subject_dir, session, 'dwi')
        runs = list_runs(session_dir)
        session_pbar = tqdm(enumerate(runs), desc=f'Structural  {subject_dir, session}', total=len(runs), leave=False)
        for index, run in session_pbar:
            streamline_file = os.path.join(session_dir, run)
                # Load tractography data
            streamlines = load_tractography(streamline_file)

            parc_file = [d for d in os.listdir(os.path.join(subject_dir,'parc')) if d.endswith('_DK_DiffusionSpace.nii.gz')]
            parc_file = os.path.join(subject_dir, 'parc', parc_file[0])
            parc_data, affine = load_parcellation(parc_file)
            conn_matrix, group, mean_length = compute_connectivity_matrix(streamlines, parc_data, affine)

            if args.no_ses:
                output_file = os.path.join(args.output_dir, f'{args.subject_id}_run-{index+1}_sc_matrix.npy')
                output_file2 = os.path.join(args.output_dir, f'{args.subject_id}_run-{index+1}_length.npy')
                np.save(output_file, conn_matrix)
                np.save(output_file2, mean_length)
                os.makedirs(os.path.join(args.output_dir, 'tvb'), exist_ok=True)
                np.savetxt(os.path.join(args.output_dir, 'tvb', f'{args.subject_id}_run-{index+1}_SC.csv'), conn_matrix, delimiter=' ', fmt='%f')
            else:
                output_file = os.path.join(args.output_dir, session,  f'{args.subject_id}_{session}_run-{index+1}_sc_matrix.npy')
                output_file2 = os.path.join(args.output_dir, session, f'{args.subject_id}_{session}_run-{index+1}_length.npy')
                np.save(output_file, conn_matrix)
                np.save(output_file2, mean_length)
                os.makedirs(os.path.join(args.output_dir, session, 'tvb'), exist_ok=True)
                np.savetxt(os.path.join(args.output_dir, session, 'tvb', f'{args.subject_id}_{session}_run-{index+1}_SC.csv'), conn_matrix, delimiter=' ', fmt='%f')


        session_pbar.close()
        print(f'Structural connectivity matrices saved to {args.output_dir}')




if __name__ == "__main__":
    main()