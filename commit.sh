#!/bin/bash

# Function to display usage
usage() {
  echo "Usage: $0 -d TRACTOFLOW_DIR -t NUM_THREADS"
  exit 1
}

# Parse command line arguments
while getopts ":d:t:" opt; do
  case ${opt} in
    d )
      TRACTOFLOW_DIR=$OPTARG
      ;;
    t )
      NUM_THREADS=$OPTARG
      ;;
    \? )
      usage
      ;;
  esac
done

# Check if required arguments are provided
if [ -z "${TRACTOFLOW_DIR}" ] || [ -z "${NUM_THREADS}" ]; then
  usage
fi

# Function to process a single subject
process_subject() {
  SUBJECT_DIR=$1
  SUBJECT_ID=$(basename ${SUBJECT_DIR})

  # Paths for subject-specific files
  TRACK_FILE="${SUBJECT_DIR}/PFT_Tracking/${SUBJECT_ID}__pft_tracking_prob_wm_seed_0.trk"
  DWI_FILE="${SUBJECT_DIR}/Resample_DWI/${SUBJECT_ID}__dwi_resampled.nii.gz"
  BVAL_FILE="${SUBJECT_DIR}/Extract_FODF_Shell/${SUBJECT_ID}__bval_fodf"
  BVEC_FILE="${SUBJECT_DIR}/Extract_FODF_Shell/${SUBJECT_ID}__bvec_fodf"
  OUT_DIR="${SUBJECT_DIR}/COMMIT"
  PEAKS_FILE="${SUBJECT_DIR}/FODF_Metrics/${SUBJECT_ID}__peaks.nii.gz"

  # Create the output directory if it does not exist
  mkdir -p ${OUT_DIR}

  export SINGULARITYENV_BINDPATH=${PWD}

  # Run the command inside the singularity container
  singularity run --bind=/d/gmi/1/timurlatypov:/d/gmi/1/timurlatypov /d/gmi/1/timurlatypov/apps/scilus_1.6.0.sif \
  scil_run_commit.py ${TRACK_FILE} ${DWI_FILE} ${BVAL_FILE} ${BVEC_FILE} ${OUT_DIR} \
    --in_peaks ${PEAKS_FILE} --b_thr 50 --nbr_dir 500 --ball_stick \
    --para_diff "1.7E-3" --iso_diff "2.0E-3" -v > ${SUBJECT_DIR}/commit.log 2>&1

}

# Loop through each subject directory in the tractoflow folder and process them in parallel
for SUBJECT_DIR in ${TRACTOFLOW_DIR}/*; do
  while [ $(jobs -r | wc -l) -ge ${NUM_THREADS} ]; do
    sleep 1
  done
  echo "Processing subject: ${SUBJECT_DIR}"
  process_subject ${SUBJECT_DIR} &
done

# Wait for all background processes to finish
wait

echo "Done."