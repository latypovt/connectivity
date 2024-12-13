#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 --tractoflow_dir <path_to_tractoflow> --fmriprep_dir  --freesurfer_subjects_dir <path_to_freesurfer_subjects> --output_folder <path_to_output_folder>"
    exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --tractoflow_dir) TRACTOFLOW_DIR="$2"; shift ;;
        --fmriprep_dir) FMRI_DIR="$2"; shift ;;
        --freesurfer_subjects_dir) FSS_DIR="$2"; shift ;;
        --output_folder) OUTPUT_FOLDER="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if all required arguments are provided
if [ -z "$TRACTOFLOW_DIR" ] || [ -z "$FSS_DIR" ] || [ -z "$OUTPUT_FOLDER" ] || [ -z "$FMRI_DIR" ]; then
    usage
fi

# Load necessary modules
module load mrtrix3src

export SINGULARITYENV_BINDPATH=${PWD},${TRACTOFLOW_DIR},${OUTPUT_FOLDER},${FSS_DIR},${MYDIR},${FMRI_DIR}

CWDIR=${PWD}

#create empty list as variable
FAILED_SUBJECTS=()


# Loop through each subject in the subject folder
for SESSION in $(ls $TRACTOFLOW_DIR); do
    SUBJECT_DIR="${TRACTOFLOW_DIR}/${SESSION}"

    #check that folders exist in all 3 directories
    if [ -d "$SUBJECT_DIR" -a -d "${FMRI_DIR}/${SESSION}" -a -d "${FSS_DIR}/${SESSION}" ]; then
        cd $CWDIR
        
        IFS='_' read -r -a array <<< "$SESSION"
        SUBJECT="${array[0]}"
        SES="${array[1]}"
        echo "Processing: $SES of $SUBJECT"
        
        # Define paths
        MGZ_ASEG="${FSS_DIR}/${SUBJECT}/mri/aparc+aseg.mgz"
        ASEG_FILE="${FSS_DIR}/${SUBJECT}/mri/aparc_aseg.nii.gz"
        LUT_FILE="${MYDIR}/apps/connectivity/FreeSurferColorLUT.txt"
        FS_A2009S_FILE="${MYDIR}/apps/connectivity/fs_default.txt"
        CONVERTED_LABEL_FILE="${SUBJECT}_${SES}_DK_T1space.nii.gz"
        T1_FILE="${FSS_DIR}/${SUBJECT}/mri/brain.mgz"
        OUTPUT_FILE="${SUBJECT}_${SES}_DK_DiffusionSpace.nii.gz"
        TARGET_SPACE_FILE="${SUBJECT_DIR}/Register_T1/${SESSION}__t1_warped.nii.gz"
        if [ ! -f "$TARGET_SPACE_FILE" ]; then
            echo "Target space file not found: $TARGET_SPACE_FILE"
            continue
        fi

        if [ ! -d "${OUTPUT_FOLDER}/${SUBJECT}/${SES}" ]; then
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT}/${SES}/dwi"
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT}/${SES}/func"
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT}/${SES}/parc"
        fi
        # Convert labels using mrtrix3
        mrconvert ${MGZ_ASEG} $ASEG_FILE --force -quiet
        labelconvert $ASEG_FILE $LUT_FILE $FS_A2009S_FILE ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/parc/$CONVERTED_LABEL_FILE --force -quiet
        python ${MYDIR}/apps/connectivity/ANTsReg.py --fixed_image $TARGET_SPACE_FILE --moving_image $T1_FILE --apply_transform ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/parc/$CONVERTED_LABEL_FILE --output_transformed ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/parc/$OUTPUT_FILE --is_label
        
        if [ ! -d ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/dwi ]; then
            mkdir -p ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/dwi
        fi

        if [ ! -d ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/func ]; then
            mkdir -p ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/func
        fi
        cd ${OUTPUT_FOLDER}
        TRACKS_PATH=${TRACTOFLOW_DIR}/${SESSION}/PFT_Tracking/${SUBJECT}_${SES}__pft_tracking_prob_wm_seed_0.trk 
        FMRI_PATH=${FMRI_DIR}/${SUBJECT}/${SES}/func/${SUBJECT}_${SES}_task-rest_space-T1w_desc-preproc_bold.nii.gz
        FMRI_CONFOUNDS_PATH=${FMRI_DIR}/${SUBJECT}/${SES}/func/${SUBJECT}_${SES}_task-rest_desc-confounds_timeseries.tsv

        if [ -f "$TRACKS_PATH" ] & [ ! -f "$FMRI_PATH" ] & [ ! -f "$FMIR_CONFOUNDS_PATH" ]; then
            cp $TRACKS_PATH ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/dwi/${SUBJECT}_${SES}_run-1__pft_tracking_prob_wm_seed_0.trk
            cp $FMRI_PATH ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/func/${SUBJECT}_${SES}_task-rest_space-T1w_desc-preproc_bold.nii.gz
            cp $FMRI_CONFOUNDS_PATH ${OUTPUT_FOLDER}/${SUBJECT}/${SES}/func/${SUBJECT}_${SES}_task-rest_desc-confounds_timeseries.tsv

            python  ${MYDIR}/apps/connectivity/functional_connectivity.py $SUBJECT $SUBJECT
            python  ${MYDIR}/apps/connectivity/structural_connectivity.py $SUBJECT $SUBJECT
        else
            echo "Required files not found: $TRACKS_PATH, $FMRI_PATH, $FMRI_CONFOUNDS_PATH"
            # append $SESSION to the FAILED_SUBJECTS 
            FAILED_SUBJECTS+=($SESSION)
            continue
        fi
        echo "Processing completed for : $SESSION"
    else
        echo "Required directories not found: $SUBJECT_DIR, ${FMRI_DIR}/${SESSION}, ${FSS_DIR}/${SESSION}"
        # append $SESSION to the FAILED_SUBJECTS 
        FAILED_SUBJECTS+=($SESSION)
    fi
done

echo "Done."
echo "Failed subjects with missing FMRI/DTI: ${FAILED_SUBJECTS[@]}"