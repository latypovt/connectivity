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
MYDIR="/d/gmi/1/timurlatypov"
# Load necessary modules
module load mrtrix3src

export SINGULARITYENV_BINDPATH=${MYDIR},${OUTPUT_FOLDER}


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
        SUBJECT_SES="${array[0]}"
        echo "Processing: $SUBJECT"

        #check if session is empty
        if [ -z "$SES" ]; then
            echo "No session number found for $SUBJECT"
            CONVERTED_LABEL_FILE="${SUBJECT}_DK_T1space.nii.gz"
            OUTPUT_FILE="${SUBJECT}_DK_DiffusionSpace.nii.gz"
            SUB_SES_FILENAME="${SUBJECT}"

        else
            SUBJECT_SES="${SUBJECT}/ses-${SES}"
            CONVERTED_LABEL_FILE="${SUB_SES_FILENAME}_DK_T1space.nii.gz"
            OUTPUT_FILE="${SUB_SES_FILENAME}_DK_DiffusionSpace.nii.gz"
            SUB_SES_FILENAME="${SUB_SES_FILENAME}"
        fi    
        
        # Define paths
        MGZ_ASEG="${FSS_DIR}/${SUBJECT}/mri/aparc+aseg.mgz"
        ASEG_FILE="${OUTPUT_FOLDER}/${SUBJECT_SES}/anat/aparc_aseg.nii.gz"
        LUT_FILE="${MYDIR}/apps/connectivity/FreeSurferColorLUT.txt"
        FS_A2009S_FILE="${MYDIR}/apps/connectivity/fs_default.txt"
        T1_FILE="${FSS_DIR}/${SUBJECT}/mri/brain.mgz"
        TARGET_SPACE_FILE="${SUBJECT_DIR}/Register_T1/${SESSION}__t1_warped.nii.gz"
        if [ ! -f "$TARGET_SPACE_FILE" ]; then
            echo "Target space file not found: $TARGET_SPACE_FILE"
            continue
        fi

        if [ ! -d "${OUTPUT_FOLDER}/${SUBJECT_SES}" ]; then
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT_SES}/anat"
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT_SES}/dwi"
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT_SES}/func"
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT_SES}/parc"
        fi
        # Convert labels using mrtrix3
        cp $MGZ_ASEG ${OUTPUT_FOLDER}/${SUBJECT_SES}/anat/aparc_aseg.mgz
        mrconvert ${OUTPUT_FOLDER}/${SUBJECT_SES}/anat/aparc_aseg.mgz $ASEG_FILE --force -quiet

        labelconvert $ASEG_FILE $LUT_FILE $FS_A2009S_FILE ${OUTPUT_FOLDER}/${SUBJECT_SES}/parc/$CONVERTED_LABEL_FILE --force -quiet
        python ${MYDIR}/apps/connectivity/ANTsReg.py --fixed_image $TARGET_SPACE_FILE --moving_image $T1_FILE --apply_transform ${OUTPUT_FOLDER}/${SUBJECT_SES}/parc/$CONVERTED_LABEL_FILE --output_transformed ${OUTPUT_FOLDER}/${SUBJECT_SES}/parc/$OUTPUT_FILE --is_label
        
        if [ ! -d ${OUTPUT_FOLDER}/${SUBJECT_SES}/dwi ]; then
            mkdir -p ${OUTPUT_FOLDER}/${SUBJECT_SES}/dwi
        fi

        if [ ! -d ${OUTPUT_FOLDER}/${SUBJECT_SES}/func ]; then
            mkdir -p ${OUTPUT_FOLDER}/${SUBJECT_SES}/func
        fi
        cd ${OUTPUT_FOLDER}

        #check number of runs in the fmri folder (all runs are stored in the same folder)
        RUNS=$(ls ${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_*space-T1w_desc-preproc_bold.nii.gz | wc -l)
        if [ $RUNS -gt 1 ]; then
            echo "Multiple runs found for $SUBJECT"

            #loop through each run
            for RUN in $(seq 1 $RUNS); do
                FMRI_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-${RUN}_space-T1w_desc-preproc_bold.nii.gz
                FMRI_CONFOUNDS_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-${RUN}_desc-confounds_timeseries.tsv
                TRACKS_PATH=${TRACTOFLOW_DIR}/${SESSION}/PFT_Tracking/${SUB_SES_FILENAME}__pft_tracking_prob_wm_seed_0.trk

                if [ -f "$TRACKS_PATH" ] & [ -f "$FMRI_PATH" ] & [ -f "$FMRI_CONFOUNDS_PATH" ]; then
                    cp $TRACKS_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/dwi/${SUB_SES_FILENAME}_run-${RUN}__pft_tracking_prob_wm_seed_0.trk
                    cp $FMRI_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-${RUN}_space-T1w_desc-preproc_bold.nii.gz
                    cp $FMRI_CONFOUNDS_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-${RUN}_desc-confounds_timeseries.tsv
                else
                    echo "Required files not found: $TRACKS_PATH, $FMRI_PATH, $FMRI_CONFOUNDS_PATH"
                    # append $SESSION to the FAILED_SUBJECTS 
                    FAILED_SUBJECTS+=($SESSION)
                fi
            done
            echo "Copying and preprocessing completed for : $SESSION"
            if [ -z "$SES" ]; then
                python  ${MYDIR}/apps/connectivity/functional_connectivity.py $SUBJECT $SUBJECT --no_ses
                python  ${MYDIR}/apps/connectivity/structural_connectivity.py $SUBJECT $SUBJECT --no_ses
            else
                python  ${MYDIR}/apps/connectivity/functional_connectivity.py $SUBJECT $SUBJECT
                python  ${MYDIR}/apps/connectivity/structural_connectivity.py $SUBJECT $SUBJECT
            fi

        elif [ $RUNS -eq 1 ]; then
            echo "Single run found for $SUBJECT"
            TRACKS_PATH=${TRACTOFLOW_DIR}/${SESSION}/PFT_Tracking/${SUB_SES_FILENAME}__pft_tracking_prob_wm_seed_0.trk 
            FMRI_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_space-T1w_desc-preproc_bold.nii.gz
            FMRI_CONFOUNDS_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_desc-confounds_timeseries.tsv

            if [ -f "$TRACKS_PATH" ] & [ -f "$FMRI_PATH" ] & [ -f "$FMRI_CONFOUNDS_PATH" ]; then
                cp $TRACKS_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/dwi/${SUB_SES_FILENAME}_run-1__pft_tracking_prob_wm_seed_0.trk
                cp $FMRI_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_space-T1w_desc-preproc_bold.nii.gz
                cp $FMRI_CONFOUNDS_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_desc-confounds_timeseries.tsv

            else
                echo "Required files not found: $TRACKS_PATH, $FMRI_PATH, $FMRI_CONFOUNDS_PATH"
                # append $SESSION to the FAILED_SUBJECTS 
                FAILED_SUBJECTS+=($SESSION)
                continue
            fi
            echo "Copying and preprocessing completed for : $SESSION"
            
            if [ -z "$SES" ]; then
                python  ${MYDIR}/apps/connectivity/functional_connectivity.py $SUBJECT $SUBJECT --no_ses
                python  ${MYDIR}/apps/connectivity/structural_connectivity.py $SUBJECT $SUBJECT --no_ses
            else
                python  ${MYDIR}/apps/connectivity/functional_connectivity.py $SUBJECT $SUBJECT
                python  ${MYDIR}/apps/connectivity/structural_connectivity.py $SUBJECT $SUBJECT
            fi

        else
            echo "Required directories not found: $SUBJECT_DIR, ${FMRI_DIR}/${SESSION}, ${FSS_DIR}/${SESSION}"
            # append $SESSION to the FAILED_SUBJECTS 
            FAILED_SUBJECTS+=($SESSION)
        fi
    else
        echo "Required directories not found: $SUBJECT_DIR, ${FMRI_DIR}/${SESSION}, ${FSS_DIR}/${SESSION}"
        # append $SESSION to the FAILED_SUBJECTS 
        FAILED_SUBJECTS+=($SESSION)
    fi
done

echo "Done."
echo "Failed subjects with missing FMRI/DTI: ${FAILED_SUBJECTS[@]}"