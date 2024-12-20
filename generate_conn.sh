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
GMI=/d/gmi/1

export SINGULARITYENV_BINDPATH=${MYDIR},${OUTPUT_FOLDER},${GMI}

CWDIR=${PWD}

#create empty list as variable
FAILED_SUBJECTS=()
FAILED_FMRI=()
FAILED_TRACTS=()

# Function to process a single subject
process_subject() {
    SESSION=$1
    SUBJECT_DIR="${TRACTOFLOW_DIR}/${SESSION}"
    IFS='_' read -r -a array <<< "$SESSION"
    SUBJECT="${array[0]}"
    SES="${array[1]}"
    SUBJECT_SES="${array[0]}"
    echo "Processing: $SUBJECT"

    #check that folders exist in all 3 directories
    if [ -d $SUBJECT_DIR -a -d ${FMRI_DIR}/${SUBJECT} -a -d ${FSS_DIR}/${SUBJECT} ]; then

        cd $CWDIR
        #check if session is empty
        if [ -z "$SES" ]; then
            CONVERTED_LABEL_FILE="${SUBJECT}_DK_T1space.nii.gz"
            OUTPUT_FILE="${SUBJECT}_DK_DiffusionSpace.nii.gz"
            SUB_SES_FILENAME="${SUBJECT}"

        else
            SUBJECT_SES="${SUBJECT}/${SES}"
            CONVERTED_LABEL_FILE="${SUB_SES_FILENAME}_DK_T1space.nii.gz"
            OUTPUT_FILE="${SUB_SES_FILENAME}_DK_DiffusionSpace.nii.gz"
            SUB_SES_FILENAME="${SUB_SES_FILENAME}"
        fi    
        
        # Define paths
        MGZ_ASEG="${FSS_DIR}/${SUBJECT}/mri/aparc+aseg.mgz"
        if [ ! -f "$MGZ_ASEG" ]; then
            echo "FreeSurfer segmentation file not found: $MGZ_ASEG - perhaps recon-all was not finished?"
            FAILED_SUBJECTS+=($SESSION)
            return
        fi
        ASEG_FILE="${OUTPUT_FOLDER}/${SUBJECT_SES}/anat/aparc_aseg.nii.gz"
        LUT_FILE="${MYDIR}/apps/connectivity/FreeSurferColorLUT.txt"
        FS_A2009S_FILE="${MYDIR}/apps/connectivity/fs_default.txt"
        T1_FILE="${FSS_DIR}/${SUBJECT}/mri/brain.mgz"
        TARGET_SPACE_FILE="${SUBJECT_DIR}/Register_T1/${SESSION}__t1_warped.nii.gz"
        if [ ! -f "$TARGET_SPACE_FILE" ]; then
            echo "Target space file not found: $TARGET_SPACE_FILE"
            return
        fi

        if [ ! -d ${OUTPUT_FOLDER}/${SUBJECT_SES} ]; then
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT_SES}/anat"
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT_SES}/dwi"
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT_SES}/func"
            mkdir -p "${OUTPUT_FOLDER}/${SUBJECT_SES}/parc"
        fi

        # Generating input data for connectivity analysis - masks, parcellations, etc.
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
        if [ ! -f ${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_space-T1w_desc-preproc_bold.nii.gz ] && [ ! -f ${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-1_space-T1w_desc-preproc_bold.nii.gz ]; then
            echo "fMRI data exists, but BOLD NIFTI is missing: ${FMRI_DIR}/${SUBJECT_SES} - perhaps fmriprep was not finished?" 
            FAILED_FMRI+=($SESSION)
            return
        fi

        RUNS_FMRI=$(ls ${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_*space-T1w_desc-preproc_bold.nii.gz | wc -l)
        if [ $RUNS_FMRI -gt 1 ]; then

            #loop through each run
            for RUN in $(seq 1 $RUNS_FMRI); do
                FMRI_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-${RUN}_space-T1w_desc-preproc_bold.nii.gz
                FMRI_CONFOUNDS_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-${RUN}_desc-confounds_timeseries.tsv
                

                if [ -f "$FMRI_PATH" ] && [ -f "$FMRI_CONFOUNDS_PATH" ]; then
                    cp $FMRI_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-${RUN}_space-T1w_desc-preproc_bold.nii.gz
                    cp $FMRI_CONFOUNDS_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-${RUN}_desc-confounds_timeseries.tsv
                    FMRI_RT=$(jq -r '.RepetitionTime' ${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-${RUN}_space-T1w_desc-preproc_bold.json)
                else
                    echo "fMRI files not found: $FMRI_PATH, $FMRI_CONFOUNDS_PATH"
                    # append $SESSION to the FAILED_SUBJECTS 
                    FAILED_FMRI+=($SESSION)
                    return
                fi
            done

        elif [ $RUNS_FMRI -eq 1 ]; then
            FMRI_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_space-T1w_desc-preproc_bold.nii.gz
            FMRI_CONFOUNDS_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_desc-confounds_timeseries.tsv
            

            if [ ! -f $FMRI_PATH ]; then
                FMRI_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-1_space-T1w_desc-preproc_bold.nii.gz
                FMRI_CONFOUNDS_PATH=${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-1_desc-confounds_timeseries.tsv
                FMRI_RT=$(jq -r '.RepetitionTime' ${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_run-1_space-T1w_desc-preproc_bold.json)
            else
                FMRI_RT=$(jq -r '.RepetitionTime' ${FMRI_DIR}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_space-T1w_desc-preproc_bold.json)
            fi 

            if [ -f "$FMRI_PATH" ] && [ -f "$FMRI_CONFOUNDS_PATH" ]; then
                cp $FMRI_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_space-T1w_desc-preproc_bold.nii.gz
                cp $FMRI_CONFOUNDS_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/func/${SUB_SES_FILENAME}_task-rest_desc-confounds_timeseries.tsv

            else
                echo "fMRI files not found: $FMRI_PATH, $FMRI_CONFOUNDS_PATH"
                FAILED_FMRI+=($SESSION)
                return
            fi   
        else
            echo "fMRI data is missing or not found:  ${FMRI_DIR}/${SESSION}"
            # append $SESSION to the FAILED_SUBJECTS 
            FAILED_FMRI+=($SESSION)
            return
        fi
                #check number of runs in the tractoflow folder (all runs are stored in different folders)
        RUNS_TRACTS=$(ls -d ${TRACTOFLOW_DIR}| wc -l)
        if [ $RUNS_TRACTS -gt 1 ]; then
            #loop through each run
            for RUN in $(seq 1 $RUNS_TRACTS); do

                TRACKS_PATH=${TRACTOFLOW_DIR}/${SUB_SES_FILENAME}_run-${RUN}/PFT_Tracking/${SUB_SES_FILENAME}_run-${RUN}__pft_tracking_prob_wm_seed_0.trk #__pft_tracking_prob_interface_seed_0.trk #

                if [ -f "$TRACKS_PATH" ] ; then
                    cp $TRACKS_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/dwi/${SUB_SES_FILENAME}_run-${RUN}__pft_tracking_prob_wm_seed_0.trk
                else
                    echo "Tracts are missing or not found: $TRACKS_PATH"
                    # append $SESSION to the FAILED_SUBJECTS 
                    FAILED_TRACTS+=($SESSION)
                    return
                fi
            done


        elif [ $RUNS_TRACTS -eq 1 ]; then
            TRACKS_PATH=${TRACTOFLOW_DIR}/${SESSION}/PFT_Tracking/${SUB_SES_FILENAME}__pft_tracking_prob_wm_seed_0.trk #__pft_tracking_prob_interface_seed_0.trk #
            

            if [ -f "$TRACKS_PATH" ]; then
                cp $TRACKS_PATH ${OUTPUT_FOLDER}/${SUBJECT_SES}/dwi/${SUB_SES_FILENAME}_run-1__pft_tracking_prob_wm_seed_0.trk
                
            else
                echo "Tracts are missing or not found: $TRACKS_PATH"
                # append $SESSION to the FAILED_SUBJECTS 
                FAILED_TRACTS+=($SESSION)
                return
            fi

        else
            echo "Tractography data is missing: ${TRACTOFLOW_DIR}/${SUB_SES_FILENAME}"
            # append $SESSION to the FAILED_SUBJECTS 
            FAILED_TRACTS+=($SESSION)
        fi

        if [ -z "$SES" ]; then
            python  ${MYDIR}/apps/connectivity/structural_connectivity.py $SUBJECT $SUBJECT --no_ses
            python  ${MYDIR}/apps/connectivity/functional_connectivity.py $SUBJECT $SUBJECT --no_ses --repetition_time $FMRI_RT

        else
            python  ${MYDIR}/apps/connectivity/structural_connectivity.py $SUBJECT $SUBJECT
            python  ${MYDIR}/apps/connectivity/functional_connectivity.py $SUBJECT $SUBJECT --repetition_time $FMRI_RT
        fi

    elif [ ! -d $SUBJECT_DIR ]; then
        echo "Error while processing $SESSION"
        echo "Tractoflow output is not found: $SUBJECT_DIR. Failed to find preprocessed data"
        # append $SESSION to the FAILED_SUBJECTS 
        FAILED_SUBJECTS+=($SESSION)
        return
    elif [ ! -d ${FMRI_DIR}/${SUBJECT} ]; then
        echo "Error while processing $SESSION"
        echo "fMRI data is missing or not found:  ${FMRI_DIR}/${SESSION}. Failed to find preprocessed data"
        # append $SESSION to the FAILED_SUBJECTS 
        FAILED_SUBJECTS+=($SESSION)
        return
    elif [ ! -d ${FSS_DIR}/${SUBJECT}* ]; then
        echo "Error while processing $SESSION"
        echo "FreeSurfer output is not found: ${FSS_DIR}/${SESSION}. Failed to find preprocessed data"
        # append $SESSION to the FAILED_SUBJECTS 
        FAILED_SUBJECTS+=($SESSION)
        return
    else
        echo "Error while processing $SESSION"
        echo "Something went wrong. Failed to find preprocessed data"
        # append $SESSION to the FAILED_SUBJECTS 
        FAILED_SUBJECTS+=($SESSION)
        return
    fi
}

# Loop through each subject in the subject folder and process them in parallel
for SESSION in $(ls $TRACTOFLOW_DIR); do
    process_subject $SESSION &
done

# Wait for all background processes to finish
wait

echo "Done."
echo "Failed subjects with missing data: ${FAILED_SUBJECTS[@]}"
echo "Failed subjects with missing fMRI data: ${FAILED_FMRI[@]}"
echo "Failed subjects with missing tractograms: ${FAILED_TRACTS[@]}"