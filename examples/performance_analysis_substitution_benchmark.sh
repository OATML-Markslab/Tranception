source activate tranception_env

# Replace the following variables based on the paths where you store models and data
export DMS_data_folder=$PATH_TO_DMS_SUBSTITUTIONS
export input_scoring_files_folder=$PATH_TO_MODEL_SCORES
export output_performance_file_folder=$PATH_TO_OUTPUT_PERFORMANCE

export model_list="all_models" 
export DMS_reference_file_path="../proteingym/ProteinGym_reference_file_substitutions.csv"

python3 ../performance_analysis_proteingym.py \
                --model_list ${model_list} \
                --input_scoring_files_folder ${input_scoring_files_folder} \
                --output_performance_file_folder ${output_performance_file_folder} \
                --DMS_reference_file_path ${DMS_reference_file_path} \
                --DMS_data_folder ${DMS_data_folder}