import pandas as pd
import numpy as np
import os
import argparse
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score, matthews_corrcoef
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def standardization(x):
    """Assumes input is numpy array or pandas series"""
    return (x - x.mean()) / x.std()

def main():
    parser = argparse.ArgumentParser(description='Tranception performance analysis')
    parser.add_argument('--model_list', default='tranception_only', type=str, help='Whether to compute the performance of tranception only Vs all models present in input scoring files Vs all scores that score all mutants (ie. exclude EVE, EVmutation, etc.) [tranception_only|all_models|all_mutants]')
    parser.add_argument('--input_scoring_files_folder', type=str, help='Name of folder where all input scores are present (expects one scoring file per DMS)')
    parser.add_argument('--output_performance_file_folder', default='./outputs/tranception_performance', type=str, help='Name of folder where to save performance analysis files')
    parser.add_argument('--DMS_reference_file_path', type=str, help='Reference file with list of DMSs to consider')
    parser.add_argument('--DMS_data_folder', type=str, help='Path to folder that contains all DMS datasets')
    parser.add_argument('--indel_mode', action='store_true', help='Whether to score sequences with insertions and deletions')
    args = parser.parse_args()
    
    mapping_protein_seq_DMS = pd.read_csv(args.DMS_reference_file_path)
    num_DMS=len(mapping_protein_seq_DMS)
    print("There are {} DMSs in mapping file".format(num_DMS))

    if not args.indel_mode:
        uniprot_Neff_lookup = mapping_protein_seq_DMS[['UniProt_ID','MSA_Neff_L_category']].drop_duplicates()
        uniprot_Neff_lookup.columns=['UniProt_ID','Neff_L_category']
        uniprot_taxon_lookup = mapping_protein_seq_DMS[['UniProt_ID','taxon']].drop_duplicates()
        uniprot_taxon_lookup.columns=['UniProt_ID','Taxon']

    if args.model_list=="tranception_only":
        score_variables = ['Tranception']
    elif args.model_list=="all_models":
        score_file = pd.read_csv(args.input_scoring_files_folder+os.sep+mapping_protein_seq_DMS["DMS_filename"].values[0])
        score_variables = [ x for x in score_file.columns if x not in ['DMS_score','DMS_score_bin','mutant']]
    elif args.model_list=="all_mutants":
        score_variables = ['Tranception_no_retrieval','Tranception_retrieval','MSA_Transformer_single','MSA_Transformer_ensemble','ESM1v_single','ESM1v_ensemble','Wavenet','RITA_s','RITA_m','RITA_l','RITA_xl']

    if not os.path.isdir(args.output_performance_file_folder):
        os.mkdir(args.output_performance_file_folder)
    performance_all_DMS={}
    output_filename={}
    for metric in ['Spearman','AUC','MCC']:
        performance_all_DMS[metric]={}
        output_filename[metric]="performance_"+metric+"_"+args.model_list+".csv"
        for i, score in enumerate(score_variables):
            performance_all_DMS[metric][score]=i
            if not args.indel_mode:
                for depth in ['1','2','3','4','5+']:
                    performance_all_DMS[metric][score+'_'+depth] = i
        performance_all_DMS[metric]['number_mutants']=-1
        if not args.indel_mode:
            performance_all_DMS[metric]['UniProt_ID']=-1
            performance_all_DMS[metric]['Neff_L_category']=-1
            performance_all_DMS[metric]['Taxon']=-1
        performance_all_DMS[metric]=pd.DataFrame.from_dict(performance_all_DMS[metric],orient='index').reset_index()
        performance_all_DMS[metric].columns=['score','score_index']

    list_DMS = mapping_protein_seq_DMS["DMS_id"]

    for DMS_id in list_DMS:
        try:
            print(DMS_id)    
            UniProt_ID = mapping_protein_seq_DMS["UniProt_ID"][mapping_protein_seq_DMS["DMS_id"]==DMS_id].values[0]
            DMS_binarization_cutoff_ProteinGym = mapping_protein_seq_DMS["DMS_binarization_cutoff"][mapping_protein_seq_DMS["DMS_id"]==DMS_id].values[0]
            DMS_filename = mapping_protein_seq_DMS["DMS_filename"][mapping_protein_seq_DMS["DMS_id"]==DMS_id].values[0]
            if not args.indel_mode:
                Neff_L_category	= mapping_protein_seq_DMS["MSA_Neff_L_category"][mapping_protein_seq_DMS["DMS_id"]==DMS_id].values[0]
                Taxon = mapping_protein_seq_DMS["taxon"][mapping_protein_seq_DMS["DMS_id"]==DMS_id].values[0]

            DMS_file = pd.read_csv(args.DMS_data_folder+os.sep+DMS_filename)
            print("Length DMS: {}".format(len(DMS_file)))
    
            if args.model_list=="tranception_only":
                tranception = pd.read_csv(args.input_scoring_files_folder + os.sep + DMS_id + ".csv")
                tranception = tranception[['mutant','avg_score']]
                tranception.columns=['mutant','Tranception']    
                merged_scores = pd.merge(DMS_file, tranception, on='mutant', how='inner')
                merged_scores.dropna(inplace=True)
            elif args.model_list=="all_models" or args.model_list=="all_mutants":
                merged_scores = pd.read_csv(args.input_scoring_files_folder + os.sep + DMS_id + ".csv")
                if args.model_list=="all_mutants": merged_scores=merged_scores[['mutant','DMS_score','DMS_score_bin']+score_variables]
                merged_scores.dropna(inplace=True)
        except:
            print("At least one scoring file missing")
            continue

        if not args.indel_mode:
            merged_scores['mutation_depth']=merged_scores['mutant'].apply(lambda x: len(x.split(":")))
            merged_scores['mutation_depth_grouped']=merged_scores['mutation_depth'].apply(lambda x: '5+' if x >=5 else str(x))

        performance_DMS = {}
        for metric in ['Spearman','AUC','MCC']:
            performance_DMS[metric]={}
        for score in score_variables:
            performance_DMS['Spearman'][score] = spearmanr(merged_scores['DMS_score'], merged_scores[score])[0]
            try:
                performance_DMS['AUC'][score] = roc_auc_score(y_true=merged_scores['DMS_score_bin'], y_score=merged_scores[score])
            except:
                print("AUC issue with {}".format(DMS_id))
                performance_DMS['AUC'][score] = np.nan
            try:
                median_cutoff=merged_scores[score].median()
                merged_scores[score+"_bin"]=merged_scores[score].map(lambda x: 1 if x >= median_cutoff else 0)
                performance_DMS['MCC'][score] = matthews_corrcoef(y_true=merged_scores['DMS_score_bin'], y_pred=merged_scores[score+"_bin"])
            except:
                print("MCC issue with {}".print(DMS_id))
                performance_DMS['MCC'][score] = np.nan
        
        if not args.indel_mode:
            for score in score_variables:
                for depth in ['1','2','3','4','5+']:
                    merged_scores_depth = merged_scores[merged_scores.mutation_depth_grouped==depth]
                    if len(merged_scores_depth) > 0:
                        performance_DMS['Spearman'][score+'_'+depth] = spearmanr(merged_scores_depth['DMS_score'], merged_scores_depth[score])[0]
                        try:
                            performance_DMS['AUC'][score+'_'+depth] = roc_auc_score(y_true=merged_scores_depth['DMS_score_bin'], y_score=merged_scores_depth[score])
                        except:
                            performance_DMS['AUC'][score+'_'+depth] = np.nan
                        try:
                            performance_DMS['MCC'][score+'_'+depth] = matthews_corrcoef(y_true=merged_scores_depth['DMS_score_bin'], y_pred=merged_scores_depth[score+"_bin"])
                        except:
                            performance_DMS['MCC'][score+'_'+depth] = np.nan
                    else:
                        performance_DMS['Spearman'][score+'_'+depth] = np.nan
                        performance_DMS['AUC'][score+'_'+depth] = np.nan
                        performance_DMS['MCC'][score+'_'+depth] = np.nan

        print("Number of mutants: {}".format(len(merged_scores['DMS_score'].values)))
        for metric in ['Spearman','AUC','MCC']:
            performance_DMS[metric]['number_mutants']=len(merged_scores['DMS_score'].values)
            if not args.indel_mode:
                performance_DMS[metric]['UniProt_ID'] = UniProt_ID
                performance_DMS[metric]['Neff_L_category'] = Neff_L_category
                performance_DMS[metric]['Taxon'] = Taxon
            performance_DMS[metric] = pd.DataFrame.from_dict(performance_DMS[metric],orient='index').reset_index()
            performance_DMS[metric].columns=['score',DMS_id]
        
            performance_all_DMS[metric]=pd.merge(performance_all_DMS[metric],performance_DMS[metric],on='score',how='left')
    
    for metric in ['Spearman','AUC','MCC']:
        performance_all_DMS[metric]=performance_all_DMS[metric].set_index('score')
        del performance_all_DMS[metric]['score_index']
        performance_all_DMS[metric]=performance_all_DMS[metric].transpose()
        performance_all_DMS[metric].loc['Average'] = performance_all_DMS[metric].mean()
        performance_all_DMS[metric].to_csv(args.output_performance_file_folder + os.sep + output_filename[metric])
    
        if not args.indel_mode:
            uniprot_metric_performance = performance_all_DMS[metric].groupby(['UniProt_ID']).mean().reset_index()
            uniprot_metric_performance = pd.merge(uniprot_metric_performance,uniprot_Neff_lookup,on='UniProt_ID', how='left')
            uniprot_metric_performance = pd.merge(uniprot_metric_performance,uniprot_taxon_lookup,on='UniProt_ID', how='left')
            uniprot_metric_performance.loc['Average'] = uniprot_metric_performance.mean()
            uniprot_metric_performance.to_csv(args.output_performance_file_folder + os.sep + 'Uniprot_aggregation_'+output_filename[metric], index=False)
            uniprot_metric_performance.groupby(['Neff_L_category']).mean().to_csv(args.output_performance_file_folder + os.sep + 'Uniprot_Neff_aggregation_'+output_filename[metric])
            uniprot_metric_performance.groupby(['Taxon']).mean().to_csv(args.output_performance_file_folder + os.sep + 'Uniprot_Taxon_aggregation_'+output_filename[metric])

if __name__ == '__main__':
    main()