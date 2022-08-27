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

def compute_bootstrap_standard_error(df, number_assay_reshuffle=10000):
    """
    Computes the non-parametric bootstrap standard error for the mean estimate of a given performance metric (eg., Spearman, AUC) across DMS assays (ie., the sample standard deviation of the mean across bootstrap samples)
    """
    model_names = df.columns
    mean_performance_across_samples = []
    for sample in range(number_assay_reshuffle):
        mean_performance_across_samples.append(df.sample(frac=1.0, replace=True).mean(axis=0)) #Resample a dataset of the same size (with replacement) then take the sample mean
    mean_performance_across_samples=pd.DataFrame(data=mean_performance_across_samples,columns=model_names)
    print(mean_performance_across_samples.head())
    return mean_performance_across_samples.std(ddof=1) #Unbiased estimate with ddof=1

def main():
    parser = argparse.ArgumentParser(description='Tranception performance analysis')
    parser.add_argument('--model_list', default='tranception_only', type=str, help='Whether to compute the performance of tranception only Vs all models present in input scoring files [tranception_only|all_models]')
    parser.add_argument('--input_scoring_files_folder', type=str, help='Name of folder where all input scores are present (expects one scoring file per DMS)')
    parser.add_argument('--output_performance_file_folder', default='./outputs/tranception_performance', type=str, help='Name of folder where to save performance analysis files')
    parser.add_argument('--DMS_reference_file_path', type=str, help='Reference file with list of DMSs to consider')
    parser.add_argument('--DMS_data_folder', type=str, help='Path to folder that contains all DMS datasets')
    parser.add_argument('--indel_mode', action='store_true', help='Whether to score sequences with insertions and deletions')
    parser.add_argument('--performance_by_depth', action='store_true', help='Whether to compute performance by mutation depth')
    args = parser.parse_args()
    
    mapping_protein_seq_DMS = pd.read_csv(args.DMS_reference_file_path)
    num_DMS=len(mapping_protein_seq_DMS)
    print("There are {} DMSs in mapping file".format(num_DMS))

    if not args.indel_mode:
        uniprot_Neff_lookup = mapping_protein_seq_DMS[['UniProt_ID','MSA_Neff_L_category']].drop_duplicates()
        uniprot_Neff_lookup.columns=['UniProt_ID','Neff_L_category']
        uniprot_taxon_lookup = mapping_protein_seq_DMS[['UniProt_ID','taxon']].drop_duplicates()
        uniprot_taxon_lookup.columns=['UniProt_ID','Taxon']
    else:
        args.performance_by_depth = False

    if args.model_list=="tranception_only":
        score_variables = ['Tranception']
    elif args.model_list=="all_models":
        score_file = pd.read_csv(args.input_scoring_files_folder+os.sep+mapping_protein_seq_DMS["DMS_filename"].values[0])
        score_variables = [ x for x in score_file.columns if x not in ['DMS_score','DMS_score_bin','mutant','mutated_sequence']]
    
    if not os.path.isdir(args.output_performance_file_folder):
        os.mkdir(args.output_performance_file_folder)
        for metric in ['Spearman','AUC','MCC']:
            os.mkdir(args.output_performance_file_folder+os.sep+metric)
    
    model_types={}
    alignment_based_models = ['EVE','DeepSequence','Wavenet','Site_Independent','EVmutation']
    for score in score_variables:
        if any(x in score for x in alignment_based_models):
            model_types[score]='Alignment-based model'
        else:
            model_types[score]='Protein language model'
        if "MSA_Transformer" in score:
            model_types[score]='Hybrid model'
        if score in ['Tranception_S_retrieval','Tranception_M_retrieval','Tranception_L_retrieval','Ensemble_Tranception_EVE']:
            model_types[score]='Hybrid model'
        
        
    model_types=pd.DataFrame.from_dict(model_types,columns=['Model type'],orient='index')

    model_details={
        'Tranception_L_no_retrieval':'Tranception Large model (700M params) without retrieval',
        'Tranception_S_retrieval':'Tranception Small model (85M params) with retrieval',
        'Tranception_M_retrieval':'Tranception Medium model (300M params) with retrieval',
        'Tranception_L_retrieval':'Tranception Large model (700M params) with retrieval',
        'EVE_single':'EVE model (single seed)',
        'EVE_ensemble':'EVE model (ensemble of 5 independently-trained models)',
        'MSA_Transformer_single':'MSA Transformer (single MSA sample)',
        'MSA_Transformer_ensemble':'MSA Transformer (ensemble of 5 MSA samples)',
        'ESM1v_single':'ESM-1v (single seed)',
        'ESM1v_ensemble':'ESM-1v (ensemble of 5 independently-trained models)',
        'Wavenet':'Wavenet model',
        'DeepSequence_single':'DeepSequence model (single seed)',
        'DeepSequence_ensemble':'DeepSequence model (ensemble of 5 independently-trained models)',
        'Site_Independent':'Site-Independent model',
        'EVmutation':'EVmutation model',
        'RITA_s':'RITA small model (85M params)',
        'RITA_m':'RITA medium model (300M params)',
        'RITA_l':'RITA large model (680M params)',
        'RITA_xl':'RITA xlarge model (1.2B params)',
        'RITA_ensemble':'Ensemble of the 4 RITA models',
        'Progen2_small':'Progen2 small model (150M params)',
        'Progen2_medium':'Progen2 medium model (760M params)',
        'Progen2_base':'Progen2 base model (760M params)',
        'Progen2_large':'Progen2 large model (2.7B params)',
        'Progen2_xlarge':'Progen2 xlarge model (6.4B params)',
        'Progen2_ensemble':'Ensemble of the 5 Progen2 models',
        'Ensemble_Tranception_EVE':'Ensemble of Tranception Large with retrieval and EVE (ensemble)',
    }
    model_details=pd.DataFrame.from_dict(model_details,columns=['Model details'],orient='index')

    model_references={
        'Tranception_L_no_retrieval':"<a href='https://proceedings.mlr.press/v162/notin22a.html'>Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A.N., Marks, D.S., & Gal, Y. (2022). Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval. ICML.</a>",
        'Tranception_S_retrieval':"<a href='https://proceedings.mlr.press/v162/notin22a.html'>Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A.N., Marks, D.S., & Gal, Y. (2022). Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval. ICML.</a>",
        'Tranception_M_retrieval':"<a href='https://proceedings.mlr.press/v162/notin22a.html'>Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A.N., Marks, D.S., & Gal, Y. (2022). Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval. ICML.</a>",
        'Tranception_L_retrieval':"<a href='https://proceedings.mlr.press/v162/notin22a.html'>Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A.N., Marks, D.S., & Gal, Y. (2022). Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval. ICML.</a>",
        'EVE_single':"<a href='https://www.nature.com/articles/s41586-021-04043-8'>Frazer, J., Notin, P., Dias, M., Gomez, A.N., Min, J.K., Brock, K.P., Gal, Y., & Marks, D.S. (2021). Disease variant prediction with deep generative models of evolutionary data. Nature.</a>",
        'EVE_ensemble':"<a href='https://www.nature.com/articles/s41586-021-04043-8'>Frazer, J., Notin, P., Dias, M., Gomez, A.N., Min, J.K., Brock, K.P., Gal, Y., & Marks, D.S. (2021). Disease variant prediction with deep generative models of evolutionary data. Nature.</a>",
        'MSA_Transformer_single':"<a href='http://proceedings.mlr.press/v139/rao21a.html'>Rao, R., Liu, J., Verkuil, R., Meier, J., Canny, J.F., Abbeel, P., Sercu, T., & Rives, A. (2021). MSA Transformer. ICML.</a>",
        'MSA_Transformer_ensemble':"<a href='http://proceedings.mlr.press/v139/rao21a.html'>Rao, R., Liu, J., Verkuil, R., Meier, J., Canny, J.F., Abbeel, P., Sercu, T., & Rives, A. (2021). MSA Transformer. ICML.</a>",
        'ESM1v_single':"<a href='https://proceedings.neurips.cc/paper/2021/hash/f51338d736f95dd42427296047067694-Abstract.html'>Meier, J., Rao, R., Verkuil, R., Liu, J., Sercu, T., & Rives, A. (2021). Language models enable zero-shot prediction of the effects of mutations on protein function. NeurIPS.</a>",
        'ESM1v_ensemble':"<a href='https://proceedings.neurips.cc/paper/2021/hash/f51338d736f95dd42427296047067694-Abstract.html'>Meier, J., Rao, R., Verkuil, R., Liu, J., Sercu, T., & Rives, A. (2021). Language models enable zero-shot prediction of the effects of mutations on protein function. NeurIPS.</a>",
        'Wavenet':"<a href='https://www.nature.com/articles/s41467-021-22732-w'>Shin, J., Riesselman, A.J., Kollasch, A.W., McMahon, C., Simon, E., Sander, C., Manglik, A., Kruse, A.C., & Marks, D.S. (2021). Protein design and variant prediction using autoregressive generative models. Nature Communications, 12.</a>",
        'DeepSequence_single':"<a href='https://www.nature.com/articles/s41592-018-0138-4'>Riesselman, A.J., Ingraham, J., & Marks, D.S. (2018). Deep generative models of genetic variation capture the effects of mutations. Nature Methods, 15, 816-822.</a>",
        'DeepSequence_ensemble':"<a href='https://www.nature.com/articles/s41592-018-0138-4'>Riesselman, A.J., Ingraham, J., & Marks, D.S. (2018). Deep generative models of genetic variation capture the effects of mutations. Nature Methods, 15, 816-822.</a>",
        'Site_Independent':"<a href='https://www.nature.com/articles/nbt.3769'>Hopf, T.A., Ingraham, J., Poelwijk, F.J., Schärfe, C.P., Springer, M., Sander, C., & Marks, D.S. (2017). Mutation effects predicted from sequence co-variation. Nature Biotechnology, 35, 128-135.</a>",
        'EVmutation':"<a href='https://www.nature.com/articles/nbt.3769'>Hopf, T.A., Ingraham, J., Poelwijk, F.J., Schärfe, C.P., Springer, M., Sander, C., & Marks, D.S. (2017). Mutation effects predicted from sequence co-variation. Nature Biotechnology, 35, 128-135.</a>",
        'RITA_s':"<a href='https://arxiv.org/abs/2205.05789'>Hesslow, D., Zanichelli, N., Notin, P., Poli, I., & Marks, D.S. (2022). RITA: a Study on Scaling Up Generative Protein Sequence Models. ArXiv, abs/2205.05789.</a>",
        'RITA_m':"<a href='https://arxiv.org/abs/2205.05789'>Hesslow, D., Zanichelli, N., Notin, P., Poli, I., & Marks, D.S. (2022). RITA: a Study on Scaling Up Generative Protein Sequence Models. ArXiv, abs/2205.05789.</a>",
        'RITA_l':"<a href='https://arxiv.org/abs/2205.05789'>Hesslow, D., Zanichelli, N., Notin, P., Poli, I., & Marks, D.S. (2022). RITA: a Study on Scaling Up Generative Protein Sequence Models. ArXiv, abs/2205.05789.</a>",
        'RITA_xl':"<a href='https://arxiv.org/abs/2205.05789'>Hesslow, D., Zanichelli, N., Notin, P., Poli, I., & Marks, D.S. (2022). RITA: a Study on Scaling Up Generative Protein Sequence Models. ArXiv, abs/2205.05789.</a>",
        'RITA_ensemble':"<a href='https://arxiv.org/abs/2205.05789'>Hesslow, D., Zanichelli, N., Notin, P., Poli, I., & Marks, D.S. (2022). RITA: a Study on Scaling Up Generative Protein Sequence Models. ArXiv, abs/2205.05789.</a>",
        'Progen2_small':"<a href='https://arxiv.org/abs/2206.13517'> Nijkamp, E., Ruffolo, J.A., Weinstein, E.N., Naik, N., & Madani, A. (2022). ProGen2: Exploring the Boundaries of Protein Language Models. ArXiv, abs/2206.13517. </a>",
        'Progen2_medium':"<a href='https://arxiv.org/abs/2206.13517'> Nijkamp, E., Ruffolo, J.A., Weinstein, E.N., Naik, N., & Madani, A. (2022). ProGen2: Exploring the Boundaries of Protein Language Models. ArXiv, abs/2206.13517. </a>",
        'Progen2_base':"<a href='https://arxiv.org/abs/2206.13517'> Nijkamp, E., Ruffolo, J.A., Weinstein, E.N., Naik, N., & Madani, A. (2022). ProGen2: Exploring the Boundaries of Protein Language Models. ArXiv, abs/2206.13517. </a>",
        'Progen2_large':"<a href='https://arxiv.org/abs/2206.13517'> Nijkamp, E., Ruffolo, J.A., Weinstein, E.N., Naik, N., & Madani, A. (2022). ProGen2: Exploring the Boundaries of Protein Language Models. ArXiv, abs/2206.13517. </a>",
        'Progen2_xlarge':"<a href='https://arxiv.org/abs/2206.13517'> Nijkamp, E., Ruffolo, J.A., Weinstein, E.N., Naik, N., & Madani, A. (2022). ProGen2: Exploring the Boundaries of Protein Language Models. ArXiv, abs/2206.13517. </a>",
        'Progen2_ensemble':"<a href='https://arxiv.org/abs/2206.13517'> Nijkamp, E., Ruffolo, J.A., Weinstein, E.N., Naik, N., & Madani, A. (2022). ProGen2: Exploring the Boundaries of Protein Language Models. ArXiv, abs/2206.13517. </a>",
        'Ensemble_Tranception_EVE':"<a href='https://proceedings.mlr.press/v162/notin22a.html'>Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A.N., Marks, D.S., & Gal, Y. (2022). Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval. ICML.</a>",
    }
    model_references=pd.DataFrame.from_dict(model_references,columns=['References'],orient='index')

    clean_names={
        'Tranception_L_no_retrieval':'Tranception L no retrieval',
        'Tranception_S_retrieval':'Tranception S',
        'Tranception_M_retrieval':'Tranception M',
        'Tranception_L_retrieval':'Tranception L',
        'EVE_single':'EVE (single)',
        'EVE_ensemble':'EVE (ensemble)',
        'MSA_Transformer_single':'MSA Transformer (single)',
        'MSA_Transformer_ensemble':'MSA Transformer (ensemble)',
        'ESM1v_single':'ESM-1v (single)',
        'ESM1v_ensemble':'ESM-1v (ensemble)',
        'Wavenet':'Wavenet',
        'DeepSequence_single':'DeepSequence (single)',
        'DeepSequence_ensemble':'DeepSequence (ensemble)',
        'Site_Independent':'Site-Independent',
        'EVmutation':'EVmutation',
        'RITA_s':'RITA S',
        'RITA_m':'RITA M',
        'RITA_l':'RITA L',
        'RITA_xl':'RITA XL',
        'RITA_ensemble':'RITA (ensemble)',
        'Progen2_small':'Progen2 S',
        'Progen2_medium':'Progen2 M',
        'Progen2_base':'Progen2 Base',
        'Progen2_large':'Progen2 L',
        'Progen2_xlarge':'Progen2 XL',
        'Progen2_ensemble':'Progen2 (ensemble)',
        'Ensemble_Tranception_EVE':'Ensemble Tranception & EVE',
    }

    performance_all_DMS={}
    output_filename={}
    for metric in ['Spearman','AUC','MCC']:
        performance_all_DMS[metric]={}
        mutation_type = "substitutions" if not args.indel_mode else "indels"
        output_filename[metric]=args.model_list+"_"+mutation_type+"_"+metric
        for i, score in enumerate(score_variables):
            performance_all_DMS[metric][score]=i
            if not args.indel_mode and args.performance_by_depth:
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
                tranception = tranception[['mutated_sequence','avg_score']]
                tranception.columns=['mutated_sequence','Tranception']    
                merged_scores = pd.merge(DMS_file, tranception, on='mutated_sequence', how='inner')
                merged_scores.dropna(inplace=True)
            elif args.model_list=="all_models":
                merged_scores = pd.read_csv(args.input_scoring_files_folder + os.sep + DMS_id + ".csv") #We assume no missing value (all models were enforced to score all mutants)
            if 'mutant' not in merged_scores: merged_scores['mutant'] = merged_scores['mutated_sequence'] #if mutant not in DMS file we default to mutated_sequence (eg., for indels)
        except:
            print("At least one scoring file missing")
            continue

        if not args.indel_mode and args.performance_by_depth:
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
                print("AUC issue with: {} for model: {}".format(DMS_id,score))
                performance_DMS['AUC'][score] = np.nan
            try:
                median_cutoff=merged_scores[score].median()
                merged_scores[score+"_bin"]=merged_scores[score].map(lambda x: 1 if x >= median_cutoff else 0)
                performance_DMS['MCC'][score] = matthews_corrcoef(y_true=merged_scores['DMS_score_bin'], y_pred=merged_scores[score+"_bin"])
            except:
                print("MCC issue with: {} for model: {}".format(DMS_id,score))
                performance_DMS['MCC'][score] = np.nan
        
        if not args.indel_mode and args.performance_by_depth:
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
        if args.indel_mode: 
            bootstrap_standard_error = pd.DataFrame(compute_bootstrap_standard_error(performance_all_DMS[metric].subtract(performance_all_DMS[metric]['Tranception_M_retrieval'],axis=0)),columns=["Bootstrap_standard_error_"+metric])
            performance_all_DMS[metric].loc['Average'] = performance_all_DMS[metric].mean() #DMS-level average = Uniprot-level average for indels
        for var in performance_all_DMS[metric]:
            if var not in ['UniProt_ID','Neff_L_category','Taxon']:
                performance_all_DMS[metric][var]=performance_all_DMS[metric][var].astype(float).round(3)
            if var in ['number_mutants']:
                performance_all_DMS[metric][var]=performance_all_DMS[metric][var].astype(int)
        performance_all_DMS[metric].to_csv(args.output_performance_file_folder + os.sep + metric + os.sep + output_filename[metric] + '_DMS_level.csv')
        if not args.indel_mode and args.performance_by_depth:
            all_columns = performance_all_DMS[metric].columns
            performance_all_DMS_html=performance_all_DMS[metric].copy()
            performance_all_DMS_html.columns=performance_all_DMS_html.columns.map(lambda x: clean_names[x] if x in clean_names else x)
            all_not_depth_columns = all_columns[[all_columns[x].split("_")[-1] not in ['1','2','3','4','5+'] for x in range(len(all_columns))]]
            all_not_depth_columns_clean = all_not_depth_columns.map(lambda x: clean_names[x] if x in clean_names else x)
            performance_all_DMS_html[all_not_depth_columns_clean].to_html(args.output_performance_file_folder + os.sep + metric + os.sep + output_filename[metric] + '_DMS_level.html')
        else:
            performance_all_DMS_html=performance_all_DMS[metric].copy()
            performance_all_DMS_html.columns = performance_all_DMS_html.columns.map(lambda x: clean_names[x] if x in clean_names else x)
            performance_all_DMS_html.to_html(args.output_performance_file_folder + os.sep + metric + os.sep + output_filename[metric] + '_DMS_level.html')
        
        if not args.indel_mode:
            uniprot_metric_performance = performance_all_DMS[metric].groupby(['UniProt_ID']).mean()#.reset_index()
            bootstrap_standard_error = pd.DataFrame(compute_bootstrap_standard_error(uniprot_metric_performance.subtract(uniprot_metric_performance['Ensemble_Tranception_EVE'],axis=0)),columns=["Bootstrap_standard_error_"+metric])
            uniprot_metric_performance = uniprot_metric_performance.reset_index()
            uniprot_metric_performance = pd.merge(uniprot_metric_performance,uniprot_Neff_lookup,on='UniProt_ID', how='left')
            uniprot_metric_performance = pd.merge(uniprot_metric_performance,uniprot_taxon_lookup,on='UniProt_ID', how='left')
            del uniprot_metric_performance['number_mutants']
            uniprot_level_average = uniprot_metric_performance.mean()
            uniprot_metric_performance.loc['Average'] = uniprot_level_average
            uniprot_metric_performance=uniprot_metric_performance.round(3)
            uniprot_metric_performance.to_csv(args.output_performance_file_folder + os.sep + metric + os.sep + output_filename[metric] + '_Uniprot_level.csv', index=False)
            
            if args.performance_by_depth:
                performance_by_depth = {}
                all_not_depth_columns = [x for x in all_not_depth_columns if x != 'number_mutants']
                for depth in ['1','2','3','4','5+']:
                    depth_columns = all_columns[[all_columns[x].split("_")[-1]==depth for x in range(len(all_columns))]]
                    performance_by_depth[depth] = uniprot_metric_performance.loc['Average',depth_columns].transpose().reset_index()
                    performance_by_depth[depth]['model_name'] = performance_by_depth[depth]['score'].map(lambda x: '_'.join(x.split('_')[:-1]))
                    performance_by_depth[depth]=performance_by_depth[depth][['model_name','Average']]
                    performance_by_depth[depth].columns = ['model_name','Depth_'+depth]
                    performance_by_depth[depth].set_index('model_name', inplace=True)
                uniprot_metric_performance = uniprot_metric_performance[all_not_depth_columns]
            performance_by_MSA_depth = uniprot_metric_performance.groupby(['Neff_L_category']).mean().transpose()
            performance_by_MSA_depth = performance_by_MSA_depth[['low','medium','high']]
            performance_by_MSA_depth.columns = ['Low_MSA_depth','Medium_MSA_depth','High_MSA_depth']
            performance_by_taxon = uniprot_metric_performance.groupby(['Taxon']).mean().transpose()
            performance_by_taxon = performance_by_taxon[['Human','Eukaryote','Prokaryote','Virus']]
            performance_by_taxon.columns = ['Human','Other Eukaryote','Prokaryote','Virus']
            
            summary_performance = pd.merge(pd.DataFrame(uniprot_level_average,columns=['Average_'+metric]), performance_by_MSA_depth,left_index=True, right_index=True,how='inner')
            summary_performance = pd.merge(summary_performance, performance_by_taxon,left_index=True, right_index=True,how='inner')
            if args.performance_by_depth:
                for depth in ['1','2','3','4','5+']:
                    summary_performance = pd.merge(summary_performance, performance_by_depth[depth],left_index=True, right_index=True,how='inner')
            final_column_order = ['Model_name','Model type','Average_'+metric,'Bootstrap_standard_error_'+metric,'Low_MSA_depth','Medium_MSA_depth','High_MSA_depth','Human','Other Eukaryote','Prokaryote','Virus','Depth_1','Depth_2','Depth_3','Depth_4','Depth_5+','Model details','References']
        else:
            del performance_all_DMS[metric]['number_mutants']
            summary_performance = pd.DataFrame(performance_all_DMS[metric].transpose().loc[:,'Average'])
            summary_performance.columns = ['Average_'+metric]
            final_column_order = ['Model_name','Model type','Average_'+metric,'Bootstrap_standard_error_'+metric,'Model details','References']
        
        summary_performance.sort_values(by='Average_'+metric,ascending=False,inplace=True)
        summary_performance.index.name = 'Model_name'
        summary_performance.reset_index(inplace=True)
        summary_performance.index = range(1,len(summary_performance)+1)
        summary_performance.index.name = 'Model_rank'
        summary_performance = pd.merge(summary_performance, bootstrap_standard_error, left_on='Model_name', right_index=True, how='left')
        summary_performance = pd.merge(summary_performance, model_types, left_on='Model_name', right_index=True, how='left')
        summary_performance = pd.merge(summary_performance, model_details, left_on='Model_name', right_index=True, how='left')
        summary_performance = pd.merge(summary_performance, model_references, left_on='Model_name', right_index=True, how='left')
        summary_performance=summary_performance.round(3)
        summary_performance['Model_name']=summary_performance['Model_name'].map(lambda x: clean_names[x] if x in clean_names else x)
        
        summary_performance=summary_performance.reindex(columns=final_column_order)
        summary_performance.to_csv(args.output_performance_file_folder + os.sep + metric + os.sep + 'Summary_performance_'+output_filename[metric]+'.csv')
        summary_performance.to_html(args.output_performance_file_folder + os.sep + metric + os.sep + 'Summary_performance_'+output_filename[metric]+'.html',escape=False)

if __name__ == '__main__':
    main()