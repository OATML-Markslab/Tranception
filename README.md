# Tranception

This is the official code repository for the paper "Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval". This project is a joint collaboration between the [Marks lab](https://www.deboramarkslab.com/) and the [OATML group](https://oatml.cs.ox.ac.uk/).

## Abstract
The ability to accurately model the fitness landscape of protein sequences is critical to a wide range of applications, from quantifying the effects of human variants on disease likelihood, to predicting immune-escape mutations in viruses and designing novel biotherapeutic proteins. Deep generative models of protein sequences trained on multiple sequence alignments have been the most successful approaches so far to address these tasks. The performance of these methods is however contingent on the availability of sufficiently deep and diverse alignments for reliable training. Their potential scope is thus limited by the fact many protein families are hard, if not impossible, to align. Large language models trained on massive quantities of non-aligned protein sequences from diverse families address these problems and show potential to eventually bridge the performance gap. We introduce Tranception, a novel transformer architecture leveraging autoregressive predictions and retrieval of homologous sequences at inference to achieve state-of-the-art fitness prediction performance. Given its markedly higher performance on multiple mutants, robustness to shallow alignments and ability to score indels, our approach offers significant gain of scope over existing approaches. To enable more rigorous model testing across a broader range of protein families, we develop ProteinGym -- an extensive set of multiplexed assays of variant effects, substantially increasing both the number and diversity of assays compared to existing benchmarks.

## Setup
You may download the Tranception repository and create a conda environment with the proper dependencies (as listed in `tranception_env.yml`) as follows:
```
git clone https://github.com/OATML-Markslab/Tranception.git
conda env create -f tranception_env.yml
```
To leverage retrieval when scoring insertions and deletions, you will need [Clustal Omega](http://www.clustal.org/omega/#Download).

## Tranception
Tranception is a novel autoregressive transformer architecture that was designed with two core principles in mind: 1) promoting specialization across attention heads 2) explicitly extracting patterns from contiguous subsequences.

To download the *Tranception Large* model checkpoint (~3.1GB unzipped):
```
curl -o Tranception_Large_checkpoint.zip https://marks.hms.harvard.edu/tranception/Tranception_Large_checkpoint.zip
unzip Tranception_Large_checkpoint.zip
rm Tranception_Large_checkpoint.zip
```

The Tranception model checkpoint is also made available through the [Huggging Face hub](https://huggingface.co/OATML-Markslab/Tranception).

When scoring with retrieval, we compute weighted pseudocounts at each position using sequence weights as per the procedure described in [Hopf et al.](https://www.nature.com/articles/nbt.3769).
Weights for all proteins in the ProteinGym benchmarks may be downloaded as follows (~68M unzipped):
```
curl -o MSA_weights.zip https://marks.hms.harvard.edu/tranception/MSA_weights.zip
unzip MSA_weights.zip
rm MSA_weights.zip
```
To compute sequence weights for new proteins, you may use the MSA_processing class under `tranception/utils/msa_utils.py`.

The `examples` folder provides several bash scripts that may be used for scoring and evaluating Tranception on the ProteinGym benchmarks.

## ProteinGym
ProteinGym is an extensive set of Deep Mutational Scanning (DMS) assays curated to enable thorough comparisons of various mutation effect predictors indifferent regimes. It is comprised of two benchmarks: 1) a substitution benchmark which consists of the experimental characterisation of ∼1.5M missense variants across 87 DMS assays 2) an indel benchmark that includes ∼300k mutants across 7 DMS assays.

Each processed file in each benchmark corresponds to a single DMS assay, and contains the following variables:
- mutant (str): describes the set of substitutions to apply on the reference sequence to obtain the mutated sequence (eg., A1P:D2N implies the amino acid 'A' at position 1 should be replaced by 'P', and 'D' at position 2 should be replaced by 'N'). Present in the the ProteinGym substitution assays only (not indels).
- mutated_sequence (str): represents the full amino acid sequence for the mutated protein.
- DMS_score (float): corresponds to the experimental measurement in the DMS assay. Across all assays, the higher the DMS_score value, the higher the fitness of the mutated protein
- DMS_score_bin (int): indicates whether the DMS_score is above the fitness cutoff (1 is fit, 0 is not fit)

Additionally, we provide reference files in the [ProteinGym folder](https://github.com/OATML-Markslab/Tranception/tree/main/proteingym) that give further details on each assay and contain in particular:
- The UniProt_ID of the corresponding protein, along with taxon and MSA depth category
- The target sequence (target_seq) used in the assay
- Details on how the DMS_score was created from the raw files and how it was binarized 

To download the substitution benchmark (~867M unzipped):
```
curl -o ProteinGym_substitutions.zip https://marks.hms.harvard.edu/tranception/ProteinGym_substitutions.zip 
unzip ProteinGym_substitutions.zip
rm ProteinGym_substitutions.zip
```

Similarly, to download the indel benchmark (~223M unzipped):
```
curl -o ProteinGym_indels.zip https://marks.hms.harvard.edu/tranception/ProteinGym_indels.zip
unzip ProteinGym_indels.zip
rm ProteinGym_indels.zip
```

The ProteinGym benchmarks are also available on the [Hugging Face Hub](https://huggingface.co/datasets/OATML-Markslab/ProteinGym).

## Fitness prediction performance

The [proteingym folder](https://github.com/OATML-Markslab/Tranception/tree/main/ProteinGym) provides detailed performance files for Tranception and baseline models on the two ProteinGym benchmarks, for the 3 focus metrics in the paper (eg., Spearman's rank, AUC, MCC).

We recommand to aggregate fitness pediction performance at the Uniprot ID level to avoid biasing results towards proteins for which several DMS assays are available in ProteinGym. The corresponding aggregated files are suffixed with "_Uniprot_level", while the non aggregated performance files are suffixed with "_DMS_level".
Furthermore, to enable fair comparison with models trained multiple-sequence alignments (eg., EVE, DeepSequence, EVmutation), we only evaluate on the subset of mutations where position coverage is deemed high enough by these models to make a prediction. The corresponding files are preffixed with "all_models_". For comprehensiveness, we also provide performance files on all possible mutants available in ProteinGym (preffixed with "all_mutants_"), comparing only with the baselines that are able to score all mutants.

Note that for the ProteinGym indel benchmark, baselines that are able to score indels do not have the aforementionned coverage constraints (ie., no distinction between "all_models_" and "all_mutants_") and there is at most one DMS per Uniprot_ID (ie., no difference between "_Uniprot_level" and "_DMS_level"). We thus only provide one set of performance metrics for that benchmark.

### ProteinGym substitution benchmark - Leaderboard
The table below provides the average Spearman's rank correlation between DMS experimental fitness measurements and fitness predictions from Tranception or other baselines on the ProteinGym substitution benchmark. Following the terminology introduced above, we report the performance at the "Uniprot" level for "All models".

Rank | Model name | Spearman | Reference
--- | --- | --- | --- |
1 | Ensemble Tranception & EVE | 0.476 | [Notin et al.](https://arxiv.org/abs/2205.13760)
2 | Tranception (w/ retrieval) | 0.451 | [Notin et al.](https://arxiv.org/abs/2205.13760)
3 | EVE | 0.448 | [Frazer et al.](https://www.nature.com/articles/s41586-021-04043-8)
4 | EVmutation | 0.427 | [Hopf et al.](https://www.nature.com/articles/nbt.3769)
5 | MSA Transformer | 0.422 | [Rao et al.](https://proceedings.mlr.press/v139/rao21a.html)
6 | DeepSequence | 0.415 | [Riesselman et al.](https://www.nature.com/articles/s41592-018-0138-4)
7 | Tranception (no retrieval) | 0.406 | [Notin et al.](https://arxiv.org/abs/2205.13760)
8 | Wavenet | 0.398 | [Shin et al.](https://www.nature.com/articles/s41467-021-22732-w)
9 | Site Independent | 0.397 | [Hopf et al.](https://www.nature.com/articles/nbt.3769)
10 | ESM-1v | 0.371 | [Meier et al.](https://proceedings.neurips.cc/paper/2021/hash/f51338d736f95dd42427296047067694-Abstract.html)

### ProteinGym indel benchmark - Leaderboard
The table below provides the average Spearman's rank correlation between DMS experimental fitness measurements and fitness predictions from Tranception or other baselines on the ProteinGym indel benchmark.

Rank | Model name | Spearman | Reference
--- | --- | --- | --- |
1 | Tranception (w/ retrieval) | 0.463 | [Notin et al.](https://arxiv.org/abs/2205.13760)
2 | Tranception (no retrieval) | 0.43 | [Notin et al.](https://arxiv.org/abs/2205.13760)
3 | Wavenet | 0.412 | [Shin et al.](https://www.nature.com/articles/s41467-021-22732-w)

## Aggregated model scoring files
The scores for all DMS assays in the ProteinGym substitution benchmark for Tranception and other baselines (eg., EVE, Wavenet, ESM-1v, MSA Transformer) may be downloaded as follows;
```
curl -o scores_all_models_proteingym_substitutions.zip https://marks.hms.harvard.edu/tranception/scores_all_models_proteingym_substitutions.zip
unzip scores_all_models_proteingym_substitutions.zip
rm scores_all_models_proteingym_substitutions.zip
```
Similarly for the indel benchmark, all scoring files may be downloaded as follows:
```
curl -o scores_all_models_proteingym_indels.zip https://marks.hms.harvard.edu/tranception/scores_all_models_proteingym_indels.zip
unzip scores_all_models_proteingym_indels.zip
rm scores_all_models_proteingym_indels.zip
```

## Multiple Sequence Alignments (MSAs)

The MSAs used to train alignment-based methods or used at inference in Tranception with retrieval and MSA Transformer may be downloaded as follows (~2.2GB unzipped):
```
curl -o MSA_ProteinGym.zip https://marks.hms.harvard.edu/tranception/MSA_ProteinGym.zip
unzip MSA_ProteinGym.zip
rm MSA_ProteinGym.zip
```

## License
This project is available under the MIT license.

## Reference
If you use Tranception, ProteinGym or other files provided through this repository (eg., aggregated model scoring files) in your work, please cite the following paper:
```
Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A., Marks, D.S., Gal, Y. (2022). Tranception: Protein Fitness Prediction with Autoregressive Transformers and Inference-time Retrieval. ICML.
```

## Links
Pre-print: https://arxiv.org/abs/2205.13760