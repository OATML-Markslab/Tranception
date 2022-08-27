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

To download the *Tranception Medium* model checkpoint (~1.4GB unzipped):
```
curl -o Tranception_Medium_checkpoint.zip https://marks.hms.harvard.edu/tranception/Tranception_Medium_checkpoint.zip
unzip Tranception_Medium_checkpoint.zip
rm Tranception_Medium_checkpoint.zip
```

To download the *Tranception Small* model checkpoint (~400MB unzipped):
```
curl -o Tranception_Small_checkpoint.zip https://marks.hms.harvard.edu/tranception/Tranception_Small_checkpoint.zip
unzip Tranception_Small_checkpoint.zip
rm Tranception_Small_checkpoint.zip
```

The Tranception model checkpoints are also made available through the [Huggging Face hub](https://huggingface.co/OATML-Markslab/Tranception_Large).

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

We recommand to aggregate fitness pediction performance at the Uniprot ID level to avoid biasing results towards proteins for which several DMS assays are available in ProteinGym. The corresponding aggregated files are suffixed with "_Uniprot_level", while the non aggregated performance files are suffixed with "_DMS_level". Note that there is at most one DMS per Uniprot_ID for the ProteinGym indel benchmark (ie., no difference between "_Uniprot_level" and "_DMS_level"). We thus only provide one set of performance metrics for that benchmark.

Note: an earlier version of our benchmark contained scores for all model baselines on the subset of mutations happening at well-covered positions in the corresponding MSAs. We since extended the scoring of alignment-based methods (eg., EVE, DeepSequence, EVmutation) to all sequence positions and noticed minimal changes to the resulting relative performance and score rankings. As a result, we are now only reporting performance on the full set of all DMS mutations for simplicity.

### ProteinGym benchmarks - Leaderboard

The full ProteinGym benchmarks performance files are also accessible via our dedicated website: https://www.proteingym.org/.
It includes leaderboards for the substitution and indel benchmarks, as well as detailed DMS-level performance files for all baselines.
The current version of the substitution benchmark includes the following baselines:

Model name | Model type | Reference
--- | --- | --- |
Site Independent | Alignment-based model | [Hopf, T.A., Ingraham, J., Poelwijk, F.J., Schärfe, C.P., Springer, M., Sander, C., & Marks, D.S. (2017). Mutation effects predicted from sequence co-variation. Nature Biotechnology, 35, 128-135.](https://www.nature.com/articles/nbt.3769)
EVmutation | Alignment-based model | [Hopf, T.A., Ingraham, J., Poelwijk, F.J., Schärfe, C.P., Springer, M., Sander, C., & Marks, D.S. (2017). Mutation effects predicted from sequence co-variation. Nature Biotechnology, 35, 128-135.](https://www.nature.com/articles/nbt.3769)
Wavenet | Alignment-based model | [Shin, J., Riesselman, A.J., Kollasch, A.W., McMahon, C., Simon, E., Sander, C., Manglik, A., Kruse, A.C., & Marks, D.S. (2021). Protein design and variant prediction using autoregressive generative models. Nature Communications, 12.](https://www.nature.com/articles/s41467-021-22732-w)
DeepSequence | Alignment-based model | [Riesselman, A.J., Ingraham, J., & Marks, D.S. (2018). Deep generative models of genetic variation capture the effects of mutations. Nature Methods, 15, 816-822.](https://www.nature.com/articles/s41592-018-0138-4)
EVE | Alignment-based model | [Frazer, J., Notin, P., Dias, M., Gomez, A.N., Min, J.K., Brock, K.P., Gal, Y., & Marks, D.S. (2021). Disease variant prediction with deep generative models of evolutionary data. Nature.](https://www.nature.com/articles/s41586-021-04043-8)
ESM-1v | Protein language model |[Meier, J., Rao, R., Verkuil, R., Liu, J., Sercu, T., & Rives, A. (2021). Language models enable zero-shot prediction of the effects of mutations on protein function. NeurIPS.](https://proceedings.neurips.cc/paper/2021/hash/f51338d736f95dd42427296047067694-Abstract.html)
MSA Transformer | Protein language model |[Rao, R., Liu, J., Verkuil, R., Meier, J., Canny, J.F., Abbeel, P., Sercu, T., & Rives, A. (2021). MSA Transformer. ICML.](http://proceedings.mlr.press/v139/rao21a.html)
RITA | Protein language model | [Hesslow, D., Zanichelli, N., Notin, P., Poli, I., & Marks, D.S. (2022). RITA: a Study on Scaling Up Generative Protein Sequence Models. ArXiv, abs/2205.05789.](https://arxiv.org/abs/2205.05789)
Progen2 | Protein language model | [Nijkamp, E., Ruffolo, J.A., Weinstein, E.N., Naik, N., & Madani, A. (2022). ProGen2: Exploring the Boundaries of Protein Language Models. ArXiv, abs/2206.13517.](https://arxiv.org/abs/2206.13517)
Tranception | Hybrid | [Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A.N., Marks, D.S., & Gal, Y. (2022). Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval. ICML.](https://proceedings.mlr.press/v162/notin22a.html)

Except for the Wavenet model (which only uses alignments to recover a set of homologous protein sequences to train on, but then trains on non-aligned sequences), all alignment-based methods are unable to score indels given the fixed coordinate system they are trained on. Similarly, the masking procedure to generate the masked-marginals for ESM-1v and MSA Transformer requires the position to exist in the wild-type sequence. All the other model architectures listed above (eg., Tranception, RITA, Progen2) are included in the indel benchmark.

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

## Protein Design

As a powerful generative model with SOTA fitness prediction capabilities, Tranception is well-suited to support the design of new proteins. To illustrate these capabilities, we built a Gradio app that enables in silico directed evolution for iterative protein redesign.
The app is available on [Hugging Face spaces](https://huggingface.co/spaces/PascalNotin/Tranception_design) (ideal for short proteins/mutation ranges) and as a [colab notebook](https://colab.research.google.com/drive/12ni4U1Na9VWwnwdpGzYHCkNC7sud9SLy?usp=sharing) (w/ GPU support each directed evolution cycle takes ~5 mins for full proteins).

## License
This project is available under the MIT license.

## Reference
If you use Tranception, ProteinGym or other files provided through this repository (eg., aggregated model scoring files) in your work, please cite the following paper:
```
Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A., Marks, D.S., Gal, Y. (2022). Tranception: Protein Fitness Prediction with Autoregressive Transformers and Inference-time Retrieval. ICML.
```

## Links
ICML proceedings: https://proceedings.mlr.press/v162/notin22a.html
Arxiv pre-print: https://arxiv.org/abs/2205.13760
HuggingFace Hub (model checkpoints): https://huggingface.co/OATML-Markslab/Tranception
ProteinGym website: https://www.proteingym.org/