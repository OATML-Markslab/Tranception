# Tranception

This is the official code repository for the paper "Tranception: Protein Fitness Prediction with Autoregressive Transformers and Inference-time Retrieval". This project is a joint collaboration between the [Marks lab](https://www.deboramarkslab.com/) and the [OATML group](https://oatml.cs.ox.ac.uk/).

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

When scoring with retrieval, we compute weighted pseudocounts at each position using sequence weights as per the exact procedure described in [Hopf et al.](https://www.nature.com/articles/nbt.3769), using the default hyperparameter values (eg., threshold_sequence_frac_gaps = 0.5 and threshold_focus_cols_frac_gaps = 0.3).
Weights for all proteins in the ProteinGym benchmarks may be downloaded as follows (~68M unzipped):
```
curl -o MSA_weights.zip https://marks.hms.harvard.edu/tranception/MSA_weights.zip
unzip MSA_weights.zip
rm MSA_weights.zip
```
To compute sequence weights for new proteins, you may use the MSA_processing class under `tranception/utils/msa_utils.py`.

We also provide weights computed with the procedure from [Hopf et al.](https://www.nature.com/articles/nbt.3769) but without performing any position-specific filtering based on minimum coverage (ie., threshold_focus_cols_frac_gaps=1.0). This is used for instance to enable [EVE models](https://www.nature.com/articles/s41586-021-04043-8) to score mutants at all positions in the input protein sequence (see more details in this [issue](https://github.com/OATML-Markslab/Tranception/issues/9)). These weights may be downloaded as follows:
```
curl -o MSA_weights_substitutions_all_positions.zip https://marks.hms.harvard.edu/tranception/MSA_weights_substitutions_all_positions.zip
curl -o MSA_weights_indels_all_positions.zip https://marks.hms.harvard.edu/tranception/MSA_weights_indels_all_positions.zip
```

### Fitness prediction with Tranception & performance evaluation on ProteinGym
To predict the fitness of mutated sequences (substitutions and indels) with Tranception, you may use the [main scoring script](https://github.com/OATML-Markslab/Tranception/blob/main/score_tranception_proteingym.py). To then evaluate the performance of these predictions against ground truth labels (eg., DMS assays from the ProteinGym benchmark) across the different metrics discussed in the paper (ie., Spearman's rank correlation, AUC, MCC) you may use the [performance analysis script](https://github.com/OATML-Markslab/Tranception/blob/main/performance_analysis_proteingym.py).
The [examples](https://github.com/OATML-Markslab/Tranception/tree/main/examples) folder provides several bash scripts illustrating how these two scripts may be used.

## ProteinGym benchmarks

ProteinGym is an extensive set of Deep Mutational Scanning (DMS) assays curated to enable thorough comparisons of various mutation effect predictors indifferent regimes. It is comprised of two benchmarks: 1) a substitution benchmark which consists of the experimental characterisation of ∼1.5M missense variants across 87 DMS assays 2) an indel benchmark that includes ∼300k mutants across 7 DMS assays.
For more details, please refer to the ProteinGym [repo](https://github.com/OATML-Markslab/ProteinGym) and [website](https://www.proteingym.org/).

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

Note: more recent versions of the scoring files (eg., with more baselines) are available on the ProteinGym [repo](https://github.com/OATML-Markslab/ProteinGym).

## Multiple Sequence Alignments (MSAs)

The MSAs used to train alignment-based methods or used at inference in Tranception with retrieval and MSA Transformer may be downloaded as follows (~2.2GB unzipped):
```
curl -o MSA_ProteinGym.zip https://marks.hms.harvard.edu/tranception/MSA_ProteinGym.zip
unzip MSA_ProteinGym.zip
rm MSA_ProteinGym.zip
```

## Protein Design

As a powerful generative model with SOTA fitness prediction capabilities, Tranception is well-suited to support the design of new proteins. To illustrate these capabilities, we built a Gradio app that enables in silico directed evolution for iterative protein redesign.
The app is available on [Hugging Face spaces](https://huggingface.co/spaces/PascalNotin/Tranception_design) (ideal for short proteins/mutation ranges) and as a [Colab notebook](https://colab.research.google.com/drive/12ni4U1Na9VWwnwdpGzYHCkNC7sud9SLy?usp=sharing) (w/ GPU support each directed evolution cycle takes ~5 mins for full proteins).

## License
This project is available under the MIT license.

## Reference
If you use Tranception, ProteinGym or other files provided through this repository (eg., aggregated model scoring files) in your work, please cite the following paper:
```
Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A., Marks, D.S., Gal, Y. (2022). Tranception: Protein Fitness Prediction with Autoregressive Transformers and Inference-time Retrieval. ICML.
```

## Links
- ICML proceedings: https://proceedings.mlr.press/v162/notin22a.html
- Arxiv pre-print: https://arxiv.org/abs/2205.13760
- Hugging Face Hub (model checkpoints): https://huggingface.co/OATML-Markslab/Tranception
- ProteinGym website: https://www.proteingym.org/
- ProteinGym repo: https://github.com/OATML-Markslab/ProteinGym
- Design app: [Hugging Face spaces](https://huggingface.co/spaces/PascalNotin/Tranception_design) or [Colab](https://colab.research.google.com/drive/12ni4U1Na9VWwnwdpGzYHCkNC7sud9SLy?usp=sharing)
