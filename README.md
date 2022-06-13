# Tranception

This is the official code repository for the paper "Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval". This project is a joint collaboration between the [Marks lab](https://www.deboramarkslab.com/) and the [OATML group](https://oatml.cs.ox.ac.uk/).

## Abstract
The ability to accurately model the fitness landscape of protein sequences is critical to a wide range of applications, from quantifying the effects of human variants on disease likelihood, to predicting immune-escape mutations in viruses and designing novel biotherapeutic proteins. Deep generative models of protein sequences trained on multiple sequence alignments have been the most successful approaches so far to address these tasks. The performance of these methods is however contingent on the availability of sufficiently deep and diverse alignments for reliable training. Their potential scope is thus limited by the fact many protein families are hard, if not impossible, to align. Large language models trained on massive quantities of non-aligned protein sequences from diverse families address these problems and show potential to eventually bridge the performance gap. We introduce Tranception, a novel transformer architecture leveraging autoregressive predictions and retrieval of homologous sequences at inference to achieve state-of-the-art fitness prediction performance. Given its markedly higher performance on multiple mutants, robustness to shallow alignments and ability to score indels, our approach offers significant gain of scope over existing approaches. To enable more rigorous model testing across a broader range of protein families, we develop ProteinGym -- an extensive set of multiplexed assays of variant effects, substantially increasing both the number and diversity of assays compared to existing benchmarks.

## Tranception
Coming soon!

## ProteinGym
ProteinGym is an extensive set of Deep Mutational Scanning (DMS) assays curated to enable thorough comparisons of various mutation effect predictors indifferent regimes. ProteinGym is comprised of two benchmarks: 1) a substitution benchmark which consists of the experimental characterisation of ∼1.5M missense variants across 87 DMS assays 2) an indel benchmark that includes ∼300k mutants across 7 DMS assays.

Each processed file in each benchmark corresponds to a single DMS assay, and contains the following three variables:
- mutant (str): 
    - for the substitution benchmark, it describes the set of substitutions to apply on the reference sequence to obtain the mutated sequence (eg., A1P:D2N implies the amino acid 'A' at position 1 should be replaced by 'P', and 'D' at position 2 should be replaced by 'N')
    - for the indel benchmark, it corresponds to the full mutated sequence
- DMS_score (float): corresponds to the experimental measurement in the DMS assay. Across all assays, the higher the DMS_score value, the higher the fitness of the mutated protein
- DMS_score_bin (int): indicates whether the DMS_score is above the fitness cutoff (1 is fit, 0 is not fit)

Additionally, we provide reference files in the [ProteinGym folder](https://github.com/OATML-Markslab/Tranception/tree/main/ProteinGym) that give further details on each assay and contain in particular:
- The UniProt_ID of the corresponding protein, along with taxon and MSA depth category
- The target sequence (target_seq) used in the assay
- Details on how the DMS_score was created from the raw files and how it was binarized 

To download the substitution benchmark:
```
curl -o ProteinGym_substitutions.zip https://marks.hms.harvard.edu/ProteinGym/ProteinGym_substitutions.zip 
unzip ProteinGym_substitutions.zip
rm ProteinGym_substitutions.zip
```

Similarly, to download the indel benchmark:
```
curl -o ProteinGym_indels.zip https://marks.hms.harvard.edu/ProteinGym/ProteinGym_indels.zip
unzip ProteinGym_indels.zip
rm ProteinGym_indels.zip
```

The [ProteinGym folder](https://github.com/OATML-Markslab/Tranception/tree/main/ProteinGym) also includes detailed performance files for Tranception and other baselines on the two ProteinGym benchmarks.

## Multiple Sequence Alignments (MSAs)

The MSAs used to train alignment-based methods or used at inference in Tranception with retrieval and MSA Transformer may be downloaded as follows:
```
curl -o MSA_ProteinGym.zip https://marks.hms.harvard.edu/ProteinGym/MSA_ProteinGym.zip
unzip MSA_ProteinGym.zip
rm MSA_ProteinGym.zip
```

## Reference
If you use Tranception or ProteinGym in your work, please cite the following paper:
```
Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A., Marks, D.S., Gal, Y. (2022). Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval. ICML.
```
