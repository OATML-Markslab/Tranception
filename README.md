# Tranception

This is the official code repository for the paper "Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval". This project is a joint collaboration between the [Marks lab](https://www.deboramarkslab.com/) and the [OATML group](https://oatml.cs.ox.ac.uk/).

## Abstract
The ability to accurately model the fitness landscape of protein sequences is critical to a wide range of applications, from quantifying the effects of human variants on disease likelihood, to predicting immune-escape mutations in viruses and designing novel biotherapeutic proteins. Deep generative models of protein sequences trained on multiple sequence alignments have been the most successful approaches so far to address these tasks. The performance of these methods is however contingent on the availability of sufficiently deep and diverse alignments for reliable training. Their potential scope is thus limited by the fact many protein families are hard, if not impossible, to align. Large language models trained on massive quantities of non-aligned protein sequences from diverse families address these problems and show potential to eventually bridge the performance gap. We introduce Tranception, a novel transformer architecture leveraging autoregressive predictions and retrieval of homologous sequences at inference to achieve state-of-the-art fitness prediction performance. Given its markedly higher performance on multiple mutants, robustness to shallow alignments and ability to score indels, our approach offers significant gain of scope over existing approaches. To enable more rigorous model testing across a broader range of protein families, we develop ProteinGym -- an extensive set of multiplexed assays of variant effects, substantially increasing both the number and diversity of assays compared to existing benchmarks.

## Tranception
Coming soon!

## ProteinGym
ProteinGym is an extensive set of Deep Mutational Scanning (DMS) assays curated to enable thorough comparisons of various mutation effect predictors indifferent regimes. ProteinGym is comprised of two benchmarks: 1) a substitution benchmark which consists of the experimental characterisation of ∼1.5M missense variants across 87 DMS assays 2) an indel benchmark that includes ∼300k mutants across 7 DMS assays.

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

## Reference
If you use Tranception or ProteinGym in your work, please cite the following paper:
```
Notin, P., Dias, M., Frazer, J., Marchena-Hurtado, J., Gomez, A., Marks, D.S., Gal, Y. (2022). Tranception: protein fitness prediction with autoregressive transformers and inference-time retrieval. ICML.
```
