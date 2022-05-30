## Tranception detailed performance on ProteinGym

We provide detailed performance file for Tranception and baselines, for the 3 focus metrics in the paper (eg., Spearman's rank, AUC, MCC).

For the substitution benchmark, we have 3 sets of files:
- "All_models" prefixed-files, aggregated at the UniProt level --> performance metrics are computed on the subset of mutants scorable by all baselines. Performance is computed for each DMS, then aggregated at UniProt_ID level.
- "All_models" prefixed-files, with no aggregation (DMS level) --> performance metrics are computed on the subset of mutants scorable by all baselines. Performance is computed for each DMS, there is no further aggregation
- "All_mutants" prefixed-files --> performance metrics are computed on all mutants available in DMS assays and reported for the subset of models able to score all mutants.

For the indel benchmark, each DMS corresponds to a unique UniProt_ID so no separate aggregation is needed. Furthermore, all models that can score indels (Tranception & Wavenet) are also able to score all mutants.