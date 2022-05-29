## ProteinGym reference files

In the reference files, we provide detailed information about all DMS assays included in ProteinGym. There are two reference files: one for the substitution benchmark and one for the indel benchmark.

The meaning of each column in the ProteinGym reference files is provided below:
- DMS_id: Uniquely identifies each DMS assay in ProteinGym. It is obtained as the concatenation of the UniProt ID of the mutated protein, the first author name and the year of publication. If there are several datasets with the same characteristics, another defining attribute of the assay is added to preserve unicity.
- multi_AA_mutants: Indicates whether the DMS contains several amino acid mutations.
- MSA_filename: Name of the Multiple Sequence Alignment that the prediction models used to make predictions for this DMS.
- region_mutated: Region of the target protein that is mutated in the DMS.
- target_seq: Sequence of the target protein (reference sequence mutated in the assay).
- seq_len: Length of the target protein sequence.
- MSA_start: Locates the beginning of the first sequence in the MSA with respect to the target sequence. For example, if the MSA covers from position 10 to position 60 of the target sequence, then MSA_start is 10.
- MSA_end: Locates the end of the first sequence in the MSA with respect to the target sequence. For example, if the MSA covers from position 10 to position 60 of the target sequence, then MSA_end is 60.
- Bitscore: Bitscore threshold used to generate the alignment divided by the length of the target protein.
- Theta: Hamming distance cutoff for sequence re-weighting.
- num_seqs: Number of sequences in the Multiple Sequence Alignment (MSA) used in this work for this DMS.
- perc_cov: Percentage of positions of the MSA that had a coverage higher than 70% (less than 30% gaps).
- num_cov: Number of positions of the MSA that had a coverage higher than 70% (less than 30% gaps).
- N_eff: The effective number of sequences in the MSA defined as the sum of the different sequence weights.
- N_eff_L: Neff / num_cov.
- num_significant: Number of evolutionary couplings that are considered significant. Significance is defined by having more than 90% probability of belonging to the log-normal distribution in a Gaussian Mixture Model of normal and log-normal distributions.
- num_significant_L: num_significant / num_cov.
- DMS_binarization_cutoff_ProteinGym: Cutoff used to divide fitness scores into binary labels.
- DMS_binarization_method: Method used to decide the binarization cutoff.
- Total_mutant_number: Number of rows of the DMS in ProteinGym.
- singles (only for substitutions): Number of single amino acid substitutions in the DMS.
- multiples (only for substitutions): Number of multiple amino acid substitutions in the DMS.

## Raw DMS assays files

We additionally provide the raw, unprocessed DMS files. The corresponding reference files provide additional details on how we processed the raw files to create ProteinGym:
- DMS_phenotype_name: Name of the column in the raw DMS that we used as fitness score.
- DMS_directionality: Direction of the correlation between the DMS_phenotype column values and the protein fitness in the raw DMS (directly or inversely correlated). In any given DMS, the directionality is 1 if higher values of the measurement are associated with higher fitness, and -1 otherwise. For simplicity, we adjusted directionality in the final ProteinGym benchmarks so that a higher value of DMS_score is always associated with higher fitness.
- DMS_mutant_column: Name of the column in the raw DMS that indicates which mutants were assayed. In the final ProteinGym, this column is always called "mutant". 

To download substitution raw DMS files:

```
curl -o substitutions_raw_DMS.zip https://marks.hms.harvard.edu/ProteinGym/substitutions_raw_DMS.zip
unzip substitutions_raw_DMS.zip
rm substitutions_raw_DMS.zip
```

Similarly, to download indel raw DMS files:

```
curl -o indels_raw_DMS.zip https://marks.hms.harvard.edu/ProteinGym/indels_raw_DMS.zip
unzip indels_raw_DMS.zip
rm indels_raw_DMS.zip
```

