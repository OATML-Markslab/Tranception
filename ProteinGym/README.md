
## ProteinGym reference file explanations

Next, we detail the meaning of each column in the DMS mapping spreadsheet. Certain columns (indicated between brackets) only apply to original datasets, and not to the final datasets, while certain other columns only apply to substitions and not to indels.

DMS_id: Maps every DMS assay in a unique manner. It comprises the UniProt ID of the mutated protein, the first author name and the year of publication. If there are several datasets with those same characteristics, another defining word is added.
DMS_filename (only for original datasets): Name of the DMS file in the set of original datasets. Note that these are not DMS filenames in ProteinGym, since ProteinGym filenames take the form of DMS_id.csv.
in_DeepSequence_benchmark: indicates whether the DMS was present in the previous DeepSequence benchmark set.
multi_AA_mutants: indicates whether the DMS contains multi amino acid mutants (more than singles).
MSA_filename: name of the Multiple Sequence Alignment that the prediction models used to make predictions for this DMS.
DMS_phenotype_name (only for original datasets): name of the column in the original dataset that we used as fitness score.
DMS_directionality (only for original datasets): Direction of the correlation between the DMS_phenotype column values and the protein fitness in the original datasets (directly or inversely correlated). In any given DMS, if the higher the phenotype value, the more fit the mutant is, then directionality is 1. On the other hand, if the lower the phenotype value, the less fit the mutation is, directionality is -1. In the final ProteinGym, all DMS_directionalities are 1 because we multiplied all fitness scores by -1 if the directionality was -1, thereby inverting its directionality.
DMS_mutant_column (only for original datasets): 
region_mutated: Region of the target protein that is mutated in the DMS.
target_seq: sequence of the target protein.
seq_len: Lenght of the target protein sequence.
MSA_start: Locates the beginning of the first sequence in the MSA with respect to the target sequence. For example, if the MSA covers from position 10 to position 60 of the target sequence, then MSA_start is 10.
MSA_end: Locates the end of the first sequence in the MSA with respect to the target sequence. For example, if the MSA covers from position 10 to position 60 of the target sequence, then MSA_start is 60.
Bitscore: Bitscore threshold used to generate the alignment divided by the length of the target protein.
Theta: Hamming distance cutoff for sequence re-weighting.
num_seqs: Number of sequences in the Multiple Sequence Alignment (MSA) used in this work for this DMS.
perc_cov: Percentage of positions of the MSA that had a coverage higher than 70% (less than 30% gaps).
num_cov: Number of positions of the MSA that had a coverage higher than 70% (less than 30% gaps).
N_eff: The effective number of sequences in the MSA defined as the sum of the different sequence weights.
N_eff_L: Neff/num_cov
num_significant: Number of evolutionary couplings that are considered significant. Significance is defined by having more than 90% probability of belonging to the log-normal distribution in a Gaussian Mixture Model of normal and log-normal distributions.
num_significant_L: num_significant/num_cov.
DMS_binarization_cutoff_ProteinGym: Cutoff used to divide fitness scores into binary labels.
DMS_binarization_method: Method used to decide the binarization cutoff.
Total_mutant_number: Number of rows of the DMS in ProteinGym.
singles (only for substitutions): number of single amino acid substitutions in the DMS.
doubles (only for substitutions): number of double amino acid substitutions in the DMS.
multiples (only for substitutions): number of multiple amino acid substitutions in the DMS.
