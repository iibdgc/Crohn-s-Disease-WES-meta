import pandas as pd

uniq_prism = pd.read_table('unique_prism_ids.txt')

for var in ['chr1_62819215_C_CT', 'chr1_62834058_TTG_T', 'chr1_62841438_G_A', 'chr1_62841450_A_G']:
	var_table = pd.read_table(var + '_carriers.tsv')
	var_table.merge(uniq_prism).to_csv(var + '_prism_carriers.tsv', sep='\t', index=False)