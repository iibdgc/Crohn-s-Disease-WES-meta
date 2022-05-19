import hail as hl

#Import mt file
print("Reading in mts...")
pop_strings = ['AFR', 'AJ', 'AMR', 'EAS', 'FIN', 'NFE', 'SAS']
cd_tables = []
ibd_tables = []
uc_tables = []

for pop_string in pop_strings:
	cd_tables.append(hl.read_matrix_table('gs://ibd-exomes/' + pop_string + '_cd+control.mt'))
	ibd_tables.append(hl.read_matrix_table('gs://ibd-exomes/' + pop_string + '_ibd+control.mt'))
	uc_tables.append(hl.read_matrix_table('gs://ibd-exomes/' + pop_string + '_uc+control.mt'))

for test in ['wald', 'lrt', 'firth']:
	for i, pop_string in enumerate(pop_strings):
		print("Conducting " + test + " tests for " + pop_string + " population...")

		print("Crohn's disease...")
		cd_result = hl.logistic_regression_rows(test=test, y=cd_tables[i].has_cd, x=cd_tables[i].GT.n_alt_alleles(), covariates=[1, cd_tables[i].PC1, cd_tables[i].PC2, cd_tables[i].PC3, cd_tables[i].PC4])
		cd_result = cd_result.order_by(cd_result.p_value)
		cd_result.export('gs://ibd-exomes/' + pop_string + '_CD_' + test + '.tsv.bgz')

		print("Inflammatory bowel disease...")
		ibd_result = hl.logistic_regression_rows(test=test, y=ibd_tables[i].has_ibd, x=ibd_tables[i].GT.n_alt_alleles(), covariates=[1, ibd_tables[i].PC1, ibd_tables[i].PC2, ibd_tables[i].PC3, ibd_tables[i].PC4])
		ibd_result = ibd_result.order_by(ibd_result.p_value)
		ibd_result.export('gs://ibd-exomes/' + pop_string + '_IBD_' + test + '.tsv.bgz')

		print("Ulcerative colitis...")
		uc_result = hl.logistic_regression_rows(test=test, y=uc_tables[i].has_uc, x=uc_tables[i].GT.n_alt_alleles(), covariates=[1, uc_tables[i].PC1, uc_tables[i].PC2, uc_tables[i].PC3, uc_tables[i].PC4])
		uc_result = uc_result.order_by(uc_result.p_value)
		uc_result.export('gs://ibd-exomes/' + pop_string + '_UC_' + test + '.tsv.bgz')