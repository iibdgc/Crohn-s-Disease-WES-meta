import hail as hl

#Import the initial ibd_mt containing all of the variants
print("Importing the variant table...")
mt = hl.read_matrix_table('gs://ibd-exomes/v35+ccdg_complete_qc.mt')

#Import PC table
print("Importing the PC table...")
pc_table = hl.import_table('gs://ibd-exomes/pc.tsv', delimiter='\t', key='SAMPLE', types={
	'PC1': hl.tfloat64, 
	'PC2': hl.tfloat64, 
	'PC3': hl.tfloat64, 
	'PC4': hl.tfloat64, 
	'PC5': hl.tfloat64, 
	'PC6': hl.tfloat64, 
	'PC7': hl.tfloat64, 
	'PC8': hl.tfloat64, 
	'PC9': hl.tfloat64, 
	'PC10': hl.tfloat64})

def read_list_data(input_file):
    sample_list = hl.import_table(input_file, delimiter='\t', key='SAMPLE_ID_IN_VCF')
    return sample_list

all_tables = []
table_files = ['gs://ibd-exomes/CD.tsv', 'gs://ibd-exomes/IBD.tsv', 'gs://ibd-exomes/UC.tsv', 'gs://ibd-exomes/CONTROL.tsv', 'gs://ibd-exomes/CUT.tsv']

print("Reading in disease status tables...")
for table_file in table_files:
	all_tables.append(read_list_data(table_file))

cd_table, ibd_table, uc_table, control_table, cut_table = all_tables[0], all_tables[1], all_tables[2], all_tables[3], all_tables[4]


def read_pop_data(input_file):
	pop_list = hl.import_table(input_file, no_header=True)
	pop_list = pop_list.key_by('f0')
	return pop_list

pop_tables = []
pop_files = ['gs://ibd-exomes/AFR_20181116.txt', 'gs://ibd-exomes/AJ_20181116.txt', 'gs://ibd-exomes/AMR_20181116.txt', 'gs://ibd-exomes/EAS_20181116.txt', 'gs://ibd-exomes/FIN_20181116.txt', 'gs://ibd-exomes/NFE_20181116.txt', 'gs://ibd-exomes/SAS_20181116.txt']

print("Reading in cohort tables...")
for pop_file in pop_files:
	pop_tables.append(read_pop_data(pop_file))

#Annotate the PCs into original file
print("Annotating PCs...")
mt = mt.annotate_cols(**pc_table[mt.s])

#We want to annotate in pHWE for the controls and filter out the ones that don't have adequate ones.
print("Annotating controls and filtering for pHWE...")
control_mt = mt.filter_cols(hl.is_defined(control_table[mt.s]), keep=True)

print("Before filtering pHWE:")
print(control_mt.count())
control_mt = hl.variant_qc(control_mt, name='variant_qc')
control_mt = control_mt.filter_rows(control_mt.variant_qc.p_value_hwe > 10^-4)
print("After filtering pHWE:")
print(control_mt.count())

#In order to verify that we do not include IBD cases as controls for the CD and UC tests, we make new tables.
print("Subsetting into test tables...")
cd_mt = mt.filter_cols(hl.is_defined(cd_table[mt.s]), keep=True)
ibd_mt = mt.filter_cols(hl.is_defined(cd_table[mt.s]) | hl.is_defined(ibd_table[mt.s]) | hl.is_defined(uc_table[mt.s]), keep=True)
uc_mt = mt.filter_cols(hl.is_defined(uc_table[mt.s]), keep=True)

print("Joining controls and cases for testing...")
cd_mt = cd_mt.union_cols(control_mt)
ibd_mt = ibd_mt.union_cols(control_mt)
uc_mt = uc_mt.union_cols(control_mt)

# print(cd_mt.count(), ibd_mt.count(), uc_mt.count())

#Exclude samples from CUT list
print("Excluding samples from cut list...")
cd_mt = cd_mt.filter_cols(hl.is_defined(cut_table[cd_mt.s]), keep=False)
ibd_mt = ibd_mt.filter_cols(hl.is_defined(cut_table[ibd_mt.s]), keep=False)
uc_mt = uc_mt.filter_cols(hl.is_defined(cut_table[uc_mt.s]), keep=False)

# print(cd_mt.count(), ibd_mt.count(), uc_mt.count())

#Annotate the cases
print("Annotating disease status...")
cd_mt = cd_mt.annotate_cols(has_cd = hl.is_defined(cd_table[cd_mt.s]))
ibd_mt = ibd_mt.annotate_cols(has_ibd = hl.is_defined(ibd_table[ibd_mt.s]))
uc_mt = uc_mt.annotate_cols(has_uc = hl.is_defined(uc_table[uc_mt.s]))

#Subsetting into populations
print("Subsetting into populations and writing tables to files...")
pop_strings = ['AFR', 'AJ', 'AMR', 'EAS', 'FIN', 'NFE', 'SAS']

for pop, pop_string in zip(pop_tables, pop_strings):
	pop_specific_cd_mt = cd_mt.filter_cols(hl.is_defined(pop[cd_mt.s]), keep=True)
	pop_specific_ibd_mt = ibd_mt.filter_cols(hl.is_defined(pop[ibd_mt.s]), keep=True)
	pop_specific_uc_mt = uc_mt.filter_cols(hl.is_defined(pop[uc_mt.s]), keep=True)

	print(pop_string, pop_specific_cd_mt.count(), pop_specific_ibd_mt.count(), pop_specific_uc_mt.count())

	pop_specific_cd_mt.write('gs://ibd-exomes/' + pop_string + '_cd+control.mt', overwrite=True)
	pop_specific_ibd_mt.write('gs://ibd-exomes/' + pop_string + '_ibd+control.mt', overwrite=True)
	pop_specific_uc_mt.write('gs://ibd-exomes/' + pop_string + '_uc+control.mt', overwrite=True)

# #for test in ['wald', 'lrt', 'firth']:
# for test in ['firth']:
# 	print("Conducting " + test + " tests...")

# 	# print("Crohn's disease...")
# 	# cd_result = hl.logistic_regression_rows(test=test, y=cd_mt.has_cd, x=cd_mt.GT.n_alt_alleles(), covariates=[1, cd_mt.PC1, cd_mt.PC2, cd_mt.PC3, cd_mt.PC4])
# 	# cd_result = cd_result.order_by(cd_result.p_value)
# 	# cd_result.export('gs://ibd-exomes/cd_' + test + '_guhan.tsv.bgz')

# 	# print("Inflammatory bowel disease...")
# 	# ibd_result = hl.logistic_regression_rows(test=test, y=ibd_mt.has_ibd, x=ibd_mt.GT.n_alt_alleles(), covariates=[1, ibd_mt.PC1, ibd_mt.PC2, ibd_mt.PC3, ibd_mt.PC4])
# 	# ibd_result = ibd_result.order_by(ibd_result.p_value)
# 	# ibd_result.export('gs://ibd-exomes/ibd_' + test + '_guhan.tsv.bgz')

# 	print("Ulcerative colitis...")
# 	uc_result = hl.logistic_regression_rows(test=test, y=uc_mt.has_uc, x=uc_mt.GT.n_alt_alleles(), covariates=[1, uc_mt.PC1, uc_mt.PC2, uc_mt.PC3, uc_mt.PC4])
# 	uc_result = uc_result.order_by(uc_result.p_value)
# 	uc_result.export('gs://ibd-exomes/uc_' + test + '_guhan.tsv.bgz')