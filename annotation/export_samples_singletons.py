import hail as hl

#Import the initial vcf
print("Importing the table...")
mt = hl.read_matrix_table('gs://ibd-exomes/v35_complete_qc.mt')
cols = mt.cols()
cols.select(cols.sample_qc.n_singleton).export('gs://ibd-exomes/v35_samples_singletons_post_qc.mt')