import hail as hl

#Import the initial vcf
print("Importing the table...")
mt = hl.read_matrix_table('gs://ibd-exomes/v35+ccdg_complete_qc.mt')

hl.export_vcf(mt, 'gs://ibd-exomes/v35_ccdg_complete_qc.vcf.bgz', parallel=None)