import hail as hl

#Import the mt file
print("Reading file...")
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow_autosomal_gnomad.mt')

#Annotate with sample-level QC metrics
print("Annotating with sample-level QC metrics...")
mt = hl.sample_qc(mt, name='sample_qc')

#Filter to have low singleton count
print("Filtering out high singleton counts...")
mt = mt.filter_cols(mt.sample_qc.n_singleton < 500)

print(mt.count())

#Filter mean GQ > 40
print("Filtering out low mean GQ...")
mt = mt.filter_cols(mt.sample_qc.gq_stats.mean > 40)

print(mt.count())

#Filter call rate > 0.9
print("Filtering out low call rate...")
mt = mt.filter_cols(mt.sample_qc.call_rate > 0.9)

print(mt.count())

#Write to file
print("Writing to file...")
mt.write('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow_autosomal_gnomad_sampleqc.mt', overwrite=True)