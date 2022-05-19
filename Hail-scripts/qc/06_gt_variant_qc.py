import hail as hl

#Import the mt file
print("Reading file...")
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow_autosomal_gnomad_sampleqc.mt')

#Genotype QC####################################################

#Set GQ to missing if low
print("Doing Genotype QC...")
mt = mt.annotate_entries(GT = hl.or_missing(mt.GQ > 20, mt.GT))

#Variant QC#####################################################

#Annotate with variant-level QC metrics
print("Annotating with variant-level QC metrics...")
mt = hl.variant_qc(mt, name='variant_qc')

#Filter for call rate > 0.95
print("Filtering out low call rate...")
mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95)

#Filter for allele count > 0
print("Filtering out low allele count...")
mt = mt.filter_rows(mt.variant_qc.AC[0] > 0)

#Filter for mean depth > 10
print("Filtering out low mean depth...")
mt = mt.filter_rows(mt.variant_qc.dp_stats.mean > 10)

print(mt.count())

#Write to file
print("Writing to file...")
mt.write('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow_autosomal_gnomad_sampleqc_gt-variantqc.mt', overwrite=True)