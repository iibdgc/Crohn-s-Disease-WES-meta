import hail as hl

#Import mt file
print("Reading file...")
mt = hl.read_matrix_table('gs://ibd-exomes/v35_var_narrow_autosomal_gnomad_sample_qc_variant_qc_ab.mt')

#Filter for pHWE > 10^-4
print("Filtering for pHWE...")
mt = mt.filter_rows(mt.variant_qc.p_value_hwe > 10^-4)

print(mt.count())

#Write to file
print("Writing to file...")
mt.write('gs://ibd-exomes/v35_var_narrow_autosomal_gnomad_sample_qc_variant_qc_ab_phwe.mt', overwrite = True)