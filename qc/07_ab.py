import hail as hl

#Import mt file
print("Reading file...")
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow_autosomal_gnomad_sampleqc_gt-variantqc.mt', overwrite=True)

#Annotate in (alternate) allele balance
print("Annotating with Allele Balance...")
mt = mt.annotate_entries(AB = mt.AD[1] / hl.sum(mt.AD))

#Annotate in "good" allele balance filters at a variant level
print("Annotating with Allele Balance filters at a variant level...")
mt = mt.annotate_rows(goodHetAlleleBalance = hl.agg.fraction(mt.GT.is_hom_ref() | (mt.GT.is_het() & ((mt.AB >= 0.3) & (mt.AB <= 0.7))) | mt.GT.is_hom_var()))
mt = mt.annotate_rows(goodAlleleBalance = hl.agg.fraction(((mt.GT.is_hom_ref()) & (mt.AB < 0.3)) | (mt.GT.is_het() & (mt.AB >= 0.3) & (mt.AB <= 0.7)) | (mt.GT.is_hom_var() & (mt.AB > 0.7))))

#Filter based on "good" alelle balance filters
print("Filtering for Allele Balance...")
mt = mt.filter_rows((mt.goodHetAlleleBalance > 0.9) & (mt.goodAlleleBalance > 0.5))

print(mt.count())

#Write to file
print("Writing to file...")
mt.write('gs://ibd-exomes/v36meta/ice.v36+ccdg_qc.mt', overwrite=True)