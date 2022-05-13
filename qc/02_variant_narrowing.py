import hail as hl

hl.init(log='/home/variant_narrowing.log')

mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/ice.v36+ccdg.mt')
#1.2M

mt = mt.naive_coalesce(1000)

#Variant narrowing:#############################################

variants = mt.rows()

#Annotate with variant type data via VEP
print("Annotating with VEP...")
variants = hl.vep(variants, "gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json")

#make a literal with all variant types of interest
vars_of_interest = hl.literal({"frameshift_variant", "inframe_deletion", "inframe_insertion", "stop_lost", "stop_gained", "start_lost", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "missense_variant", "synonymous_variant"})

#Filter for all variants that have a defined most severe consequence that is part of the literal containing all variant types of interest
print("Filtering for variants of interest...")
variants = variants.filter(hl.is_defined(variants.vep.most_severe_consequence))
variants = variants.annotate_globals(x = vars_of_interest)
variants = variants.filter(variants.x.contains(variants.vep.most_severe_consequence))

print("Variants left:")
print(variants.count())

#Annotate back in
mt = mt.annotate_rows(vep = variants[mt.row_key].vep)
mt = mt.filter_rows(hl.is_defined(mt.vep))

print("Writing to file...")
mt.write('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow.mt', overwrite=True)