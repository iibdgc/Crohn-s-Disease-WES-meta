import hail as hl

#Import the mt file
print("Reading file...")
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow_autosomal.mt')

#gnomad quality filtering - read in tables
print("Reading in gNOMAD tables...")
mt2 = hl.read_table('gs://gnomad-public/release/2.1/ht/exomes/gnomad.exomes.r2.1.sites.ht/')
mt3 = hl.read_table('gs://gnomad-public/release/2.1/ht/genomes/gnomad.genomes.r2.1.sites.ht/')

print("Annotating gNOMAD exomes and genomes fields...")
mt = mt.annotate_rows(gnomad_exomes = mt2[mt.row_key])
mt = mt.annotate_rows(gnomad_genomes = mt3[mt.row_key])

mt = mt.annotate_rows(GEFilter = (mt.gnomad_exomes.filters.length() > 0) | mt.gnomad_exomes.filters.contains("RF"))
mt = mt.annotate_rows(GGFilter = (mt.gnomad_genomes.filters.length() > 0) | mt.gnomad_genomes.filters.contains("RF"))

print("Filtering for PASSing variants...")
mt = mt.filter_rows(mt.GGFilter | mt.GEFilter, keep=False)

print(mt.count())

#Write to file
print("Writing to file...")
mt.write('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow_autosomal_gnomad.mt', overwrite=True)