import hail as hl

#Import the mt file
print("Reading file...")
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow.mt')

#Filter out all rows that have chr X, Y, or MT
print("Filtering out X, Y, MT...")
mt = mt.filter_rows((mt.locus.contig != 'X') & (mt.locus.contig != 'Y') & (mt.locus.contig != 'MT'))

print(mt.count())

#Write to file
print("Writing to file...")
mt.write('gs://ibd-exomes/v36meta/ice.v36+ccdg_varnarrow_autosomal.mt', overwrite=True)