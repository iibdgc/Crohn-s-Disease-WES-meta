import hail as hl

#Import the initial vcf
print("Importing the table...")
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/v36_082119.mt')
print(mt.count())

positions = [62834058, 62841438, 62841450, 62819215]
refs = ["TTG", "G", "A", "C"]
alts = ["T", "A", "G", "CT"]

for position, ref, alt in zip(positions, refs, alts):
	var_mt = mt.filter_rows((mt.locus == hl.locus("chr1", position, reference_genome="GRCh38")) & (mt.alleles == hl.array([ref, alt])))
	var_mt = var_mt.filter_cols(hl.agg.count_where(var_mt.GT.is_non_ref()) > 0)
	print(var_mt.count())
	cols = var_mt.cols()
	cols.select().export('gs://ibd-exomes/v36meta/chr1_' + str(position) + '_' + ref + '_' + alt + '_carriers.tsv')