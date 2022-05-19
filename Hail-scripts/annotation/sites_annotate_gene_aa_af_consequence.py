import hail as hl

#Import mt file
print('Reading in vcfs...')
#disease_strings = ['CD', 'IBD', 'UC']
disease_strings = ['IBD']
disease_tables = []

for disease_string in disease_strings:
	disease_tables.append(hl.import_table('gs://ibd-exomes/' + disease_string + '_summ_stats_sites_only.tsv', types={
	'POS': hl.tint32}))

gnomad_ht = hl.read_table('gs://gnomad-public/release/2.1/ht/exomes/gnomad.exomes.r2.1.sites.ht')

for disease_string, disease_table in zip(disease_strings, disease_tables):
	disease_table = disease_table.key_by(locus = hl.locus(disease_table.CHROM, disease_table.POS), alleles = hl.array([disease_table.REF, disease_table.ALT]))

	disease_table = disease_table.join(gnomad_ht)
	print('Annotating with VEP....')
	disease_table = hl.vep(disease_table, "gs://hail-common/vep/vep/vep85-loftee-gcloud.json")
	print('Annotating specific fields...')
	disease_table = disease_table.annotate(gene = disease_table.vep.transcript_consequences.map(lambda x: x.gene_symbol))
	disease_table = disease_table.annotate(amino_a = disease_table.vep.transcript_consequences.map(lambda x: x.hgvsp))
	disease_table = disease_table.annotate(global_af = disease_table.freq[0].AF)
	disease_table.select(disease_table.gene, disease_table.amino_a, disease_table.global_af).export('gs://ibd-exomes/' + disease_string + '_gnomad_freq.tsv')