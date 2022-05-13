import hail as hl

# hailctl dataproc start guhancluster --vep GRCh38 --zone us-west1-b --num-workers 149 --worker-machine-type n1-highmem-16 --master-machine-type n1-highmem-16 --worker-boot-disk-size 100

# Read in MatrixTable
hl.init(log='/home/complete_qc.log')

#Import the initial vcf
print("Importing the VCFs...")
vcf = hl.import_vcf('gs://ibd-exomes/v36meta/raw_vcfs/*', force_bgz=True, reference_genome='GRCh38', find_replace=('nul', '.'))
print(vcf.count())
vcf = vcf.naive_coalesce(10000)
vcf = vcf.checkpoint('gs://ibd-exomes/v36meta/v36_082119_checkpt.mt', overwrite=True)

print("Splitting multiallelic sites...")
#Convert the file into mt (hail accessible) format
vcf = hl.split_multi_hts(vcf)
print(vcf.count())

#Write to file
print("Writing to file...")
vcf = vcf.checkpoint('gs://ibd-exomes/v36meta/v36_082119.mt', overwrite=True)

#Combine with CCDG
print("Importing CCDG...")
ccdg = hl.import_vcf('gs://ibd-exomes/v36meta/ccdg.vcf.bgz', reference_genome='GRCh38', force_bgz=True)
print(ccdg.count())
print("Writing to file...")
ccdg.checkpoint('gs://ibd-exomes/v36meta/ccdg.mt', overwrite=True)

print("Combining the datasets...")
mt = vcf.select_entries('GT', 'AD', 'DP', 'GQ', 'PL').union_cols(ccdg)
print(mt.count())

print("Writing combined dataset to file...")
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_082119.mt', overwrite=True)

# mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/v36+ccdg_082119.mt')
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/v36_082119.mt')
# mt = mt.naive_coalesce(10000)
print(mt.count())

#NFE AF Filter for Survival Analysis if needed##################
rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38) 

#Contains NFE AF
nfe_ht = hl.read_table('gs://daly-finland-konrad/sum_stats_enrichments/fin_enriched_exomes_37.ht/')
nfe_ht = nfe_ht.naive_coalesce(1000)
nfe_ht = nfe_ht.annotate(new_locus = hl.liftover(nfe_ht.locus, 'GRCh38'))
nfe_ht = nfe_ht.filter(hl.is_defined(nfe_ht.new_locus))
nfe_ht = nfe_ht.key_by(locus = nfe_ht.new_locus, alleles = nfe_ht.alleles)
mt = mt.annotate_rows(nfe_ht = nfe_ht[mt.row_key])

make a literal with all variant types of interest
vars_of_interest = hl.literal({"frameshift_variant", "inframe_deletion", "inframe_insertion", "stop_lost", "stop_gained", "start_lost", "splice_acceptor_variant", "splice_donor_variant", "splice_region_variant", "missense_variant", "synonymous_variant"})

#Variant narrowing:#############################################
variants = mt.rows()

#Annotate with variant type data via VEP
print("Annotating with VEP...")
variants = hl.vep(variants, "gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json")

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

variants = None

print("Checkpointing...")
# mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_varnarrow.mt', overwrite=True)
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36_varnarrow.mt', overwrite=True)

#Autosomal Filtration###########################################

print("Filtering out X, Y, MT...")
mt = mt.filter_rows((mt.locus.contig != 'X') & (mt.locus.contig != 'Y') & (mt.locus.contig != 'MT'))

print("Checkpointing...")
# mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_varnarrow_autosomal.mt', overwrite=True)
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36_varnarrow_autosomal.mt', overwrite=True)

print(mt.count())

#gNOMAD#########################################################

#gnomad quality filtering - read in tables######################
print("Reading in gNOMAD tables...")
mt2 = hl.read_table('gs://gnomad-public/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht/')
mt3 = hl.read_table('gs://gnomad-public/release/2.1.1/liftover_grch38/ht/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.ht/')

print("Annotating gNOMAD exomes and genomes fields...")
mt = mt.annotate_rows(gnomad_exomes = mt2[mt.row_key])
mt = mt.annotate_rows(gnomad_genomes = mt3[mt.row_key])

mt = mt.annotate_rows(GEFilter = (mt.gnomad_exomes.filters.length() > 0) | mt.gnomad_exomes.filters.contains("RF"))
mt = mt.annotate_rows(GGFilter = (mt.gnomad_genomes.filters.length() > 0) | mt.gnomad_genomes.filters.contains("RF"))

print("Filtering for PASSing variants...")
mt = mt.filter_rows(mt.GGFilter | mt.GEFilter, keep=False)

print("Checkpointing...")
# mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_varnarrow_autosomal_gnomad.mt', overwrite=True)
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36_varnarrow_autosomal_gnomad.mt', overwrite=True)
print(mt.count())

#Sample QC######################################################

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

#Filtering out Agilent Samples
print("Filtering out Agilent samples...")
agilent = hl.import_table('gs://ibd-exomes/v36meta/agilent_samples_removed.txt', no_header=True, delimiter='\t')
agilent = agilent.annotate(sample = agilent.f0).key_by('sample')
mt = mt.filter_cols(hl.is_defined(agilent[mt.s]), keep=False)

print(mt.count())

print("Checkpointing...")
# mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_varnarrow_autosomal_gnomad_sampleqc.mt', overwrite=True)
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36_varnarrow_autosomal_gnomad_sampleqc.mt', overwrite=True)
print(mt.count())

#Genotype QC####################################################

#Set GQ to missing if low
print("Doing Genotype QC...")
mt = mt.annotate_entries(GT = hl.or_missing(mt.GQ > 20, mt.GT))

Variant QC#####################################################

#Annotate with variant-level QC metrics
print("Annotating with variant-level QC metrics...")
mt = hl.variant_qc(mt, name='variant_qc')

#Filter for call rate > 0.95
print("Filtering out low call rate...")
mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95)

print(mt.count())

#Filter for allele count > 0
print("Filtering out low allele count...")
mt = mt.filter_rows(mt.variant_qc.AC[0] > 0)

print(mt.count())

#Filter for mean depth > 10
print("Filtering out low mean depth...")
mt = mt.filter_rows(mt.variant_qc.dp_stats.mean > 10)

print(mt.count())

print("Checkpointing...")
# mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_varnarrow_autosomal_gnomad_sampleqc_gt-variantqc.mt', overwrite=True)
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36_varnarrow_autosomal_gnomad_sampleqc_gt-variantqc.mt', overwrite=True)

# Filters for different calls in gnoMAD exomes and genomes######
gnomad_filters = hl.import_table('gs://ibd-exomes/v36meta/gnomadfilters.tsv', types={'p_value': hl.tfloat64, 'odds_ratio': hl.tfloat64, 'genomes_AF_NFE': hl.tfloat64, 'exomes_AF_NFE': hl.tfloat64, 'V': hl.tstr})

print("Gnomad Filters Table count:")
print(gnomad_filters.count())
gnomad_filters = gnomad_filters.annotate(CHROM = gnomad_filters.V.split(':')[0], POS = hl.int(gnomad_filters.V.split(':')[1]), REF = gnomad_filters.V.split(':')[2], ALT = gnomad_filters.V.split(':')[3])
gnomad_filters = gnomad_filters.annotate(locus = hl.locus(gnomad_filters.CHROM, gnomad_filters.POS), alleles = hl.array([gnomad_filters.REF, gnomad_filters.ALT]))
gnomad_filters = gnomad_filters.annotate(new_locus = hl.liftover(gnomad_filters.locus, 'GRCh38'))
gnomad_filters = gnomad_filters.filter(hl.is_defined(gnomad_filters.new_locus))

print("Gnomad Filters Table count after liftover:")
print(gnomad_filters.count())

gnomad_filters = gnomad_filters.key_by(locus = gnomad_filters.new_locus, alleles = gnomad_filters.alleles)

print("Annotating gNOMAD exome/genome filter...")
mt = mt.annotate_rows(gnomad_exome_genome_filter = gnomad_filters[mt.row_key])

# Filter out those called significantly differently in exomes and genomes
mt = mt.annotate_rows(GE_fail=(hl.is_defined(mt.gnomad_exome_genome_filter) & ((mt.gnomad_exome_genome_filter.p_value < 0.00001) & ((mt.gnomad_exome_genome_filter.odds_ratio < 0.6667) | (mt.gnomad_exome_genome_filter.odds_ratio > 1.5)))))
mt = mt.filter_rows(mt.GE_fail, keep=False)

print("Checkpointing...")
# mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg_varnarrow_autosomal_gnomad_sampleqc_gt-variantqc_ge.mt', overwrite=True)
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36_varnarrow_autosomal_gnomad_sampleqc_gt-variantqc_ge.mt', overwrite=True)

print(mt.count())

#AB QC###########################################################

#Annotate in (alternate) allele balance
print("Annotating with Allele Balance...")
mt = mt.annotate_entries(AB = mt.AD[1] / hl.sum(mt.AD))

#Annotate in "good" allele balance filters at a variant level
print("Annotating with Allele Balance filters at a variant level...")
mt = mt.annotate_rows(goodHetAlleleBalance = hl.agg.fraction(mt.GT.is_hom_ref() | (mt.GT.is_het() & ((mt.AB >= 0.3) & (mt.AB <= 0.7))) | mt.GT.is_hom_var()))
#mt = mt.annotate_rows(goodAlleleBalance = hl.agg.fraction(((mt.GT.is_hom_ref()) & (mt.AB < 0.3)) | (mt.GT.is_het() & ((mt.AB >= 0.3) & (mt.AB <= 0.7))) | (mt.GT.is_hom_var() & (mt.AB > 0.7))))

#Filter based on "good" alelle balance filters
print("Filtering for Allele Balance...")
mt = mt.filter_rows((mt.goodHetAlleleBalance > 0.9))# & (mt.goodAlleleBalance > 0.5))

print(mt.count())

#Write to file
print("Writing to file...")
# mt.write('gs://ibd-exomes/v36meta/v36+ccdg_qc.mt', overwrite=True)
mt.write('gs://ibd-exomes/v36meta/v36_qc.mt', overwrite=True)