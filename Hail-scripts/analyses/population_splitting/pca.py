import hail as hl

def pc_project(
        mt: hl.MatrixTable,
        loadings_ht: hl.Table,
        loading_location: str = "loadings",
        af_location: str = "pca_af"
) -> hl.Table:
    """
    Projects samples in `mt` on pre-computed PCs.
    :param MatrixTable mt: MT containing the samples to project
    :param Table loadings_ht: HT containing the PCA loadings and allele frequencies used for the PCA
    :param str loading_location: Location of expression for loadings in `loadings_ht`
    :param str af_location: Location of expression for allele frequency in `loadings_ht`
    :return: Table with scores calculated from loadings in column `scores`
    :rtype: Table
    """
    n_variants = loadings_ht.count()
    mt = mt.annotate_rows(
        pca_loadings=loadings_ht[mt.row_key][loading_location],
        pca_af=loadings_ht[mt.row_key][af_location]
    )
    mt = mt.filter_rows(hl.is_defined(mt.pca_loadings) & hl.is_defined(mt.pca_af) &
                        (mt.pca_af > 0) & (mt.pca_af < 1))
    gt_norm = (mt.GT.n_alt_alleles() - 2 * mt.pca_af) / hl.sqrt(n_variants * 2 * mt.pca_af * (1 - mt.pca_af))
    mt = mt.annotate_cols(scores=hl.agg.array_sum(mt.pca_loadings * gt_norm))
    return mt.cols().select('scores')

hl.init(log='/home/pca.log')
#Import matrix table
print("Reading initial MT...")
mt = hl.read_matrix_table('gs://ibd-exomes/v36meta/v36+ccdg_082119.mt')
print(mt.count())
mt = hl.sample_qc(mt, name='sample_qc')

#Import 1000G
print("Reading in 1KG data...")
mt_1kg = hl.experimental.load_dataset(name='1000_Genomes_autosomes', version='phase_3', reference_genome='GRCh38')
print(mt_1kg.count())

print("Merging datasets...")

mt = mt.select_entries('GT')

mt = mt.select_cols(mt.sample_qc.call_rate, mt.sample_qc.n_called, mt.sample_qc.n_not_called, 
	mt.sample_qc.n_hom_ref, mt.sample_qc.n_het, mt.sample_qc.n_hom_var, mt.sample_qc.n_non_ref, mt.sample_qc.n_singleton, 
	mt.sample_qc.n_snp, mt.sample_qc.n_insertion, mt.sample_qc.n_deletion, mt.sample_qc.n_transition, mt.sample_qc.n_transversion, 
	mt.sample_qc.n_star, mt.sample_qc.r_ti_tv, mt.sample_qc.r_het_hom_var, mt.sample_qc.r_insertion_deletion)

mt_1kg = mt_1kg.select_cols(mt_1kg.sample_qc.call_rate, mt_1kg.sample_qc.n_called, mt_1kg.sample_qc.n_not_called,
	mt_1kg.sample_qc.n_hom_ref, mt_1kg.sample_qc.n_het, mt_1kg.sample_qc.n_hom_var, mt_1kg.sample_qc.n_non_ref, mt_1kg.sample_qc.n_singleton,
	mt_1kg.sample_qc.n_snp, mt_1kg.sample_qc.n_insertion, mt_1kg.sample_qc.n_deletion, mt_1kg.sample_qc.n_transition, mt_1kg.sample_qc.n_transversion, 
	mt_1kg.sample_qc.n_star, mt_1kg.sample_qc.r_ti_tv, mt_1kg.sample_qc.r_het_hom_var,mt_1kg.sample_qc.r_insertion_deletion)

mt = mt.union_cols(mt_1kg)
print(mt.count())
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg+1kg.mt', overwrite=True)

#Filter on call rate
print("Performing variant QC and filtering on call rate...")
mt = hl.variant_qc(mt, name='variant_qc')
mt = mt.filter_rows(mt.variant_qc.call_rate > 0.95)
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg+1kg_call_rate.mt', overwrite=True)

#Filter on gnomAD exome NFE MAF, Konrad's table contains NFE AF
# rg37 = hl.get_reference('GRCh37')
# rg38 = hl.get_reference('GRCh38')
# rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38) 

# print("Obtaining NFE AF in gnomAD...")
# nfe_ht = hl.read_table('gs://daly-finland-konrad/sum_stats_enrichments/fin_enriched_exomes_37.ht/')
# # nfe_ht = nfe_ht.naive_coalesce(1000)
# nfe_ht = nfe_ht.annotate(new_locus = hl.liftover(nfe_ht.locus, 'GRCh38'))
# nfe_ht = nfe_ht.filter(hl.is_defined(nfe_ht.new_locus))
# nfe_ht = nfe_ht.key_by(locus = nfe_ht.new_locus, alleles = nfe_ht.alleles)
# nfe_ht = nfe_ht.checkpoint('gs://ibd-exomes/v36meta/fin_enriched_exomes_38.ht', overwrite=True)

#Lifted over version, generated with above code
nfe_ht = hl.read_table('gs://ibd-exomes/v36meta/fin_enriched_exomes_38.ht')

mt = mt.annotate_rows(nfe_ht = nfe_ht[mt.row_key])
print("Filtering gnomAD NFE AF > 0.01...")
mt = mt.filter_rows(mt.nfe_ht.nfe.AF > 0.01)
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg+1kg_nfe.mt', overwrite=True)
print(mt.count())

print("Performing VEP...")
mt = hl.vep(mt, "gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json")
print("Filtering for synonymous variants...")
mt = mt.filter_rows(mt.vep.most_severe_consequence == "synonymous_variant")
print(mt.count())
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg+1kg_filtered.mt', overwrite=True)

#LD PRUNE
pruned_variants = hl.ld_prune(mt.GT)
print("Pruning...")
mt = mt.filter_rows(hl.is_defined(pruned_variants[mt.row_key]))
mt = mt.checkpoint('gs://ibd-exomes/v36meta/v36+ccdg+1kg_pruned.mt', overwrite=True)
print(mt.count())

#Pull out 1KG samples
print("Grabbing 1KG samples...")
kg_ht = hl.import_table('gs://ibd-exomes/v36meta/1kg_samples.tsv.bgz')
kg_ht = kg_ht.key_by('s')
kg_mt = mt.filter_cols(hl.is_defined(kg_ht[mt.s]), keep=True)

print(kg_mt.count())

#PCA on the 1kG
print("Performing PCA...")
eigenvalues, scores, loadings = hl.hwe_normalized_pca(kg_mt.GT, k=10, compute_loadings=True)
kg_mt = kg_mt.annotate_rows(pca_af=hl.agg.mean(kg_mt.GT.n_alt_alleles()) / 2)
loadings = loadings.annotate(pca_af=kg_mt.rows()[loadings.key].pca_af)

#Project loadings onto all samples
print("Projecting loadings onto all samples...")
projections = pc_project(mt, loadings)

print("Exporting results...")
projections.export('gs://ibd-exomes/v36meta/1kg_pca_projections_all_samples.tsv.bgz')