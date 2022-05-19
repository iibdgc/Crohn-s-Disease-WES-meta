
import os
import subprocess
import hail as hl

pop_strings = ['AFR', 'AJ', 'AMR', 'EAS', 'FIN', 'SAS']
# 'NFE' later

for pop_string in pop_strings:
	for disease in ['cd', 'ibd', 'uc']:
		print('Reading in ' + pop_string + ' ' + disease + '+controls file...')
		mt = hl.read_matrix_table('gs://ibd-exomes/' + pop_string + '_' + disease + '+control.mt')

		print('Doing Variant QC...')
		mt = hl.variant_qc(mt, name='variant_qc')

		print('Making case and control tables...')
		if disease == 'cd':
			case_mt = mt.filter_cols(mt.has_cd, keep=True)
			control_mt = mt.filter_cols(mt.has_cd, keep=False)

			print('Annotating CaAC...')
			mt = mt.annotate_rows(caac = hl.agg.filter(mt.has_cd, hl.agg.sum(mt.GT.n_alt_alleles())))

			print('Annotating CoAC...')
			mt = mt.annotate_rows(coac = hl.agg.filter(~mt.has_cd, hl.agg.sum(mt.GT.n_alt_alleles())))

			print('Annotating CaNAC...')
			mt = mt.annotate_rows(canac = hl.agg.filter(mt.has_cd, hl.agg.sum(mt.GT.is_hom_ref() * 2 + mt.GT.is_het_ref())))

			print('Annotating CoNAC...')
			mt = mt.annotate_rows(conac = hl.agg.filter(~mt.has_cd, hl.agg.sum(mt.GT.is_hom_ref() * 2 + mt.GT.is_het_ref())))
		elif disease == 'ibd':
			case_mt = mt.filter_cols(mt.has_ibd, keep=True)
			control_mt = mt.filter_cols(mt.has_ibd, keep=False)

			print('Annotating CaAC...')
			mt = mt.annotate_rows(caac = hl.agg.filter(mt.has_ibd, hl.agg.sum(mt.GT.n_alt_alleles())))

			print('Annotating CoAC...')
			mt = mt.annotate_rows(coac = hl.agg.filter(~mt.has_ibd, hl.agg.sum(mt.GT.n_alt_alleles())))

			print('Annotating CaNAC...')
			mt = mt.annotate_rows(canac = hl.agg.filter(mt.has_ibd, hl.agg.sum(mt.GT.is_hom_ref() * 2 + mt.GT.is_het_ref())))

			print('Annotating CoNAC...')
			mt = mt.annotate_rows(conac = hl.agg.filter(~mt.has_ibd, hl.agg.sum(mt.GT.is_hom_ref() * 2 + mt.GT.is_het_ref())))
		elif disease == 'uc':
			case_mt = mt.filter_cols(mt.has_uc, keep=True)
			control_mt = mt.filter_cols(mt.has_uc, keep=False)

			print('Annotating CaAC...')
			mt = mt.annotate_rows(caac = hl.agg.filter(mt.has_uc, hl.agg.sum(mt.GT.n_alt_alleles())))

			print('Annotating CoAC...')
			mt = mt.annotate_rows(coac = hl.agg.filter(~mt.has_uc, hl.agg.sum(mt.GT.n_alt_alleles())))

			print('Annotating CaNAC...')
			mt = mt.annotate_rows(canac = hl.agg.filter(mt.has_uc, hl.agg.sum(mt.GT.is_hom_ref() * 2 + mt.GT.is_het_ref())))

			print('Annotating CoNAC...')
			mt = mt.annotate_rows(conac = hl.agg.filter(~mt.has_uc, hl.agg.sum(mt.GT.is_hom_ref() * 2 + mt.GT.is_het_ref())))

		print('Doing Variant QC on cases and controls separately...')
		case_mt = hl.variant_qc(case_mt, name='variant_qc')
		control_mt = hl.variant_qc(control_mt, name='variant_qc')

		print('Annotating case and control pHWE...')
		mt = mt.annotate_rows(case_phwe = case_mt.index_rows(mt.row_key).variant_qc.p_value_hwe)
		mt = mt.annotate_rows(control_phwe = control_mt.index_rows(mt.row_key).variant_qc.p_value_hwe)

		

		print('Filtering on CaAC and CoAC...')
		mt = mt.filter_rows((mt.caac > 0) | (mt.coac > 0))

		print('Performing Fisher Exact Test...')
		mt = mt.annotate_rows(fet = hl.fisher_exact_test(hl.int32(mt.caac), hl.int32(mt.canac), hl.int32(mt.coac), hl.int32(mt.conac)))

		print('Annotating Standard Error...')
		mt = mt.annotate_rows(OR = (((mt.caac + 0.5) * (mt.conac + 0.5)) / ((mt.coac + 0.5) * (mt.canac + 0.5))))
		mt = mt.annotate_rows(se = 1.0/(mt.caac + 0.5) + 1.0/(mt.canac + 0.5) + 1.0/(mt.coac + 0.5) + 1.0/(mt.conac + 0.5))

		print('Ordering by p-value...')
		mt = mt.annotate_rows(fet_p_value = mt.fet.p_value)
		sorted_ht = mt.rows().order_by('fet_p_value')

		print('Exporting to table...')
		sorted_ht.select(locus = sorted_ht.locus, alleles = sorted_ht.alleles, P = sorted_ht.fet_p_value, OR = sorted_ht.OR, SE = sorted_ht.se, CaAC = sorted_ht.caac, CaNAC = sorted_ht.canac, CoAC = sorted_ht.coac, CoNAC = sorted_ht.conac, maf = sorted_ht.variant_qc.AF, call_rate = sorted_ht.variant_qc.call_rate, mean_dp = sorted_ht.variant_qc.dp_stats.mean, case_phwe = sorted_ht.case_phwe, control_phwe = sorted_ht.control_phwe).export('gs://ibd-exomes/' + pop_string + '_' + disease + '_FET_results.tsv.bgz')