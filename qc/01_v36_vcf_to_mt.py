import hail as hl

hl.init(log='/home/vcf_import.log')

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