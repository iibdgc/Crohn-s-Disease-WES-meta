# Large-scale sequencing identifies multiple genes and rare variants associated with Crohn’s disease susceptibility

## Synopsis
IBD exome-wide assocaition statistics from the meta-analysis of two callsets: Nextera and Twist. We used the fixed-effect meta-analysis in [METAL](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2922887/). Column headers are self-explanatory. Allele2 is the tested/effect allele. 

## Code availability
Computer code used in this study:
- Hail (quality control, variants effect annotation; https://hail.is) Analysis scripts are available in folder "Hail-scripts" of this repository.
- SAIGE (variant-based and gene-based association test; https://github.com/weizhouUMICH/SAIGE)
- METAL (meta-analysis; https://genome.sph.umich.edu/wiki/METAL_Documentation)
- PLINK (ancestry assignment, IBD relatedness QC and sample heterozygosity QC; plink1.9, https://www.cog-genomics.org/plink/; sample level QC pipeline, https://github.com/Annefeng/PBK-QC-pipeline)

## Citation
[Sazonovs, Stevens, Venkataraman, Yuan et al., Nat Genet, 2022](https://doi.org/10.1038/s41588-022-01156-2)
