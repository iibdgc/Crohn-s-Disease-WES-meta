import hail as hl
import pandas as pd
import matplotlib.pyplot as plt
from pprint import pprint
from bokeh.io import output_notebook, show, save
from bokeh.layouts import gridplot
from bokeh.models import Span
import hail.expr.aggregators as agg
from bokeh.plotting import figure, output_file
import numpy as np

#Import summary stats
print("Importing the summary stat table...")
#ss_table = hl.import_table('/Users/Guhan/Downloads/AFR_IBD_wald.tsv', delimiter='\t', key='locus', types={'alleles': hl.tstr, 'beta': hl.tfloat64, 'standard_error': hl.tfloat64, 'z_stat': hl.tfloat64, 'p_value': hl.tfloat64, 'fit': hl.tstr})

ss_df = pd.read_csv('~/Downloads/AFR_IBD_wald.tsv', sep='\t')


p = hl.plot.qq(list(ss_df[ss_df.p_value.notna()]['p_value']))
#p = hl.plot.qq(ss_table.p_value)
show(p)
