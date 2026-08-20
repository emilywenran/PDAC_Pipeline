import gzip
import pandas as pd

wgcna_file = "/Users/emily/thesis/PDAC/05_results/tables/phase8b_wgcna_module_assignments.tsv.gz"
with gzip.open(wgcna_file, 'rt') as f:
    df = pd.read_csv(f, sep='\t')

print("WGCNA Module Gene Counts:")
print(df['module'].value_counts())
