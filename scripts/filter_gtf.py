import pyranges as pr
import pandas as pd
from pathlib import Path

#FILTER GENES 
gtf = pr.read_gtf("/path/ref/mouse_pri_assembly_annotation.gtf", as_df=True)
gene_df = gtf[gtf["Feature"] == 'gene']
gene_df = gene_df[["Chromosome", "Feature", "Source", "gene_id", "Start", "End", "Strand"]]
# Rename the columns as needed
gene_df.columns = ["chromosome", "feature","source","gene_id","start", "end", "strand"]
gene_df = gene_df.reset_index(drop=True)
filepath = Path("/path/ref/m_genes.csv")
filepath.parent.mkdir(parents=True, exist_ok=True)
gene_df.to_csv(filepath)

input_file = '/path/ref/m_genes.csv'
print('saved genes')

# Read the CSV file into a pandas DataFrame
df = pd.read_csv(input_file)

# Separate the DataFrame into two based on the 'strand' column
forward_df = df[df['strand'] == '+']
reverse_df = df[df['strand'] == '-']

# Replace 'forward.csv' and 'reverse.csv' with the desired filenames for the separated files
forward_csv_filename = '/path/m_reverse_genes.csv'
reverse_csv_filename = '/path/ref/m_forward_genes.csv'

# Save the separated DataFrames to new CSV files
forward_df.to_csv(forward_csv_filename, index=False)
reverse_df.to_csv(reverse_csv_filename, index=False)
print('saved forward and reverse genes')


# FILTER EXONS

exon_df = gtf[gtf["Feature"] == "exon"]
# Select the required columns
exon_df = exon_df[["Chromosome", "Feature", "gene_id","exon_number","transcript_name", "Start", "End", "Strand"]]
# Rename the columns as needed
exon_df.columns = ["chromosome", "feature","gene_id","exon number", "transcript", "start", "end", "strand"]
exon_df = exon_df.reset_index(drop=True)
exon_df.to_csv("/path/ref/m_exons.csv")
print("saved exons")

input_file1 = '/path/ref/m_exons.csv'
forward = pd.read_csv(input_file1)
forward_df = forward[forward['strand'] == '+']
forward_df.to_csv('/path/ref/m_forward_ex.csv', index=False)
print('saved forward exons')

input_file2 = '/path/ref/m_exons.csv'
reverse = pd.read_csv(input_file2)
reverse_df = reverse[reverse['strand'] == '-']
reverse_df.to_csv('/path/ref/m_reverse_ex.csv', index=False)
print("saved reverse exons")

df = pd.read_csv('/home/natalia/dev/BRIDGe2023-main/dev/forward_sort_ex.csv')

length_col = []
for index, exon_row in df.iterrows():
    start = exon_row['start']
    end = exon_row['end']
    length_of_exon = abs(start-end)

    length_col.append(length_of_exon)
df['length of exon'] = length_col

#Sort to get accurate next exon information
df.sort_values(by=['chromosome', 'gene_id','transcript', 'exon number'], inplace=True)

grouped = df.groupby(['gene_id'])

# Calculate 'next exon start' and 'next exon end' using shift
df['next exon number'] = grouped['exon number'].shift(-1)
df['next transcript'] = grouped['transcript'].shift(-1)
df['next exon start'] = grouped['start'].shift(-1)
df['next exon end'] = grouped['end'].shift(-1)

max_exon_info = df.groupby("gene_id").agg({"exon number": "max",  "start": "min", "end": "max"}).reset_index()

# Merge the max_exon_number DataFrame back into the original DataFrame
df = df.merge(max_exon_info, on="gene_id", suffixes=("", "_"))

# Rename the new column to indicate it's the last exon number
df.rename(columns={"exon number_": "last exon number", "start_":"first exon start", "end_": "last exon end"}, inplace=True)

# Save the updated DataFrame back to the CSV file
df.to_csv('/path/forward_sort_ex.csv', index=False)

# Get distance from the first exon to the start of the gene and the last exon to the end of the gene
df = pd.read_csv('/path/forward_sort_ex.csv')
gdf = pd.read_csv('/path/m_forward_genes.csv')
grouped_df = df.groupby('gene_id')

dist_to_start = []
dist_to_end = []

for index, gene_row in gdf.iterrows():
    start = gene_row['start']
    end = gene_row['end']
    gene_id = gene_row['gene_id']

    # Find the first exon start and last exon end for the current 'gene_id'
    first_exon_start = grouped_df.get_group(gene_id)['start'].min()
    last_exon_end = grouped_df.get_group(gene_id)['end'].max()

    # Calculate the distances to the gene start and gene end
    dist_to_gene_start = abs(start - first_exon_start)
    dist_to_gene_end = abs(end - last_exon_end)

    dist_to_start.append(dist_to_gene_start)
    dist_to_end.append(dist_to_gene_end)

# Add the 'distance to gene start' and 'distance to gene end' columns to the 'df' DataFrame
gdf['distance to gene start'] = dist_to_start
gdf['distance to gene end'] = dist_to_end

merged_df = pd.merge(df, gdf[['gene_id', 'distance to gene start', 'distance to gene end']], on='gene_id', how='left')

merged_df.to_csv('/path/forward_sort_ex.csv', index=False)
