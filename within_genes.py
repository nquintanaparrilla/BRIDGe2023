import pandas as pd
import os

chunk_size = 1200  # Set your desired chunk size
input_file = '/path/to/filtered/reads/sc_mb_reads.csv'
output_folder = '/path/to/folder/results/sc_mb'

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

reader = pd.read_csv(input_file, chunksize=chunk_size)
header = True  # Include header in each chunk

for i, chunk in enumerate(reader):
    output_file = os.path.join(output_folder,f'sc_mb_{i}.csv')
    chunk.to_csv(output_file, index=False, header=header)

# Assign genes to each read
def within_gene(chrom, read_start, read_end, forward_genes_df):
    same_gene = forward_genes_df[(forward_genes_df['chromosome'] == chrom) & 
                                 (forward_genes_df['start'] <= read_start) &
                                 (forward_genes_df['end'] >= read_end)]
    
    if not same_gene.empty:
        gene_id = same_gene.iloc[0]['gene_id']
        return gene_id
    else:
        return None
    
folder_path = '/path/results/sc_mb'
file_list = os.listdir(folder_path)
output_folder = '/path/results/sc_mb_assigned'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
for i, file in enumerate(file_list):

    file_path = os.path.join(folder_path, file)
    forward_genes_df = pd.read_csv('/path/m_forward_genes.csv')
    reads_df = pd.read_csv(folder_path) 

    result_columns = ['chromosome', 'read start', 'read end', 'reads','length of read', 'cigar string','gene_id']
    result_data = []                                       
    for index, read_row in reads_df.iterrows():
        chrom = read_row['reference']
        read_start = read_row['start']
        read_end = read_row['end']
        reads = read_row['reads']
        length_of_read = read_row['length of read']
        cigar = read_row['cigar string']

        gene_id = within_gene(chrom, read_start, read_end, forward_genes_df)

        result_data.append([chrom, read_start, read_end, reads, length_of_read, cigar, gene_id])

    result_df = pd.DataFrame(result_data, columns=result_columns)
    output_file_path = os.path.join(output_folder, f'/path/results/sc_mb_assigned/sc_mh_{i}.csv')
    result_df.to_csv(output_file_path, index=False)

