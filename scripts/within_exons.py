import pandas as pd
import os
def categorize_intervals(chrom, read_start, read_end, gene_id, cigar, reads, length_of_read,forward_exons_df):

    amb_result = []
    intr_result = []
    uns_result = []
    spl_result = []
    utr_3_result = []
    utr_5_result = []
    status = []

    if pd.notna(gene_id):
        same_gene = forward_exons_df[(forward_exons_df['chromosome'] == chrom) &
                                    (forward_exons_df['gene_id'] == gene_id)]

        if not same_gene.empty:
            same_exon = same_gene[same_gene['transcript'] == same_gene['next transcript']]
            if not same_exon.empty:
                same_exon_number = same_exon[(same_exon['exon number']) < (same_exon['last exon number'])]
                if not same_exon_number.empty:
                    between_exons = same_exon_number[(same_exon_number['end'] < read_start) & (same_exon_number['next exon start'] > read_start)]
                    if not between_exons.empty:
                        #for i in range(len(intronic)):
                        intronic = between_exons[(between_exons['next exon start'] > read_end)]
                        if not intronic.empty:
                            intr_exon_1 = intronic.iloc[0]['transcript']
                            intr_exon_2 = intronic.iloc[0]['next transcript']
                            intr_exon_number1 = intronic.iloc[0]['exon number']
                            intr_exon_number2 = intronic.iloc[0]['next exon number']
                            intr_result.append([[intr_exon_1, intr_exon_2], [intr_exon_number1, intr_exon_number2]])
                            status.append("I")
                        unspliced = between_exons[(between_exons['next exon start'] < read_end) & (between_exons['next exon end'] > read_end)]
                        if not unspliced.empty:
                            #for i in range(len(unspliced)):
                            uns_exon_1 = unspliced.iloc[0]['transcript']
                            uns_exon_2 = unspliced.iloc[0]['next transcript']
                            uns_exon_number1 = unspliced.iloc[0]['exon number']
                            uns_exon_number2 = unspliced.iloc[0]['next exon number']
                            uns_result.append([[uns_exon_1, uns_exon_2], [uns_exon_number1, uns_exon_number2]])
                            status.append("U")
                            if 'N' in cigar:
                                status.append('S')
                    spliced = same_exon_number[(same_exon['start'] <= read_start) & (same_exon_number['end'] >= read_start) & (same_exon_number['next exon start'] <= read_end) & (same_exon_number['next exon end'] >= read_end) |
                                            (same_exon['start'] <= read_start) & (same_exon_number['end'] >= read_start) & (same_exon_number['next exon start'] <= read_end) & (same_exon_number['next exon end'] <= read_end)]
                    if not spliced.empty:
                        #for i in range(len(spliced)):
                        spl_1 = spliced.iloc[0]['transcript']
                        spl_2 = spliced.iloc[0]['next transcript']
                        spl_number_1 = spliced.iloc[0]['exon number']
                        spl_number_2 = spliced.iloc[0]['next exon number']
                        spl_result.append([[spl_1, spl_2], [spl_number_1, spl_number_2]])
                        status.append("S")
                                
                    unspliced = same_exon_number[(same_exon_number['end'] >= read_start) & (same_exon_number['start'] <= read_start) & (same_exon_number['next exon start'] > read_end) & (same_exon_number['end'] < read_end)]
                    if not unspliced.empty:
                        #for i in range(len(unspliced)):
                        uns_exon_1 = unspliced.iloc[0]['transcript']
                        uns_exon_2 = unspliced.iloc[0]['next transcript']
                        uns_exon_number1 = unspliced.iloc[0]['exon number']
                        uns_exon_number2 = unspliced.iloc[0]['next exon number']
                        uns_result.append([[uns_exon_1, uns_exon_2], [uns_exon_number1, uns_exon_number2]])
                        status.append("U")
                        if 'N' in cigar:
                            status.append('S')
                    if 'N' in cigar:
                        status.clear()
                        status.append("S")
            ambiguous = same_gene[(same_gene['start'] <= read_start) & (same_gene['end'] >= read_end)]
            if not ambiguous.empty:
                #for i in range(len(ambiguous)):            
                amb_exon_name = ambiguous.iloc[0]['transcript']
                amb_exon_number = ambiguous.iloc[0]['exon number']
                amb_result.append([amb_exon_name, amb_exon_number])
                status.clear()
                status.append('A')
            if 'N' in cigar:
                        status.clear()
                        status.append("S")
    else:
        same_chrom = forward_exons_df[(forward_exons_df['chromosome'] == chrom)]
        if not same_chrom.empty:
            utr_5 = same_chrom[(same_chrom['first exon start'] >= read_end)]
            if not utr_5.empty:
                utr_5_exon_name = utr_5.iloc[0]['transcript']
                utr_5_exon_number = utr_5.iloc[0]['exon number']
                min_dist_5 = min(abs(utr_5['first exon start'] - read_end).min(), abs(utr_5['last exon end'] - read_start).min())
                dist_to_gene_start = utr_5.iloc[0]['distance to gene start']
                if min_dist_5 > (1000+dist_to_gene_start):
                    status.append("None")
                else:
                    utr_5_result.append([utr_5_exon_name, utr_5_exon_number])
                    utr_5_result.append(f"MinDist5:{min_dist_5}")
                    status.append("5T")
            else:
                min_dist_5 = float('inf')
            utr_3 = same_chrom[(same_chrom['last exon end'] <= read_start)]
            if not utr_3.empty:
                utr_3_exon_name = utr_3.iloc[0]['transcript']
                utr_3_exon_number = utr_3.iloc[0]['exon number']
                min_dist_3 = min(abs(utr_3['first exon start'] - read_start).min(), abs(utr_3['last exon end'] - read_start).min())
                utr_3_result.append([utr_3_exon_name, utr_3_exon_number])
                utr_3_result.append(f"MinDist3:{min_dist_3}")
                status.append("3T")
                if "None" in status:
                    status.clear()
                    status.append("3T")
                elif "5T" in status:
                    utr_3_result.clear()
                    status.clear()
                    status.append("5T")
            if 'N' in cigar:
                utr_3_result.clear()
                utr_5_result.clear()
                status.clear()
                status.append('S')

    return [chrom, (read_start, read_end), reads, length_of_read, gene_id, amb_result, intr_result, uns_result, spl_result,utr_3_result,utr_5_result, status]

folder_path = '/path/sc_mb'
file_list = os.listdir(folder_path)
output_folder = '/path/sc_mh_results'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
for i, file in enumerate(file_list):
    result_columns = ['chromosome', 'read interval', 'reads', 'length of read','gene id','ambiguous', 'intronic', 'unspliced', 'spliced', "3'UTR","5'UTR",'status' ]
    result_data = []

    file_path = os.path.join(folder_path, file)
    within_df = pd.read_csv(file_path)
    forward_exons_df = pd.read_csv('/path/m_forward_sort_ex.csv')

    for index, read_interval_row in within_df.iterrows():
        chrom = read_interval_row['chromosome']
        read_start = read_interval_row['read start']
        read_end = read_interval_row['read end']
        reads = read_interval_row['reads']
        gene_id = read_interval_row['gene_id']
        cigar = read_interval_row['cigar string']
        length_of_read = read_interval_row['length of read']

        result_data.append(categorize_intervals(chrom, read_start, read_end, gene_id,cigar, reads, length_of_read,forward_exons_df))

    result_df = pd.DataFrame(result_data, columns=result_columns)
    output_file_path = os.path.join(output_folder, f'/path/sc_mh_results/results_{i}.csv')
    result_df.to_csv(output_file_path, index=False)
