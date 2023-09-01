# This version is also working in chunks of 1200

import pandas as pd
import os
def one_exon(forward_exons_df, chrom, read_start, read_end): # genes with one exon
    same_chrom = forward_exons_df[(forward_exons_df['chromosome'] == chrom) & (forward_exons_df['exon number'] == 1) & (forward_exons_df['last exon number']) == 1]
    ambiguous = same_chrom[(same_chrom['start'] <= read_start) & (same_chrom['end'] >= read_end)]
    if not ambiguous.empty:
        row = ambiguous.iloc[0]
        return row['transcript'], 'A'
    unspliced = same_chrom[(same_chrom['start'] <= read_start) & (same_chrom['end'] >= read_start) & (same_chrom['end'] < read_end) |
                           (same_chrom['start'] > read_start) & (same_chrom['start'] <= read_end) & (same_chrom['end'] >= read_end)]
    if not unspliced.empty:
        row = unspliced.iloc[0]
        return row['transcript'], 'U'
    return None, None

def middle_exons(forward_exons_df, chrom, read_start, read_end):
    same_chrom = forward_exons_df[(forward_exons_df['chromosome'] == chrom) & (forward_exons_df['exon number'] < forward_exons_df['next exon number']) & (forward_exons_df['last exon number'] != 1)]
    intronic = same_chrom[(same_chrom['end'] < read_start)  & (same_chrom['next exon start'] > read_start) & (same_chrom['next exon start'] > read_end)]
    if not intronic.empty:
        row = intronic.iloc[0]
        return row['transcript'], row['exon number'], row['next exon number'], 'I'
    ambiguous = ambiguous = same_chrom[(same_chrom['start'] <= read_start) & (same_chrom['end'] >= read_end)]
    if not ambiguous.empty:
        row = ambiguous.iloc[0]
        return row['transcript'], row['exon number'], row['next exon number'], 'A'
    unspliced1 = same_chrom[(same_chrom['start'] <= read_start) & (same_chrom['end'] >= read_start) & (same_chrom['end'] < read_end) & (same_chrom['next exon start'] > read_end)]
    if not unspliced1.empty:
        row = unspliced1.iloc[0]
        return row['transcript'], row['exon number'], row['next exon number'], 'U'
    
    unspliced2 = same_chrom[(same_chrom['start'] > read_start) & (same_chrom['start'] <= read_end) & (same_chrom['end'] >= read_end)]
    if not unspliced2.empty:
        row = unspliced2.iloc[0]
        return row['transcript'], row['exon number'], row['next exon number'], 'U'
    
    unspliced3 = same_chrom[(same_chrom['start'] <= read_start) & (same_chrom['end'] >= read_start) & (same_chrom['next exon start'] <= read_end) & (same_chrom['next exon end'] >= read_end)]
    if not unspliced3.empty:
        row = unspliced3.iloc[0]
        return row['transcript'], row['exon number'], row['next exon number'], 'U'
    
    return None, None, None, None

def spliced(cigar):
    status = 'S'
    if 'N' in cigar:
        return status
    return None

def utr_5(forward_exons_df, chrom, read_end):
    same_chrom = forward_exons_df[(forward_exons_df['chromosome'] == chrom)]
    first_exon = same_chrom[(same_chrom['exon number'] == 1) & (same_chrom['start'] > read_end)]
    if not first_exon.empty:
        min_dist_5 = min(abs(first_exon['first exon start'] - read_end).min(), abs(first_exon['last exon end'] - read_end).min())
        dist_to_gene_start = first_exon.iloc[0]['distance to gene start']
        if dist_to_gene_start == 0:
            dist_to_gene_start += 100 # made the cuttoff smaller
            if min_dist_5 <= dist_to_gene_start:
                return min_dist_5, f"MinDist5:{min_dist_5}", '5T'
            else:
                return min_dist_5, f"MinDist5:{min_dist_5}", 'F5'
        else:
            if min_dist_5 <= dist_to_gene_start:
                return min_dist_5, f"MinDist5:{min_dist_5}", '5T'
            else:
                return min_dist_5, f"MinDist5:{min_dist_5}", 'F5'
    return None, None, None
def utr_3(forward_exons_df, chrom, read_start):
    same_chrom = forward_exons_df[(forward_exons_df['chromosome'] == chrom)]
    last_exon = same_chrom[(same_chrom['exon number'] == same_chrom['last exon number']) & (same_chrom['end'] < read_start)]
    if not last_exon.empty:
        min_dist_3 = min(abs(last_exon['first exon start'] - read_start).min(), abs(last_exon['last exon end'] - read_start).min())
        dist_to_gene_end = last_exon.iloc[0]['distance to gene end']
        if dist_to_gene_end == 0:
            dist_to_gene_end += 360 #3.6 times bigger than 5'UTR
            if min_dist_3 <= dist_to_gene_end:
                return min_dist_3, f"MinDist3:{min_dist_3}", '3T'
            else:
                return min_dist_3, f"MinDist3:{min_dist_3}", 'F3'
        else:
            if min_dist_3 <= dist_to_gene_end:
                return min_dist_3, f"MinDist3:{min_dist_3}", '3T'
            else:
                return min_dist_3, f"MinDist3:{min_dist_3}", 'F3'
    return None, None, None

def closest_utr(min_dist_3, min_dist_5, utr_distance_3, utr_distance_5, utr_status_5, utr_status_3):
    min_dist_3 = min_dist_3 or float('inf')
    min_dist_5 = min_dist_5 or float('inf')
    if min_dist_3 < min_dist_5:
        return utr_distance_3, utr_status_3
    else:
        return utr_distance_5, utr_status_5
def categorize_intervals(chrom, read_start, read_end, cigar, reads, length_of_read, forward_exons_df):
    transcript0, one_exon_status = one_exon(forward_exons_df, chrom, read_start, read_end)
    transcript1, exon_number0, exon_number1, middle_exons_status = middle_exons(forward_exons_df, chrom, read_start, read_end)
    spliced_status = spliced(cigar)
    min_dist_5, utr_distance_5, utr_status_5 = utr_5(forward_exons_df, chrom, read_end)
    min_dist_3, utr_distance_3, utr_status_3 = utr_3(forward_exons_df, chrom, read_start)

    status = []
    result = []

    if one_exon_status:
        result.append(transcript0)
        status.append(one_exon_status)
        #print(index,'found within an exon')
    elif middle_exons_status:
        result.append([transcript1, exon_number0, exon_number1])
        status.append(middle_exons_status)
        #print(read_start, 'is intronic, ambiguous or unspliced')
    elif spliced_status:
        result.clear()
        status.clear()
        status.append(spliced_status)
        #print(read_start, 'is spliced')  
    else:
        if utr_status_5:
            result.append(utr_distance_5)
            status.append(utr_status_5)
        if utr_status_3:
            result.append(utr_distance_3)
            status.append(utr_status_3)
        distance, right_status = closest_utr(min_dist_3, min_dist_5, utr_distance_3, utr_distance_5, utr_status_5, utr_status_3)
        if right_status:
            result.clear()
            status.clear()
            result.append(distance)
            status.append(right_status)
            #print(index, 'found a closer UTR')

    return [chrom, (read_start, read_end), reads, length_of_read, result, status]
result_columns = ['chromosome', 'read interval', 'reads', 'length of read','result','status' ]
result_data = []

reads_df = pd.read_csv('/path/to/file.csv')
forward_exons_df = pd.read_csv('path/to/m_forward_sort_ex.csv')

for index, read_interval_row in reads_df.iterrows():
    chrom = read_interval_row['chromosome']
    read_start = read_interval_row['read start']
    read_end = read_interval_row['read end']
    reads = read_interval_row['reads']
    cigar = read_interval_row['cigar string']
    length_of_read = read_interval_row['length of read']

    result_data.append(categorize_intervals(chrom, read_start, read_end, cigar, reads, length_of_read,forward_exons_df))

result_df = pd.DataFrame(result_data, columns=result_columns)
result_df.to_csv('path/to/results.csv', index=False)
print('categorization is done')

