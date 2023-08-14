import pysam    #pysam can be used to parse through bam files
import pandas as pd
from pathlib import Path
import os

def get_read_location_info(read):
    reference = read.reference_name
    start = read.reference_start
    end = read.reference_end
    flag = read.flag
    cigar = read.cigarstring
    length_of_read = abs(start - end)
    return reference, start, end, flag, cigar, length_of_read

def main(bamfile_path):
    save = pysam.set_verbosity(0)
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    if not os.path.exists(bamfile_path + ".bai"):   #access bai file (bam file index)
        bamfile.index(bamfile_path)
    pysam.set_verbosity(save)

    data = {
        'reference': [],
        'start': [],
        'end': [],
        'flag': [],
        'strand': [],
        'cigar string': [],
        'length of read': []  
    }

    for read in bamfile.fetch():
        reference, start, end, flag, cigar, length_of_read = get_read_location_info(read)
        data['reference'].append(reference)
        data['start'].append(start)
        data['end'].append(end)
        data['flag'].append(flag)
        data['cigar string'].append(cigar)
        data['length of read'].append(length_of_read)

        # Determine the strand and add it to the 'strand' column
        if flag & 0x40:  # Check if the read is a forward read (0x40 = 64 in decimal)
            data['strand'].append('+')
        elif flag & 0x80:  # Check if the read is a reverse read (0x80 = 128 in decimal)
            data['strand'].append('-')
        else:
            data['strand'].append(None)  # Handle any other cases (optional)

    bamfile.close()

    df = pd.DataFrame(data)
    return df

if __name__ == "__main__":
    bamfile_path = "/path/mhcL1_Aligned.sortedByCoord.out.bam"
    df = main(bamfile_path)

df['reads'] = df.groupby(['reference', 'start', 'end', 'strand'])['reference'].transform('size')
df.drop_duplicates(subset=['reference', 'start', 'end', 'strand', 'reads'], inplace=True)
filepath = Path("/path/sc_mh_reads.csv")
filepath.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(filepath, index=False)

df = pd.read_csv("/path/sc_mh_reads.csv")

reverse_df = df[df['strand'] == '+']
forward_df = df[df['strand'] == '-']

forward_csv_filename = '/path/sc_mb_reads.csv'
reverse_csv_filename = '/path/barcodes_sc_mh_reads.csv'

reverse_df.to_csv(forward_csv_filename, index=False)
forward_df.to_csv(reverse_csv_filename, index=False)
