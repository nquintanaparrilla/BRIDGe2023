from ncls import NCLS
from gtfparse import read_gtf

# query each read (fragment) on the tree to find which transcript it belongs to, and use the indexes

# take the start and end sites of the transcripts

import pysam
import tqdm
import pandas as pd
bamfile = "/fs/cbcb-scratch/zaza/Projects/ext_r1/star_out/GSE148504/Aligned.out.bam"
bamfile = pysam.AlignmentFile(bamfile, "rb")
rdr = bamfile.fetch(until_eof=True)
starts_query = []
ends_query = []
indexes_query = []
i = 0
for read in tqdm.tqdm(rdr):
    if read.is_mapped:
        """
        print(read.reference_start)
        print(read.reference_end)
        print(read.reference_id)
        """
        starts_query.append(read.reference_start)
        ends_query.append(read.reference_end)
        indexes_query.append(i) # the read number
        i+=1
        if i > 100:
            break # try only 100 reads

starts_query = pd.Series(starts_query)
ends_query = pd.Series(ends_query)
indexes_query = pd.Series(indexes_query)
query_df = pd.DataFrame({"Start": starts_query.values, "End": ends_query.values}, index=indexes_query.values)
print(query_df)
# exit()

# cannot reuse the index
# https://github.com/pyranges/ncls/issues/32
df = read_gtf("indices/refdata-gex-GRCh38-2020-A/genes/genes.gtf") # this is the longest step
df_exons = df[df["feature"] == "exon"]
# print(df_exons)
# this tree is for the transcripts = genes
ncls = NCLS(df_exons["start"], df_exons["end"], df_exons.index)
# builds decently quickly, but for hacking would be nice to use a jupyter notebook
print(ncls)
print(df_exons.index)
# exit()

# everything done in C/Cython; faster
l_idxs, r_idxs = ncls.all_overlaps_both(starts_query.values, ends_query.values, indexes_query.values)
print(l_idxs)
print(query_df.loc[l_idxs])

# genes.gtf contains the transcript location
# compute distance from mapping site to the transcript ending site
# the mapping site distance is in the sequence alignment

# find downstream polyA tracts, compute distance from each polyA tract to mapping site