# use pyranges
# ncls is the backend for that

import pyranges as pr
df = pr.read_gtf("indices/refdata-gex-GRCh38-2020-A/genes/genes.gtf", as_df=True)
df_exons = df[df["Feature"] == "exon"]
# print(df_exons)
transcripts_gr = pr.PyRanges(df=df_exons)

bamfile = "/fs/cbcb-scratch/zaza/Projects/ext_r1/star_out/GSE148504/Aligned.out.top.sam"
reads_gr = pr.read_bam(bamfile) # by default filters unmapped reads, https://pyranges.readthedocs.i
o/en/latest/autoapi/pyranges/index.html#pyranges.read_bam
# query each read (fragment) on the tree to find which transcript it belongs to, and use the indexe
s

print(transcripts_gr)
print(reads_gr)

# query must be a dictionary
grs = { "transcripts": transcripts_gr, "reads": reads_gr }

# do we want the set intersection or count_overlaps? I think intersection, right?
# no we want count_overlaps
res = pr.count_overlaps(grs) # .sort(by="transcripts")
# res = reads_gr.intersect(transcripts_gr)
print(res)
exit()




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
for read in tqdm.tqdm(rdr):
    starts_query.append(read.reference_start)
    ends_query.append(read.reference_end)
    indexes_query.append(read.reference_id) # numerical
    break # try only 1 read
starts_query = pd.Series(starts_query)
ends_query = pd.Series(ends_query)
indexes_query = pd.Series(indexes_query)
query_df = pd.DataFrame({"Start": starts_query.values, "End": ends_query.values}, index=indexes_query.values)
print(query_df)
# exit()

# cannot reuse the index
# https://github.com/pyranges/ncls/issues/32


# print(df_exons)
# this tree is for the transcripts = genes
# ncls = NCLS(df_exons["start"], df_exons["end"], df_exons.index)
# builds decently quickly, but for hacking would be nice to use a jupyter notebook
# print(ncls)
print(df_exons.index)
exit()

# everything done in C/Cython; faster
l_idxs, r_idxs = ncls.all_overlaps_both(starts_query.values, ends_query.values, indexes_query.values)
print(l_idxs)
print(query_df.loc[l_idxs])

# genes.gtf contains the transcript location
# compute distance from mapping site to the transcript ending site
# the mapping site distance is in the sequence alignment

# find downstream polyA tracts, compute distance from each polyA tract to mapping site