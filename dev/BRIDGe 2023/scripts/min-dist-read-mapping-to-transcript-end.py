# for TLEN \=\= 0 and REF !\= *, get minimum distance of the read mapping site to the transcript end position of the same chromosome((

from gtfparse import read_gtf

df = read_gtf("indices/refdata-gex-GRCh38-2020-A/genes/genes.gtf")
df_transcripts = df[df["feature"] == "transcript"]
# print(head(df_transcripts))