import pyranges as pr
import pandas as pd

# read the gene names from the tsv

gene_names = pd.read_csv("selected-genes.tsv", header=None)[0]
print(gene_names)
# exit()
gr = pr.read_gtf("indices/refdata-gex-GRCh38-2020-A/genes/genes.gtf")

# TODO: assume gene_names is an array
# subset the gtf by genes matching the gene names

gr = gr[gr["Feature"] == "gene" & gr["gene_name"] in gene_names]
# sort to make get_sequence fasta (get it?)
gr = gr.sort()

# TODO the header should be gene_id + "-U" for unspliced transcripts
seq = pr.get_sequence(gr, "indices/refdata-gex-GRCh38-2020-A/fasta/genome.fa")
write_fasta(seq, "indices/refdata-gex-GRCh38-2020-A/fasta/genome-selected-genes.fasta")

def write_fasta(seq, f):
    with open(f, "w") as fw:
        nchars = 60
        for row in seq.itertuples():
            s = '\n'.join([row.Sequence[i:i+nchars]
                          for i in range(0, len(row.Sequence), nchars)])
            fw.write(f'>{row.transcript}\n{s}\n')

# append a string to the row.transcript name to indicate it is unspliced

# TODO the header should be transcript_id for spliced transcripts
seq_fasta_file = "indices/refdata-gex-GRCh38-2020-A/fasta/genome.fa"
seq_transcript = pr.get_transcript_sequence(gr, seq_fasta_file)

concatted = pr.concat([seq, seq_transcript])

transcripts_fasta_file = "indices/refdata-gex-GRCh38-2020-A/fasta/genome-selected-genes-transcripts.fasta"
write_fasta(seq, transcripts_fasta_file)

# 1. with seq: run starSOLO then call findSNR
# 2. with seq_transcript: run star aligner only