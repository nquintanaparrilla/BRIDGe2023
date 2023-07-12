import tqdm
import dnaio
import argparse

parser = argparse.ArgumentParser(description="Split reads for read1 analysis.")
parser.add_argument("ids", type=str, nargs="+",
                     help="list of sample ids to split")

args = parser.parse_args()
# TODO cli args to list multiple samples
# reads_ids = ["SRR8039659"]
reads_ids = args.ids
print(reads_ids)
# reads_list = ["fastq/{}_1.fastq.gz" for i in reads_ids]

# open the CSV file to get more info
PE_datasets = "PE datasets - Sheet2.tsv"
import pandas as pd
tsv = pd.read_csv(PE_datasets, sep="\t")

for read_id in reads_ids:
    read1 = "fastq/samples/{}_1.fastq.gz".format(read_id)
    read2 = "fastq/samples/{}_2.fastq.gz".format(read_id)
    print("Splitting", read1, "and", read2)
    # continue
    # TODO: better to extract the sample id from the path instead
    # no need to recompute this with every read!
    # moving it here makes it >50x faster!
    barcode_i = 12
    if tsv.loc[tsv['SRR'] == read_id]["chromium version"].item() == 'v2':
        barcode_i = 10
    umi_i = barcode_i + 16
    polyT_i = umi_i + 30

    with dnaio.open(read1, read2) as reader, \
         dnaio.open(read1.replace(".fastq.gz", "_barcode.fastq.gz"),
                    read1.replace(".fastq.gz", "_umi.fastq.gz"),
                    read1.replace(".fastq.gz", "_polyt.fastq.gz"),
                    read1.replace(".fastq.gz", "_rest.fastq.gz"),
                    read2.replace(".fastq.gz", "_rest.fastq.gz"),
                    mode="w") as writer1, \
         dnaio.open(read1.replace(".fastq.gz", "_barcode_small.fastq.gz"),
         read1.replace(".fastq.gz", "_umi_small.fastq.gz"),
                    read1.replace(".fastq.gz", "_polyt_small.fastq.gz"),
                    read1.replace(".fastq.gz", "_rest_small.fastq.gz"),
                    read2.replace(".fastq.gz", "_rest_small.fastq.gz"),
                    mode="w") as writer2:

        for record1, record2 in tqdm.tqdm(reader):
            barcode = record1[:barcode_i]
            umi = record1[barcode_i:umi_i]
            polyT = record1[umi_i:polyT_i]
            rest = record1[polyT_i:] # will serve as read_1 for alignment
            if len(rest) >= 10:
                writer1.write(barcode, umi, polyT, rest, record2)
                # TODO: also write the corresponding read pair in read_2
            else:
                # write small reads to another file
                writer2.write(barcode, umi, polyT, rest, record2)
    print("Done splitting", read1, "and", read2)
# TODO: paired-end data
# TODO: align with Star