import argparse
import pysam
import numpy as np
import matplotlib.pyplot as plt
import os

parser = argparse.ArgumentParser('description='"Count and summarize non-zero length reads.")
parser.add_argument("bam", type=str, help="analyze this bam file")
parser.add_argument("id", type=str, help="sample id")
parser.add_argument("savedir", type=str, help="where to save the results")

args = parser.parse_args()
# TODO cli args to list multiple samples
# TODO: maybe this should have been a Jupyter notebook so we don't continually read/reread the data

# reads_ids = ["SRR8039659"]
bamfile = args.bam
sid = args.id
savedir = args.savedir

print(args.bam)

bamfile = pysam.AlignmentFile(bamfile, "rb")

zerotlens = []
tlens = []
import tqdm
i=0
rdr = bamfile.fetch(until_eof=True)
for read in tqdm.tqdm(rdr):
    # check tlen > 0
    # print(read.reference_name)
#    print(read.is_read1)
#    print(next(rdr).is_read2)
#    print(read.is_read1)
#    print(read.is_read2)
#    print(read.get_aligned_pairs())
#    if i > 10: exit()
#    i+=1
    if read.template_length > 0:
#         if read.template_length < max(len(
        tlens.append(read.template_length)
        # TODO: which fields are these?
#    elif read.template_length == 0 and read.next_reference_name() == "*":
#        zerotlens.append([read.next_reference_name(), read.get_reference_positions()])

# exit()

if len(tlens) == 0:
    exit("tlens should not be empty")

tlens = np.array(tlens)
np.save("tlens.npy", tlens)
# counts = np.bincount(tlens)
# print(np.argmax(counts))
# exit()(
tlens_80_100 = tlens[(80 < tlens) & (tlens < 100)]
tlens_100_500 = tlens[(100 < tlens) & (tlens < 500)]
tlens_small = tlens[tlens < 500]
# print(tlens.size)
# print(tlens_small.size)
tlens_80_100_ratio = tlens_80_100.size / tlens.size
tlens_100_500_ratio = tlens_100_500.size / tlens.size
tlens_small_ratio = tlens_small.size / tlens.size
tlens_medium = tlens[(500 < tlens) & (tlens < 3000)]
tlens_medium_ratio = tlens_medium.size / tlens.size
tlens_large = tlens[3000 < tlens]
tlens_large_ratio = tlens_large.size / tlens.size

print(tlens_small_ratio)
print(tlens_medium_ratio)
print(tlens_large_ratio)
with open(os.path.join(savedir, "ratios-{}.txt".format(sid)), "w") as file:
    file.writelines(["tlens_small_ratio: {}".format(tlens_small_ratio),
                     "tlens_medium_ratio: {}".format(tlens_medium_ratio),
                     "tlens_large_ratio: {}".format(tlens_large_ratio),
                     ])
    


def plot_fig(nums, xlab, ylab, suptitle, title, save):
    plt.figure()
    plt.hist(nums, 100)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.suptitle(suptitle)
    plt.savefig(os.path.join(savedir, save))
    plt.close()

plot_fig(tlens, "Length", "Count", "Template lengths - " + sid, None, "tlen-{}.png".format(sid))
plot_fig(tlens_80_100, "Length", "Count", "Template lengths (80-100) - " + sid, tlens_small_ratio, "tlen-{}-80-100.png".format(sid))
plot_fig(tlens_100_500, "Length", "Count", "Template lengths (100-500) - " + sid, tlens_small_ratio, "tlen-{}-100-500.png".format(sid))
plot_fig(tlens_small, "Length", "Count", "Template lengths (small) - " + sid, tlens_small_ratio, "tlen-{}-small.png".format(sid))
plot_fig(tlens_medium, "Length", "Count", "Template lengths (medium) - " + sid, tlens_medium_ratio,  "tlen-{}-medium.png".format(sid))
plot_fig(tlens_large, "Length", "Count", "Template lengths (large) - " + sid, tlens_large_ratio, "tlen-{}-large.png".format(sid))

# TODO: process zerotlens

# categories:
# - can keep only positive, read1 and read2 have the same absolute value
# - non-zero/zero
# - 0-500, 500-3000, >3000
# - also report proportions of subset vs total