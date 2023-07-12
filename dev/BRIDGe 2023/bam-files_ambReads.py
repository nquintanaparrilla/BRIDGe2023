import pandas
import pysam
import os
import numpy

save = pysam.set_verbosity(0)
bf = pysam.AlignmentFile("/Users/nataliaquintana/mappingAligned.sortedByCoord.out.bam", "rb")
if not os.path.exists("//Users/nataliaquintana/mappingAligned.sortedByCoord.out.bam.bai"):
    pysam.index("/Users/nataliaquintana/mappingAligned.sortedByCoord.out.bam")
pysam.set_verbosity(save)

print (bf.references)
print (bf.lengths)

