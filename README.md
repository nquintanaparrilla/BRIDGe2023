# BRIDGe REU 2023

## Project: Categorizing splicing ambiguity status RNA-sequencing reads in single cell and single nucleus data
Mentor: Dr. Rob Patro

### Purpose of Research
RNA-sequencing(RNA-seq) technologies like 10x Genomics provide the necessary tools and resources for researchers to use freely. Like other RNA-seq technologies, it follows a standard pipeline that goes from sample preparation, where we would obtain the raw data in a fastq file, to the quantification and finally the analysis. However, these tools don’t provide an intermediary between the quantification and the analysis. There isn’t any way to look at the data directly and estimate the splicing status of all of the reads. For this reason, millions of reads are not being accounted for nor do they possess the right splicing status because of splicing ambiguity. Splicing ambiguity occurs because of ambiguous reads. These reads are classified as ambiguous because we do not know whether they are spliced or unspliced. Therefore, it is important to try to categorize these reads to build a preprocessing tool that can better estimate the splicing status of all of the reads that are not being accounted for and provide accurate data.

### Methods

STAR aligner was used to create a genome index with ENCODE’s latest basic primary assembly fasta sequence and basic primary genome annotation of the common mouse and map the reads with very few parameters. 
Datasets downloaded from 10x Genomics Resources:
| Tissue Samples | Link |
| ----------- | ----------- |
| 1k Heart Cells from an E18 mouse	| https://www.10xgenomics.com/resources/datasets/10-k-heart-cells-from-an-e-18-mouse-v-3-chemistry-3-standard-3-0-0 |
| 5k Adult Mouse Heart Nuclei Isolated with Chromium Nuclei Isolation Kit	| https://www.10xgenomics.com/resources/datasets/5k-adult-mouse-heart-nuclei-isolated-with-chromium-nuclei-isolation-kit-3-1-standard |
| 5k Mouse E18 Combined Cortex, Hippocampus and Subventricular Zone Cells	| https://www.10xgenomics.com/resources/datasets/5-k-mouse-e-18-combined-cortex-hippocampus-and-subventricular-zone-cells-3-1-standard-6-0-0 |
| 5k Adult Mouse Brain Nuclei Isolated with Chromium Nuclei Isolation Kit	| https://www.10xgenomics.com/resources/datasets/5k-adult-mouse-brain-nuclei-isolated-with-chromium-nuclei-isolation-kit-3-1-standard |

The alignment process is important to procure a BAM file, and with samtools, create a BAI file (bam file index).


### Results

### Future work
