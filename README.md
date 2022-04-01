# SLANG

SLANG (Simple Long-read loci Assembly of Nanopore data for Genotyping) is a Python script pipeline for locus assembly, orthology estimation and SNP extraction of Nanopore-sequenced multi-locus data. 

There is no installation required for SLANG, only Python needs to be installed and either the Notebook version (.ipynb) or the Python script version (.py) needs to be downloaded and executed. SLANG was written with Python 3.7.6, but more recent versions work as well.

SLANG is separated into three steps:
1. Within-samples clustering (locus assembly): Input reads are clustered with vsearch according to read similarity. Similar reads are assumed to have the same origin and therefore belong to the same locus. Clustered reads are then mapped to their consensus sequence for further removal of highly erroneous reads, potential paralogs and reads of other loci (undersplitted loci).
2. Among-samples clustering (orthology inference): Orthology is then inferred by similarity clustering with vsearch of the passing consensus sequences from within-samples clustering. If multiple consensus sequences of one sample cluster together (potential paralogs or oversplitted loci), they will be removed. Only samples will be eligible for SNP calling, if only single consensus sequences cluster with those of other samples.
3. Reference-based SNP calling: SNPs are then called by mapping reads to the consensus sequence made during among-samples clustering.

The resulting SNPs are outputted in the VCF format and can be used for downstream phylogenetic inference, population genetic or phylogeographical studies.

![SLANG](https://user-images.githubusercontent.com/94844710/161264733-1c37a2c5-4bfe-4893-ba3a-68b2b108c625.png)

Here, a Python Notebook (.ipynb) and a Python script (.py) version of SLANG are offered, depending on usage preference. For running the notebook version, [Jupyter](https://jupyter.org/) is required.

### Dependencies:
- [VSEARCH](https://github.com/torognes/vsearch)
- [samtools/BCFtools/HTSlib](http://www.htslib.org/)
- [minimap2](https://github.com/lh3/minimap2)
- [vcflib](https://github.com/vcflib/vcflib)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)  
- [matplotlib](https://matplotlib.org/)
- [natsort](https://pypi.org/project/natsort/)

Installation with the package manager [conda](https://docs.conda.io/en/latest/) is recommended.

### Input and Usage:

SLANG assembles loci, infers orthology and calls SNPs only from your preprocessed (demultiplexed and quality -and length-filtered) Nanopore reads in FASTQ format. Before running SLANG, one FASTQ file of each sample is required. The name of the FASTQ file must include the Nanopore Barcode used for this sample, e.g. `bc01_samplename.fastq`. All FASTQ files should then be placed into the same directory.

Afterwards, the run can be configured in the SLANG script by opening it with a text editor. In the header of the SLANG script, input and parameters for the run can be set:

| parameter | explanation | example                |
| ------------- |:-------------:|:----------------------:|
| analysis_dir  | Defines a directory for your SLANG analysis. The directory will be made and if it already exists, you will be notified. Needs to be given in `''` and end with `/`. |'/home/mySLANGanalysis/'|
| filtered_reads_dir | Input path to your directory containing the preprocessed FASTQ reads. Needs to be given in `''` and end with `/`. | '/home/myFASTQreads/'
| barcodes | Insert all barcodes used in the analysis. Needs to be given in `''` and as a Python list in `[]`, separated by `,` | '[bc01, bc02, bc03]'
| within_samples_clustering_ct | Defines the similarity threshold for the within-samples clustering (locus assembly). Needs to be given as a value between 0.00 (0% similarity) and 1.00 (100% similarity) in `''`. | '0.75' |
| in_between_samples_clustering_ct | Defines the similarity threshold for the among-samples clustering (orthology inference). Needs to be given as a value between 0.00 (0% similarity) and 1.00 (100% similarity) in `''`. This value should be higher than the within-samples clustering similarity threshold. | '0.90'
| minimum_depth | Defines the minimum read depth for the analysis. Variant calling requires at least five reads per allele. | 10 |
| minimum_sample_per_locus | Defines how many samples need to share the same locus for the locus to pass filters and be analyzed. The value must be at least 2. | 2 |
| threads | Defines how many threads of your system should be used for certain tasks of SLANG, e.g. vsearch clustering and minimap2 mapping. The higher, the faster the analysis. Needs to be given in `''` | '10' |

With these parameters edited, SLANG can be executed through a command line interface by entering `python SLANG.py` and executing, if using the Python script version. 
When the Python Notebook version is used, execute the column of the parameters, check the printed output for errors, then execute the column below with the main script.
