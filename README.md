# SLANG

SLANG (Simple Long-read loci Assembly of Nanopore data for Genotyping) is a Python script pipeline for locus assembly, orthology estimation and SNP extraction of Nanopore-sequenced multi-locus data. 

There is no installation required for SLANG, only Python needs to be installed and either the Notebook version (.ipynb) or the Python script version (.py) needs to be downloaded and executed. SLANG was written with Python 3.7.6, but more recent versions work as well.

SLANG is separated into three steps:
1. Within-samples clustering (locus assembly): Input reads are clustered with vsearch according to read similarity. Similar reads are assumed to have the same origin and therefore belong to the same locus. Clustered reads are then mapped to their consensus sequence for further removal of highly erroneous reads, potential paralogs and reads of other loci (undersplitted loci).
2. Among-samples clustering (orthology inference): Consensus sequences of the assembled loci are then clustered according to their similarity among all samples. Orthologous loci are therefore grouped together. Multiple consensus sequences of a single sample within a cluster (paralogous or oversplitted loci) are filtered. Clusters not meeting the minimum samples per locus count are also filtered.
3. Reference-based SNP calling: SNPs are then called by mapping reads sample per sample to the consensus sequence made during among-samples clustering. 

The resulting SNPs are outputted in the VCF format and can be used for downstream phylogenetic inference, population genetic or phylogeographical studies. In order to get a SNP matrix, the CLI tool 'SLANG_vcf_to_SNPs.py' can be used.

![SLANG](https://user-images.githubusercontent.com/94844710/161264733-1c37a2c5-4bfe-4893-ba3a-68b2b108c625.png)

Here, a Python Notebook (.ipynb) and a Python script (.py) version of SLANG are offered, depending on usage preference. For running the notebook version, [Jupyter](https://jupyter.org/) is required.

### Dependencies:
- [VSEARCH](https://github.com/torognes/vsearch): 
- [samtools/BCFtools/HTSlib](http://www.htslib.org/)
- [minimap2](https://github.com/lh3/minimap2)
- [vcflib](https://github.com/vcflib/vcflib)
- [numpy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)  
- [matplotlib](https://matplotlib.org/)
- [natsort](https://pypi.org/project/natsort/)
- [nei_vcf](https://github.com/TankredO/nei_vcf) (if Nei-Li genetic distances need to be calculated from the VCF output of SLANG)

Installation with the package manager [conda](https://docs.conda.io/en/latest/) is recommended.
Command lines for installing the dependencies:
- `conda install -c bioconda vsearch`
- `conda install -c bioconda samtools`
- `conda install -c bioconda bcftools`
- `conda install -c bioconda htslib`
- `conda install -c bioconda minimap2`
- `conda install -c bioconda vcflib`
- `conda install -c anaconda numpy`
- `conda install -c anaconda pandas`
- `conda install -c conda-forge matplotlib`
- `conda install -c anaconda natsort`

### Input and Usage:

SLANG assembles loci, infers orthology and calls SNPs only from your preprocessed (demultiplexed and quality -and length-filtered) Nanopore reads in FASTQ format. Before running SLANG, one FASTQ file of each sample is required. The name of the FASTQ file must include the Barcode used for this sample, e.g. `bc01_samplename.fastq`. All FASTQ files should then be placed into the same directory.

Afterwards, the run can be configured in the SLANG script by opening it with a text editor. In the header of the SLANG script, input and parameters for the run can be set:

| parameter | explanation | example   |
| ------------- |:-------------:|:----------------------:|
| analysis_dir  | Defines a directory for your SLANG analysis. The directory will be made and if it already exists, you will be notified. Needs to be given in `''` and end with `/`. |'/home/mySLANGanalysis/'|
| filtered_reads_dir | Input path to your directory containing the preprocessed FASTQ reads. Needs to be given in `''` and end with `/`. | '/home/myFASTQreads/'
| barcodes | Insert all barcode identifiers used in the analysis. Needs to be given in `''` and as a Python list in `[]`, separated by `,` | '[bc01, bc02, bc03]'
| within_samples_clustering_ct | Defines the similarity threshold for the within-samples clustering (locus assembly). Needs to be given as a value between 0.00 (0% similarity) and 1.00 (100% similarity) in `''`. | '0.75' |
| in_between_samples_clustering_ct | Defines the similarity threshold for the among-samples clustering (orthology inference). Needs to be given as a value between 0.00 (0% similarity) and 1.00 (100% similarity) in `''`. This value should be higher than the within-samples clustering similarity threshold. | '0.90'
| minimum_depth | Defines the minimum read depth for the analysis. Variant calling requires at least five reads per allele. | 10 |
| minimum_sample_per_locus | Defines how many samples need to share the same locus for the locus to pass filters and be analyzed. The value must be at least 2. | 2 |
| threads | Defines how many threads of your system should be used for certain tasks of SLANG, e.g. vsearch clustering and minimap2 mapping. The higher, the faster the analysis. Needs to be given in `''` | '10' |

With these parameters edited, SLANG can be executed through a command line interface by entering `python SLANG.py` and executing, if using the Python script version. 
When the Python Notebook version is used, execute the column of the parameters, check the printed output for errors, then execute the column below with the main script.

### Usage example: Adjusting within-samples and among-samples clustering similarity thresholds
An important part of the analysis is the choice of the within-samples and among-samples clustering similarity threshold. With these, loci and orthology are defined, therefore a careful selection is necessary in order to attain robust results. If similarity thresholds are too high, there is an increased risk of oversplitting the alleles of a locus into multiple clusters, and if similarity thesholds are too low, the risk of multiple loci clustered together (locus undersplitting) is increased. The goal is to find a similarity threshold, where both under- and oversplitting is minimized. While being optional and choosing the default values (0.75 and 0.90) or other non-extreme values can provide satisfactory results, it is recommended to adjust the similarity thresholds for each individual dataset. In the following, one possible way of determining similarity thresholds is explained.

To determine a valid within-samples clustering threshold, run multiple instances of the script with the within-samples similarity threshold changed in incremental steps (e.g. 0.05). Run the analysis until the `my_analysis/within_samples_clustering/unmapped_reads_log.csv` is produced (a notification is printed to the terminal), then count the number of clusters with any number of unmapped reads from it. An unmapped read could indicate a highly erroneous read, a paralogous read or a read belonging to another locus (locus undersplitting). For the reason of minimizing locus undersplitting, the number of clusters with unmapped reads should be low. Next, add up the total cluster count of all samples. The total cluster count is outputted by the script for each sample after within-samples clustering is finished. Locus oversplitting should be increased with higher total cluster count. Finally, multiply the total cluster count and the number of clusters with unmapped reads. This formula favors similarity thresholds, where high number of total clusters are formed, but with as less clusters with unmapped reads as possible. For this reason, the similarity threshold with the highest value should represent a tradeoff between locus over- and undersplitting.

After determining a within-samples clustering threshold, again run multiple instances of the script with the previously determined within-samples clustering similarity threshold, but this time with the among-samples similarity threshold changed in incremental steps (e.g. 0.05). The value should be higher than the within-samples similarity threshold. Here, only clusters will pass if each sample is represented by a single consensus sequence. High similarity thresholds lead to many single sample clusters, which are filtered. Low similarity thresholds will result in clusters, where multiple consensus sequences of the same sample are clustered together, which are also filtered out (indication of oversplitted loci clustering back together). As a consequence, an optimal among-samples clustering similarity threshold filters out as few clusters as possible. The VCF file resulting from this analysis can then be used for downstream analysis.
