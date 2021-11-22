# SLANG v1.0

# parameters
# please fill out manually

analysis_dir = 'enter/path/here/with/dir/name' # dir will be created
filtered_reads_dir = '/enter/path/to/preprocessed/fastq/reads' # dir need to end with '/'
barcodes = ['bc01', 'bc02', 'bc03', ect.] # need to be given as a python list
within_samples_clustering_ct = '0.75' # needs to be lower than the among-samples clustering threshold
in_between_samples_clustering_ct = '0.90' # needs to be higher than the within-samples clustering threshold
minimum_depth = 10
minimum_sample_per_locus = 2 # needs to be at least 2
threads = '1' # needs to be given in ''

####################################################

print('Analysis will be done in:', analysis_dir)
print('Filtered reads will be taken from:', filtered_reads_dir)
print('There are', len(barcodes), 'samples. Barcodes are:', barcodes)
print('Clustering within-samples will have a clustering threshold of', within_samples_clustering_ct)
print('Clustering among-samples will have a clustering threshold of', in_between_samples_clustering_ct)
print('Depth will be at a minimum of', minimum_depth)
print('There need to be at least', minimum_sample_per_locus, 'samples to accept a locus')
print('Clustering will use', threads, 'thread(s)')

import numpy as np
import matplotlib.pyplot as plt
from natsort import natsorted, ns
import glob
import os
import pandas as pd
import subprocess
from itertools import islice
import shutil
import sys
from datetime import datetime
import time
import ntpath

start_time = time.time()


def path_leaf(path):
    '''Returns only the filename as a string if a path is given'''
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

def get_infos_from_cluster_file(cluster_file_path):
    '''Extracts barcode number and read id from vsearch cluster file'''
    infos = []
    with open(cluster_file_path) as cf:
        for line in cf:
            if line.startswith('>'):
                entry_name, entry_id, entry_nseq = line.strip().split('=')
                entry_name = entry_name[1:].replace('_centroid', '')
                entry_id = entry_id.split(';')[0]                
                infos.append([entry_name, entry_id])    
        return infos

def write_infos_to_file(infos, info_file, mode='x'):
    '''Writes barcode number and read id from vsearch cluster file, extracted by get_infos_from_cluster_file(), to output file'''
    with open(info_file, mode) as f:
        f.write('barcode\tread_id\n')
        for entry in infos:
            f.write(f'{entry[0]}\t{entry[1]}\n')
            
def take_four(it):
    '''While there are elements in input, take four elements'''
    fastq_entries = []
    
    _it = iter(it)
    while True:
        try:
            yield next(_it), next(_it), next(_it), next(_it)
        except StopIteration:
            # no more elements in the iterator
            return

def consout_to_dict(infile):
    '''Parses vsearch consout fasta output and returns it as a dictionary'''
    with open(infile) as cons_file:
        test_dict = {}
        
        cons_file = cons_file.read().split('>')                                    # splits the fasta file by '>'
        while '' in cons_file:                                                     # removes empty list elements
            cons_file.remove('')
        cons_file = ['>' + i for i in cons_file]                                   # restores fasta header by readding '>'
        
        for i in cons_file:                                                        # loop through consout fasta entries
            if i.startswith('>'):
                centroid_read_id = i[i.find('_centroid=')+10:i.find(';seqs=')]     # extracts centroid read id
                sequence = i[i.find('\n')+1:]                                      # extracts the DNA sequence
                test_dict[centroid_read_id] = sequence                             # adds fasta entry to dictionary
                                                                                   # dictionary key = centroid read id
        return test_dict                                                           # dictionary value = DNA sequence

    
# create directory for full analysis

if os.path.exists(analysis_dir):                                                             # check if directory already exists
    userinput = input('Chosen directory already exists. Do you wish to overwrite? (y/n)')    # directory exists, ask user to overwrite
    print('\n')
    if userinput == 'n':
        print('Directory was not overwritten')
    if userinput == 'y':                                                                     # user chose yes
        shutil.rmtree(analysis_dir)                                                          # remove directory with content
        if not os.path.exists(analysis_dir):                                                 # if directory now does not exist:
            os.makedirs(analysis_dir)                                                        # make the directory
            print(analysis_dir, 'was created')
else:
    if not os.path.exists(analysis_dir):                                                     # if directory does not exist:
        os.makedirs(analysis_dir)                                                            # make the directory
        print(analysis_dir, 'was created')

        
# create directory for within-sample-clustering

within_sampe_clustering_dir = analysis_dir + '/within_sample_clustering'    # name and path within-sample-clustering directory
if not os.path.exists(within_sampe_clustering_dir):
    os.makedirs(within_sampe_clustering_dir)                               # make the directory
    print(within_sampe_clustering_dir, 'was created')
else:
    print(within_sampe_clustering_dir, 'seems to already exist')


# create paths for all defined barcodes

for i in range(len(barcodes)):
    wsc_sample_dir = within_sampe_clustering_dir + '/' + barcodes[i]
    os.makedirs(wsc_sample_dir)
    print(wsc_sample_dir, 'was created')

print('\n')


# vsearch within samples clustering

loopcounter = 0

samplename_args = filtered_reads_dir + '/*.fastq'
samplename = glob.glob(samplename_args)             # list all fastqs in input read directory


for i in range(len(barcodes)):                                                            # loop through each sample/barcode
    
    wsc_sample_dir = within_sampe_clustering_dir + '/' + barcodes[i]                      # define working directory for each sample/barcode
    
    msaout = wsc_sample_dir + '/msaout_' + barcodes[i]                                    # define the vsearch msaout path and name
    consout = wsc_sample_dir + '/consout_' + barcodes[i]                                  # define the vsearch consout path and name
    clusters_dir = wsc_sample_dir + '/clusters'                                           # define the vsearch clusters path and directory name
    
    if not os.path.exists(clusters_dir):                                                  # if the clusters directory does not yet exist
        os.makedirs(clusters_dir)                                                         # create it
    
    clusters = clusters_dir + '/' + barcodes[i] + '_cluster_'                             # define the vsearch clusters names
    stdout_path = wsc_sample_dir + '/vsearch_stdout'                                      # define the vsearch stdout
    stderr_path = wsc_sample_dir + '/vsearch_stderr'                                      # define the vsearch stderr
    
    for j in samplename:
        if barcodes[i] in path_leaf(j):
            vsearch_input = j

    #vsearch_input = filtered_reads_dir + barcodes[i] + '_q7_200to1000bp.fastq'            # define the path to the correct sample.fastq
    
    vsearch_wsc_args = ['vsearch', 
                        '--cluster_fast', vsearch_input, 
                        '--id', within_samples_clustering_ct, 
                        '--msaout', msaout, 
                        '--consout', consout, 
                        '--clusters', clusters, 
                        '--threads', threads]                                             # list of the vsearch command line input formatted for the subprocess package
    
    if os.path.exists(stdout_path):                                                       # remove vsearch stdout and stderr if it already exists from previous runs
        os.remove(stdout_path)
    if os.path.exists(stderr_path):
        os.remove(stderr_path)
    
    loopcounter = loopcounter + 1  
    print('Sample', loopcounter, 'of', len(barcodes), ': Clustering', path_leaf(vsearch_input), '(', barcodes[i], ') at vsearch --id', within_samples_clustering_ct, '...', end='\r')
   
    with open(stdout_path, 'x') as stdout_file:                                                                   # create and open files to write stderr and stdout to
        with open(stderr_path, 'x') as stderr_file:
            run_vsearch_wsc = subprocess.run(vsearch_wsc_args, stdout = stdout_file, stderr = stderr_file)        # run vsearch

            if run_vsearch_wsc.returncode == 1:                                                                   # print a warning message if vsearch process was unsuccessful
                print('Warning. Something went wrong concerning the within-sample-clustering of', barcodes[i])
                print(run_vsearch_wsc, '\n')
                break

print('\n')

# count reads in every cluster = clustercount

loopcounter = 0

for i in range(len(barcodes)):                                                       # loop through every sample
    loopcounter = loopcounter + 1  
    
    wsc_sample_dir = within_sampe_clustering_dir + '/' + barcodes[i]                 # define working directory for the current sample/barcode
    
    clustercount_dir = wsc_sample_dir + '/clustercount_' + barcodes[i] + '.csv'      # define path and name for the clustercount.csv

    if os.path.exists(clustercount_dir):                                             # if the clustercount.csv already exists (from a previous run)
        os.remove(clustercount_dir)                                                  # remove it
    if not os.path.exists(clustercount_dir):                                         # if it does not exist yet
        with open(clustercount_dir, 'x') as f:                                       # create and open the clustercount.csv
            f.write('cluster\tcount\n')                                              # write a header with 2 columns: cluster(name) and (read)count

    clusters_in_dir = wsc_sample_dir + '/clusters' + '/*'                            # define the vsearch clusters directory
    list_clusters = glob.glob(clusters_in_dir)                                       # list those cluster fasta files
    
    # inform about the current progress of the script
    print('Sample', loopcounter, 'out of', len(barcodes), ': Counting reads per cluster of', barcodes[i], '(', len(list_clusters), 'clusters )', end='\r')
    
    for j in list_clusters:                                # loop through each vsearch clusters file
        with open(j) as f:                                 # open a cluster fasta file from the list
            n = 0                                          # start a counter at zero
            clustercount = []                              # placeholder variable for the (read)count of the current cluster
            for line in f:                                 # loop through the cluster fasta file line per line
                if line[0] == '>':                         # if it is a header line starting with '>'
                    n = n + 1                              # count + 1
            clustercount.append(n)                         # add the final (read)count number to the placeholder
            
            with open(clustercount_dir, 'a') as f:         # open the clustercount.csv file
                f.write(path_leaf(j))                      # write in the 'cluster' column the clustername
                f.write('\t')
                f.write(str(clustercount[0]))              # write in the 'count' column the (read)count number
                f.write('\n')
    
    clustercount_df = pd.read_csv(clustercount_dir, sep='\t')                        # open the finished clustercount.csv in pandas
    clustercount_hist = clustercount_dir.replace('.csv', '_hist.csv')                # name and path the clustercount histogram output
    clustercount_df.groupby('count').count().to_csv(clustercount_hist, sep="\t")     # write a histogram file counting how many clusters have how many reads

print('\n')    


# prepare a list of just the normal clustercount.csv without the _hist versions

path_to_clustercount_csvs = within_sampe_clustering_dir + '/*/clustercount*.csv'     # gather all clustercount files (also lists the hist versions)
list_of_clustercount_files = glob.glob(path_to_clustercount_csvs)

remove_hist_files = []

for i in list_of_clustercount_files:
    if 'hist' in path_leaf(i):
        remove_hist_files.append(i)

list_of_clustercount_files = [x for x in list_of_clustercount_files if x not in remove_hist_files]

        
# loop through all (non-histogram (removed above)) clustercount files for statistics to print

for i in list_of_clustercount_files:
    print(i)
    df = pd.read_csv(i, sep='\t')
    print('Total clusters:\t\t\t\t', len(df))
    print('Total reads:\t\t\t\t', df['count'].sum())
    print('Mean reads in cluster:\t\t\t', round(df['count'].mean(), 2), '+-', round(df['count'].std(), 2))
    print('Median reads in cluster:\t\t', round(df['count'].median(), 2))
    print('Singleton clusters:\t\t\t', df['count'].isin([1]).sum(), 
          '(', (round(df['count'].isin([1]).sum() / len(df), 4) * 100), '% )')
    print('Clusters with', minimum_depth, 'or more reads:\t\t', df[df['count'] >= minimum_depth].count()[0],
          '(', (round(df[df['count'] >= minimum_depth].count()[0] / len(df), 4) * 100), '% )')
    print('\n')


# add an identifer to the consensus fasta sequences in consout vsearch output

list_of_consout = []
for i in range(len(barcodes)):
    wsc_sample_dir = within_sampe_clustering_dir + '/' + barcodes[i]
    consout = wsc_sample_dir + '/consout_' + barcodes[i]
    list_of_consout.append(consout)

loopcounter = 0

for i in list_of_consout:                                                 # loop through the consout files
    loopcounter = loopcounter + 1
    print('Sample', loopcounter, 'out of', len(list_of_consout), ': Writing the barcode in the consout fasta headers for', i, end='\r')
    for j in barcodes:                                                    # loop through barcodes
        if j in path_leaf(i):                                             # if a barcode matches the barcode in the consout file name:
            with open(i) as f:                                            # open the respective consout file
                outfile_path = i + '_with_barcode.fasta'                  # name the output with barcode consout file after the respective consout
                if os.path.exists(outfile_path):
                    os.remove(outfile_path)                               # remove any output from previous runs
                if not os.path.exists(outfile_path):
                    with open(outfile_path, 'x') as outfile:              # create and open the consout_filtered.fasta file
                        for line in f:
                            if line.startswith('>'):
                                replacer = '>' + j + '_'
                                line = line.replace('>', replacer)        # add the respective barcode to the fasta header behind the '>'
                            outfile.write(line)

print('\n')


# filter out low-depth clusters in the vsearch consout file   

consout_dirs = within_sampe_clustering_dir + '/*/consout*_with_barcode.fasta'
list_of_consout = glob.glob(consout_dirs)

loopcounter = 0

for i in list_of_consout:
    
    loopcounter = loopcounter + 1
    for bc in barcodes:
        if bc in path_leaf(i):
            current_sample = bc
    print('Sample', loopcounter, 'out of', len(list_of_consout), ': Filtering sample', current_sample, 'for clusters >=', minimum_depth, 'reads...')
    
    with open(i) as f:
        # parse the fasta
        content = f.read()                                    # loads fasta in readable form 
        content_split = content.split('>')                    # create a list of the fasta, splitting entries by '>'
        while '' in content_split:
            content_split.remove('')                          # remove emtpy entries in list
        content_split = [">" + j for j in content_split]      # restore the fasta format by readding '>' at the beginning of the header
        
        # write the filtered consout file
        outfile_path = i.replace('_with_barcode.fasta', '_filtered.fasta')              # name the output 'consout*_filtered.fasta'
        if os.path.exists(outfile_path):
            os.remove(outfile_path)                                                     # remove any output from previous runs 
        if not os.path.exists(outfile_path):
            with open(outfile_path, 'x') as outfile:                                    # create and open 'consout*_filtered.fasta' file
                for j in content_split:
                    if int(j[j.find(';seqs=')+6 : j.find('\n')]) >= minimum_depth:      # if seqs= is higher than or exactly the defined minimum depth:
                        outfile.write(j)                                                # write to 'consout*_filtered.fasta'

            # how many clusters are left?
            with open(outfile_path) as outfile:      # reopen the newly written consout_*_filtered.fasta
                n = 0
                conscount = []
                for line in outfile:
                    if line[0] == '>':
                        n = n + 1                    # count + 1 if there is a '>' = one fasta line
                conscount.append(n)                  # append the number to the conscount variable
                
                # print details to within-sample clustering depth filtering
                print('Filtered out', (len(content_split) - conscount[0]), 'low-depth clusters', 'out of', len(content_split), '.', conscount[0], 'clusters remain for', current_sample)              
    
    os.remove(i) # remove the '_with_barcode.fasta' file as it is now redundant
                            
print('\n')


# Move clusters with 'depth' >= 'minimum_depth' to 'depth_filtered_clusters_dir'

loopcounter = 0

for bc in barcodes:
    unfiltered_clusters_args = within_sampe_clustering_dir + '/' + bc + '/clusters/*'
    unfiltered_clusters = glob.glob(unfiltered_clusters_args)
    depth_filtered_clusters_dir = within_sampe_clustering_dir + '/' + bc + '/depth_filtered_clusters'
    
    loopcounter = loopcounter + 1
    print(loopcounter, 'out of', len(barcodes), ': Moving', bc, 'clusters with >=', minimum_depth, 'depth to', depth_filtered_clusters_dir, end='\r')
    
    # create directory for depth-filtered clusters
    
    if os.path.exists(depth_filtered_clusters_dir):
        shutil.rmtree(depth_filtered_clusters_dir)
    if not os.path.exists(depth_filtered_clusters_dir):
        os.makedirs(depth_filtered_clusters_dir)
    
    
    # Move clusters with 'depth' >= 'minimum_depth' to 'depth_filtered_clusters_dir'

    for i in unfiltered_clusters:
        with open(i) as f:
            depth = 0
            for line in f.read():
                if line.startswith('>'):
                    depth = depth + 1
            if depth >= minimum_depth:
                shutil.move(i, depth_filtered_clusters_dir) 

print('\n')

          
# map reads to consensus sequence
          
list_of_consout_args = within_sampe_clustering_dir + '/*/consout*filtered.fasta'
list_of_consout = glob.glob(list_of_consout_args)


loopcounter = 0

for i in list_of_consout:
    # select the corresponding /within_samples_clustering/*/depth_filtered_clusters directory
    
    for bc in barcodes:
        if bc in path_leaf(i):
            clusters_dir_args = within_sampe_clustering_dir + '/' + bc + '/depth_filtered_clusters/*'
            clusters_dir = glob.glob(clusters_dir_args)
            
            consout_ref_files = within_sampe_clustering_dir + '/' + bc + '/clusters_reads_to_consout_mapping'
            
            if os.path.exists(consout_ref_files):
                shutil.rmtree(consout_ref_files)
            if not os.path.exists(consout_ref_files):
                os.makedirs(consout_ref_files)
    
            
    # select centroid ids and search for them in the clusters directory to identify corresponding reads
            
    with open(i) as f:
        loopcounter = loopcounter + 1

        # parse consout.fasta
        content = f.read().split('>')
        while '' in content:
            content.remove('')
        content = ['>' + j for j in content]
        
        loopcounter2 = 0
        
        for k in content:
            #print(k)
            
            loopcounter2 = loopcounter2 + 1
            print(loopcounter2, 'out of', len(content), 
                  '- Currently mapping reads to their correspondong consout.fasta entry:', path_leaf(i), 
                  '(', loopcounter, 'out of', len(list_of_consout), 'consouts )', end='\r')
            
            if k.startswith('>'):
                header = k[k.find('_centroid=')+10 : k.find(';seqs=')]         # extract centroid read id from consout fasta entry
                #print(header)
                
                consout_ref = consout_ref_files + '/' + header                 # name reference file
                #print(consout_ref)
                
                with open(consout_ref, 'x') as consout_ref_f:
                    consout_ref_f.write(k)                                     # write reference to file
                
                
                # find corresponding 'within_sample_clustering/*/clusters/cluster_*'
                
                for j in clusters_dir:
                    with open(j) as cluster:
                        if header in cluster.read():
                            #print('Found header in', j)
                    
                            
                            # mapping to reference
                            
                            stdout_path = consout_ref_files + '/' + path_leaf(j) + '.sam'        # mapping output .sam file
                            if os.path.exists(stdout_path):
                                os.remove(stdout_path)
                            with open(stdout_path, 'x') as stdout_outfile:
                                mapping = subprocess.run(['minimap2', 
                                                          '-ax', 'map-ont', 
                                                          '-t', threads, 
                                                          consout_ref,                 # reference path
                                                          j],                          # input fastq path
                                                          stdout = stdout_outfile      # minimap2 .sam output
                                                        )
                            if mapping.returncode == 1:
                                print('Error. Subprocess for minimap2 mapping returncode = 1. Please check your output.')
                                break
         
            else:
                print('Error.', k, 'is not in correct FASTA format.')
                break


# log mapping data (clusters with unmapped reads)
                
unmapped_reads_log_dir_wsc = within_sampe_clustering_dir + '/unmapped_reads_log.csv'

if os.path.exists(unmapped_reads_log_dir_wsc):                                                          # check if missing reads in the mapping file exists
    os.remove(unmapped_reads_log_dir_wsc)                                                               # remove it
if not os.path.exists(unmapped_reads_log_dir_wsc):                                                      # if file does not exist
    with open(unmapped_reads_log_dir_wsc, 'x') as f:                                                    # make it
        f.write('barcode\tcluster\tunmapped reads\ttotal reads in cluster\tPercentage unmapped reads\n')     # and write the header
    print('Unmapped reads are logged to:', unmapped_reads_log_dir_wsc)

print('\n')

list_of_all_sam_args = within_sampe_clustering_dir + '/*/clusters_reads_to_consout_mapping/*.sam'
list_of_all_sam = glob.glob(list_of_all_sam_args)                                                      # list the sam files made from minimap2 mapping above

loopcounter = 0

for i in list_of_all_sam:
    loopcounter = loopcounter + 1
    print('Now counting unmapped reads in all .sam files...', loopcounter, 'out of', len(list_of_all_sam), end='\r')
    
    for bc in barcodes:
        if bc in i:
            barcode = bc
    
    with open(i) as f:
        unmapped_reads_count = 0
        total_reads = []
        for line in f:
            if not line.startswith('@'):
                line = line.split('\t')
                if line[0] not in total_reads:
                    total_reads.append(line[0])
                if line[2] == '*':
                    unmapped_reads_count = unmapped_reads_count + 1
        #print(unmapped_reads_count, 'out of', len(total_reads), 'reads are unmapped')
   
    with open(unmapped_reads_log_dir_wsc, 'a') as unmapped_reads_f:
        unmapped_reads_f.write(barcode)
        unmapped_reads_f.write('\t')
        unmapped_reads_f.write(path_leaf(i))
        unmapped_reads_f.write('\t')
        unmapped_reads_f.write(str(unmapped_reads_count))
        unmapped_reads_f.write('\t')
        unmapped_reads_f.write(str(len(total_reads)))
        unmapped_reads_f.write('\t')
        percentage_unmapped_reads = round( ( (unmapped_reads_count / len(total_reads) ) * 100), 2)
        unmapped_reads_f.write(str(percentage_unmapped_reads))
        unmapped_reads_f.write('\n')

print('\n')


# open the unmapped_reads_log.csv for statistics

df = pd.read_csv(unmapped_reads_log_dir_wsc, sep='\t')
df['mapped reads'] = df['total reads in cluster'] - df['unmapped reads'] # adding a column with mapped reads

print('depth-filtered clusters:\t\t\t', len(df))
print('clusters with unmapped reads:\t\t\t', len(df[df['unmapped reads'] > 0]),
      '(', round(len(df[df['unmapped reads'] > 0]) / len(df), 4) * 100, '% of depth-filtered clusters)',
      '( mean sample:', len(df[df['unmapped reads'] > 0]) / len(barcodes), ')')
print('clusters with less than', minimum_depth, 'unmapped reads:\t', len(df[df['unmapped reads'] < minimum_depth]),
      '(', round(len(df[df['unmapped reads'] < minimum_depth]) / len(df), 4) * 100, '% of depth-filtered clusters)')
print('clusters with', minimum_depth, 'or more unmapped reads:\t', len(df[df['unmapped reads'] >= minimum_depth]),
      '(', round((len(df[df['unmapped reads'] >= minimum_depth]) / len(df)) * 100, 2), '% of depth-filtered clusters)')
print('clusters with', minimum_depth, 'or more mapped reads:\t\t', len(df[df['mapped reads'] >= minimum_depth]),
      '(', round(len(df[df['mapped reads'] >= minimum_depth]) / len(df), 4) * 100, '% of depth-filtered clusters )')
print('\n')

# create a pandas dataframe with clustering and mapping data for each sample for the user to assess an optimal clustering threshold for within-samples-clustering

wsc_data = {}                       # empty placeholder dictionary
wsc_df = pd.DataFrame(wsc_data)     # create the dataframe by using an empty dictionary
wsc_df['samples'] = barcodes        # add samples as rows

# add 'total clusters' and 'clusters with >= minimum_depth' as a pandas dataframe column
total_clusters_column = {}       # placeholder dictionary to add entries to
clusters_mindepth_column = {}

for i in barcodes:                                                # loop through barcodes
    for j in list_of_clustercount_files:                          # loop through the clustercount files
        if i in path_leaf(j):                                     # if current barcode matches one of the clustercount files
            df = pd.read_csv(j, sep='\t')                         # load the clustercount file as a pandas dataframe
            tmp_dict_1 = {len(df) : i}                            # combine barcode with its total cluster value
            total_clusters_column.update(tmp_dict_1)              # add dictionary entry to placeholder dictionary
            
            tmp_df = df[df['count'] >= minimum_depth].count()[0]  # count "clusters with >= minimum_depth reads" for barcode i
            tmp_dict_2 = {tmp_df : i}                             # combine barcode i with its "clusters with >= minimum_depth reads" value
            clusters_mindepth_column.update(tmp_dict_2)           # add dictionary entry to placeholder dictionary

# add 'total clusters' and 'clusters with >= minimum_depth' as a pandas dataframe column
total_clusters_column = []       # placeholder dictionary to add entries to
clusters_mindepth_column = []

test_barcodes = ['bc09']

for i in barcodes:                                                # loop through barcodes
    for j in list_of_clustercount_files:                          # loop through the clustercount files
        if i in path_leaf(j):                                     # if current barcode matches one of the clustercount files
            df = pd.read_csv(j, sep='\t')                         # load the clustercount file as a pandas dataframe
            total_clusters_column.append(len(df))                 # count total clusters by length of the clustercount file
            clusters_mindepth_column.append(df[df['count'] >= minimum_depth].count()[0])  # count "clusters with >= minimum_depth reads" for barcode i
            
# add mapping data to pandas dataframe as a column
df = pd.read_csv(unmapped_reads_log_dir_wsc, sep='\t')                      # load the 'unmapped_reads_log_dir_wsc' as 'df'
df['mapped reads'] = df['total reads in cluster'] - df['unmapped reads']    # adding a column with mapped reads

depth_filtered_clusters_column = []         # empty placeholder dictionaries to add entries to
unmapped_reads_clusters_column = []
less_than_min_depth_unmapped_column = []
mincov_or_more_unmapped_column = []
mincov_or_more_mapped_column = []

for i in barcodes:
    df_i = df[df['barcode'] == i]
    depth_filtered_clusters_column.append(len(df_i))
    unmapped_reads_clusters_column.append(len(df_i[df_i['unmapped reads'] > 0]))
    less_than_min_depth_unmapped_column.append(len(df_i[df_i['unmapped reads'] < minimum_depth]))
    mincov_or_more_unmapped_column.append(len(df_i[df_i['unmapped reads'] >= minimum_depth]))
    mincov_or_more_mapped_column.append(len(df_i[df_i['mapped reads'] >= minimum_depth]))   

    
# add completed (placeholder) dictionaries as columns to pandas dataframe
wsc_df['total clusters'] = total_clusters_column

column_name_arg = 'clusters with >=' + str(minimum_depth) + ' reads'
wsc_df[column_name_arg] = clusters_mindepth_column

wsc_df['depth-filtered clusters'] = depth_filtered_clusters_column

wsc_df['clusters with unmapped reads'] = unmapped_reads_clusters_column

column_name_arg = 'clusters with <' + str(minimum_depth) + ' unmapped reads'
wsc_df[column_name_arg] = less_than_min_depth_unmapped_column

column_name_arg = 'clusters with >=' + str(minimum_depth) + ' unmapped reads'
wsc_df[column_name_arg] = mincov_or_more_unmapped_column

column_name_arg = 'clusters with >=' + str(minimum_depth) + ' mapped reads'
wsc_df[column_name_arg] = mincov_or_more_mapped_column


# save the pandas dataframe as a csv
wsc_optimization_file = within_sampe_clustering_dir + '/clustering_mapping_statistics.csv'

if os.path.exists(wsc_optimization_file):
    os.remove(wsc_optimization_file)
    
wsc_df.to_csv(wsc_optimization_file, sep='\t')
print('An in-depth .csv file was placed here:', wsc_optimization_file)
print('\n')

# remove consout consensus sequence entries with < 'minimum_depth' mapped reads

df = pd.read_csv(unmapped_reads_log_dir_wsc, sep='\t')                     # open unmapped reads csv in pandas
df['mapped reads'] = df['total reads in cluster'] - df['unmapped reads']   # adding a column with mapped reads
df_low_mapped_reads = df[df['mapped reads'] < minimum_depth]               # reducing the dataset to cluster with < 'minimum_depth' mapped reads

loopcounter = 0

for bc in barcodes:
    loopcounter = loopcounter + 1
    print('Sample', loopcounter, 'out of', len(barcodes), ': Filter out consensus sequences in', bc, 'consout with <', minimum_depth, 'mapped reads', end='\r')
    
    df_bc_low_mapped_reads = df_low_mapped_reads[df_low_mapped_reads['barcode'] == bc]    # get data for current barcode
    list_low_mapped_reads = df_bc_low_mapped_reads['cluster'].tolist()                    # list cluster with < 'minimum_depth' mapped reads for current barcode
    
    
    # list centroid read ids from sam files (will be used to filter out consensus sequences in vsearch consout with < 'minimum_depth' mapped reads)
    consout_rm_list = []          # placeholder list for appending centroid read ids of clusters with < 'minimum_depth' mapped reads to
    
    for i in list_low_mapped_reads:
        path_to_sam = within_sampe_clustering_dir + '/' + bc + '/clusters_reads_to_consout_mapping/' + i       # get full path to sam file
        with open(path_to_sam) as f:
            content = f.read()
            rm_centroid_read_id = content[content.find('_centroid=')+10 : content.find(';seqs=')]              # get centroid read id whose entry in consout is to be removed
            consout_rm_list.append(rm_centroid_read_id)
     
    
    # parse vsearch consout.fasta
    consout_path = within_sampe_clustering_dir + '/' + bc + '/consout_' + bc + '_filtered.fasta'
    with open(consout_path) as consout:
        consout = consout.read().split('>')
        while '' in consout:
            consout.remove('')
        consout = ['>' + l for l in consout]
    
    
    # remove consensus sequences with < 'minimum_depth' mapped reads from consout.fasta
    to_remove = []                                  # placeholder list with consensus sequences to remove
    
    for consensus in consout:
        for read_id in consout_rm_list:
            if read_id in consensus:
                    to_remove.append(consensus)     # append consensus sequence to remove to list above

    for consensus in to_remove:
        consout.remove(consensus)                   # safely remove consensus sequence from consout.fasta
        
        
    # overwrite the consout.fasta with the filtered consout consensus sequences
    if os.path.exists(consout_path):
        os.remove(consout_path)
    if not os.path.exists(consout_path):
        with open(consout_path, 'x') as f:
            for entry in consout:
                f.write(entry)
                
print('\n')

# remove unmapped reads from clusters.fasta in /depth_filtered_clusters/*

# load df to identify clusters with unmapped reads
df = pd.read_csv(unmapped_reads_log_dir_wsc, sep='\t')

loopcounter = 0

for bc in barcodes:
    loopcounter = loopcounter + 1
    
    # create (temporary) directories for clusters.fasta with only mapped reads
    mapped_reads_only_dir = within_sampe_clustering_dir + '/' + bc + '/mapped_reads_only_clusters'
    if os.path.exists(mapped_reads_only_dir):
        shutil.rmtree(mapped_reads_only_dir)
    if not os.path.exists(mapped_reads_only_dir):
        os.makedirs(mapped_reads_only_dir)
    
    
    # list the sam files with unmapped reads
    df_bc = df[df['barcode'] == bc]                             # gather entries of sample with bc = bc
    df_bc_unmapped = df_bc[df_bc['unmapped reads'] > 0]         # only include cluster with unmapped reads
    df_bc_unmapped_list = df_bc_unmapped['cluster'].tolist()    # list those clusters with unmapped reads
    for i, j in enumerate(df_bc_unmapped_list):
        j = within_sampe_clustering_dir + '/' + bc + '/clusters_reads_to_consout_mapping/' + j    # complete path of .sam
        df_bc_unmapped_list[i] = j                                                                # replace sam with pathed sam

    # list all depth-filtered-clusters    
    depth_filtered_clusters_args = within_sampe_clustering_dir + '/' + bc + '/depth_filtered_clusters/*_cluster_*'
    depth_filtered_clusters = glob.glob(depth_filtered_clusters_args)
    
    
    loopcounter2 = 0  
    
    # list all unmapped reads of cluster 'k'
    for k in df_bc_unmapped_list:                        # loop through the sam files
        unmapped_reads = []                              # placeholder list to append read ids of unmapped reads to
        with open(k) as f:                               # open sam file
            for line in f:                               # loop through sam file line per line
                if not line.startswith('@'):             # ignore header
                    line = line.split('\t')              # split lines by tab-del
                    if line[2] == '*':                   # gather unmapped entries
                        unmapped_reads.append(line[0])   # and append to list of unmapped reads of current cluster
            
            
        # navigate to cluster in depth_filtered_clusters, open cluster, list cluster, remove unmapped entries, write only mapped reads to new file in new directory
        k = k.replace('.sam', '')                        # remove the .sam from the element, so the file name = cluster name

        for cluster in depth_filtered_clusters:                            # loop through all depth-filtered clusters
            if path_leaf(cluster) == path_leaf(k):                         # find the matching cluster by their names
                
                loopcounter2 = loopcounter2 + 1
                
                marked_unmapped_reads = []
                with open(cluster) as f:                                   # open the cluster if file name matches

                    
                    # parse consout.fasta
                    parsed_fasta = f.read().split('>')
                    while '' in parsed_fasta:
                        parsed_fasta.remove('')
                    parsed_fasta = ['>' + l for l in parsed_fasta]

                for entry in parsed_fasta:                                 # append all unmapped reads to a list
                    for unmapped_read in unmapped_reads:
                        if unmapped_read in entry:
                            marked_unmapped_reads.append(entry)

                mapped_reads_only = list(set(parsed_fasta)^set(marked_unmapped_reads))            # remove unmapped reads from the parsed cluster.fasta

                mapped_reads_only_cluster = mapped_reads_only_dir + '/' + path_leaf(cluster)      # define path and name for new file with only mapped reads

                with open(mapped_reads_only_cluster, 'x') as f:                                   # write only the mapped reads into new file
                    for entry in mapped_reads_only:
                        f.write(entry)
                        
                        
                # move mapped reads only cluster to 'depth-filtered_clusters' and replace its old cluster with unmapped reads
                new_file_path = within_sampe_clustering_dir + '/' + bc + '/depth_filtered_clusters/' + path_leaf(mapped_reads_only_cluster)
                shutil.move(mapped_reads_only_cluster, new_file_path)

                print('Sample', loopcounter, 'out of', len(barcodes), ': Filter out unmapped reads from', bc, 'clusters:', 
                      loopcounter2, 'out of', len(df_bc_unmapped_list), 'clusters with unmapped reads', end='\r')
     
    
# delete 'mapped_reads_only_dir'
                                            
print('\n')

# count reads per cluster, remove clusters with less than 10 mapped reads

loopcounter = 0

for bc in barcodes:
    loopcounter = loopcounter + 1
    print('Sample', loopcounter, 'out of', len(barcodes), ': Remove', bc, 'clusters with less than', minimum_depth, 'mapped reads', end='\r')
    
    depth_filtered_clusters_args = within_sampe_clustering_dir + '/' + bc + '/' + 'depth_filtered_clusters/*_cluster_*'
    depth_filtered_clusters = glob.glob(depth_filtered_clusters_args)

    low_unmapped_count = []                           # placeholder list for paths to clusters with less than 'minimum_depth' mapped reads
    
    for i in depth_filtered_clusters:                 # loop through all clusters in '/depth_filtered_clusters/'
        with open(i) as f:                            # open the cluster.fasta to count entries by counting '>'
            n = 0
            for line in f:
                if line[0] == '>':
                    n = n + 1
        
        if n < minimum_depth:
            low_unmapped_count.append(i)             # append clusters with less than 'minimum_depth' mapped reads to placeholder list
    
    for i in low_unmapped_count:
        os.remove(i)                                 # remove the clusters with less than 'minimum_depth' mapped reads
 
print('\n')

# create the directory for among-samples clustering

in_between_samples_clustering_dir = analysis_dir + '/among_samples_clustering'         # path the directory into the analysis directory

print('Creating a working directory for the "among-samples clustering":', in_between_samples_clustering_dir, '\n')

if not os.path.exists(in_between_samples_clustering_dir):
    os.makedirs(in_between_samples_clustering_dir)                                         # create the directory if it does not exist
    
    
# list all consout

list_of_filtered_consout = []                                                                                             # empty list to append to from loop below

for i in range(len(barcodes)):                                                                                            # loop through barcodes
    filtered_consout = within_sampe_clustering_dir + '/' + barcodes[i] + '/consout_' + barcodes[i] + '_filtered.fasta'    # get path of consout_*_filtered.fasta
    list_of_filtered_consout.append(filtered_consout)                                                                     # append path of consout*filtered.fasta to list


# cat all consout_*_filtered.fasta files

cat_consout = in_between_samples_clustering_dir + '/cat_filtered_consout_wsc.fasta'

if os.path.exists(cat_consout):
    os.remove(cat_consout)

print('Concatenate', len(barcodes), 'samples to', cat_consout, '\n')
    
with open(cat_consout,'wb') as wfd:          # create and open in binary mode the cat_filtered_consout.fasta to write in
    for f in list_of_filtered_consout:       # loop through the filtered consouts
        with open(f,'rb') as fd:             # open a filtered consout
            shutil.copyfileobj(fd, wfd)      # copy content of filtered consout to the cat_filtered_consout.fasta


# create directory for among-samples cluster files            

msaout = in_between_samples_clustering_dir + '/msaout'                                # define the vsearch msaout path and name
consout = in_between_samples_clustering_dir + '/consout'                              # define the vsearch consout path and name
clusters_dir = in_between_samples_clustering_dir + '/clusters'                        # define the vsearch clusters path and directory name

if not os.path.exists(clusters_dir):                                                  # if the clusters directory does not yet exist
    os.makedirs(clusters_dir)                                                         # create it

clusters = clusters_dir + '/wsc_cluster_'                                             # define the vsearch clusters names
stdout_path = in_between_samples_clustering_dir + '/vsearch_stdout'                   # define the vsearch stdout
stderr_path = in_between_samples_clustering_dir + '/vsearch_stderr'                   # define the vsearch stderr

vsearch_input = cat_consout                                                           # define the path to the correct sample

vsearch_wsc_args = ['vsearch', 
                    '--cluster_fast', vsearch_input, 
                    '--id', in_between_samples_clustering_ct, 
                    '--msaout', msaout, 
                    '--consout', consout, 
                    '--clusters', clusters, 
                    '--threads', threads]                                             # list of the vsearch command line input formatted for the subprocess package

if os.path.exists(stdout_path):                                                       # remove vsearch stdout and stderr if it already exists from previous runs
    os.remove(stdout_path)
if os.path.exists(stderr_path):
    os.remove(stderr_path)

print('among-samples clustering: clustering filtered within-sample consout of all samples')
print('vsearch --id:', in_between_samples_clustering_ct)
print('consout:', consout, '\n')

with open(stdout_path, 'x') as stdout_file:                                                                   # create and open files to write stderr and stdout to
    with open(stderr_path, 'x') as stderr_file:
        run_vsearch_wsc = subprocess.run(vsearch_wsc_args, stdout = stdout_file, stderr = stderr_file)        # run vsearch
    
        if run_vsearch_wsc.returncode == 0:                                                                   # print subprocess.run stdout if vsearch process was successful
            print('Done. Details are written to', stderr_path, '\n')
            #print(run_vsearch_wsc, '\n')
            with open(stderr_path) as f:
                print(f.read())
        if run_vsearch_wsc.returncode == 1:                                                                   # print a warning message if vsearch process was unsuccessful
            print('Warning. Something went wrong concerning the within-sample-clustering of', barcodes[i])
            print(run_vsearch_wsc, '\n')


# count reads in every cluster = clustercount

clusters_in_dir = clusters_dir + '/*'
list_clusters = glob.glob(clusters_in_dir)                                       # list all clusters in /in_between_samples_clustering/clusters


# create clustercount.csv

clustercount_dir = in_between_samples_clustering_dir + '/clustercount.csv'       # name and path the clustercount.csv

if os.path.exists(clustercount_dir):
    os.remove(clustercount_dir)                                                  # remove clustercount.csv if it already exists
if not os.path.exists(clustercount_dir):
    with open(clustercount_dir, 'x') as f:
        f.write('cluster\tcount\n')                                              # create and open clustercount.csv if it does not yet exist and add header columns 'cluster' and 'count'

for i in list_clusters:                                                          # loop through the listed clusters
    with open(i) as f:                                                           # open the cluster
        n = 0
        n_reads = []
        for line in f:
            if line[0] == '>':
                n = n + 1                                                        # if a line starts with '>' (= fasta header), count + 1
        n_reads.append(n)                                                        # append count to n_reads variable
        
        with open(clustercount_dir, 'a') as f:                                   # open the clustercount.csv for appending n_reads of cluster i
            f.write(path_leaf(i))
            f.write('\t')
            f.write(str(n_reads[0]))
            f.write('\n')

print('Clustercount.csv created in', clustercount_dir, '\n')

df = pd.read_csv(clustercount_dir, sep='\t')
print('Total reads:\t\t\t', df['count'].sum())
print('Mean reads in cluster:\t\t', round(df['count'].mean(), 2), '+-', round(df['count'].std(), 2))
print('Singleton clusters:\t\t', df['count'].isin([1]).sum())
print('% Singleton clusters:\t\t', (round(df['count'].isin([1]).sum() / df['count'].sum(), 4) * 100))
print('\n')


# create clustercount_hist.csv

clustercount_df = pd.read_csv(clustercount_dir, sep='\t')                        # open clustercount.csv in pandas
clustercount_hist = clustercount_dir.replace('.csv', '_hist.csv')                # define clustercount_hist.csv path and name
clustercount_df.groupby('count').count().to_csv(clustercount_hist, sep="\t")     # group clustercount.csv by column 'count' and output to clustercount_hist.csv

print('Clustercount_hist.csv created in', clustercount_hist)


# filter out loci cluster with multiple occurences of a sample

df_filtered = df.loc[(df['count'] >= minimum_sample_per_locus) & (df['count'] <= len(barcodes))]        # create a dataframe filtering out below minimum sample per locus and above total
                                                                                                        # amout of sample = len(barcodes)
    
df_filtered_list = df_filtered['count'].tolist()                                                        # puts the cluster names in df_filtered to a list

# loop through list and count occurences of barcodes, if it is higher than 1, remove from list.
# the remaining clusters are completely filtered, move on to read_aligner.ipynb

print('\n')
print("Total clusters:\t\t\t", len(df.index))
print("Singleton clusters:\t\t", len(df.loc[(df['count'] < minimum_sample_per_locus)]))
print("clusters >", len(barcodes), 'samples:\t\t', len(df.loc[(df['count'] > len(barcodes))]))
print("Potentially usable clusters:\t", len(df_filtered.index))  # clusters/loci left after initial filtering

print('\n')



df_filtered_list = df_filtered['cluster'].tolist()

print('prefiltered loci:', len(df_filtered_list))


# create a csv to log occurence of samples in each loci cluster

loci_info_path = in_between_samples_clustering_dir + '/loci_info.csv'

if os.path.exists(loci_info_path):      # remove loci_info.csv if it already exists (from previous runs)
    os.remove(loci_info_path)

with open(loci_info_path, 'x') as f:    # create the loci_info.csv and prepare the header
    f.write('cluster\t')
    for i in barcodes:
        f.write(i)
        f.write('\t')
    f.write('\n')

    
# count sample occurence and log to csv    
    
for i in df_filtered_list:
    path_i = clusters_dir + '/' + i
    #print(path_i)
    with open(loci_info_path, 'a') as outfile:
        outfile.write(i)
        outfile.write('\t')
        with open(path_i) as f:
            content = f.read()
            for j in barcodes:
                count_arg = '>' + j
                count_barcodes = content.count(count_arg)
                #print(j, count_barcodes)
                outfile.write(str(count_barcodes))
                outfile.write('\t')
        outfile.write('\n')

print('\n')
        
    
# loop through list and count occurences of barcodes, if it is higher than 1, remove from list.

print('Now removing Loci with multiple occurence of a sample consout...')

to_remove_list = []

for i in df_filtered_list:
    path_i = clusters_dir + '/' + i
    #print(path_i)
    with open(path_i) as f:
        content = f.read()
        for j in barcodes:
            count_arg = '>' + j # specifies search string just in case the centroid read-id in the header randomly contains 'bcXX'
            count_barcodes = content.count(count_arg) # used 'count_arg' instead of just 'j' for the reason explained above
            #print(count_barcodes)
            if count_barcodes >= 2:
                print(i, 'is marked for removal due to multiple occurence of the same sample in loci', end='\r')
                to_remove_list.append(i)

df_filtered_list = [x for x in df_filtered_list if x not in to_remove_list]            # removes the list entries of 'to_remove_list' in 'df_filtered_list'

print('\n')
print('Loci remaining after filtering multiple sample occurences:', len(df_filtered_list))
#print(df_filtered_list)


# remove the tabstop at the end of lines

df_loci_csv = pd.read_csv(loci_info_path, delimiter='\t')
df_loci_csv.drop(df_loci_csv.columns[len(df_loci_csv.columns)-1], axis=1, inplace=True)
df_loci_csv.to_csv(loci_info_path, sep='\t', index=False)

print('\nA table documenting occurence of samples per loci is placed to', loci_info_path)


# use df_filtered_list to search the filtered clusters and move them to a new directory

clusters_filtered_dir = in_between_samples_clustering_dir + '/clusters_filtered'

if os.path.exists(clusters_filtered_dir):
    shutil.rmtree(clusters_filtered_dir)
if not os.path.exists(clusters_filtered_dir):
    os.makedirs(clusters_filtered_dir)

for i in df_filtered_list:
    src = clusters_dir + '/' + i
    dst = clusters_filtered_dir + '/' + i
    shutil.move(src, dst)
    #print(src, 'moved to', dst)
    
print('The filtered loci have been moved to', clusters_filtered_dir, '\n')


# Prepare 'infos' files, which contain sample/bc and its centroid read id, homologized by the among-samples clustering (= a locus)

infos_dir = in_between_samples_clustering_dir + '/infos'         # name and path the infos file directory

if os.path.exists(infos_dir):                                    # remove the infos file directory with content if it already exists
    shutil.rmtree(infos_dir)
if not os.path.exists(infos_dir):                                # make the infos file directory if it does not exist
    os.makedirs(infos_dir)
    
filtered_clusters_list_arg = clusters_filtered_dir + '/*'
filtered_clusters_list = glob.glob(filtered_clusters_list_arg)   # make a list with all filtered clusters

for i in range(len(filtered_clusters_list)):
    infos_dir_output = infos_dir + '/infos_' + path_leaf(filtered_clusters_list[i])
    write_infos_to_file(get_infos_from_cluster_file(filtered_clusters_list[i]), infos_dir_output)      # write infos file for each filtered cluster

info_file_list_arg = infos_dir + '/infos*'
info_file_list = glob.glob(info_file_list_arg)                   # list all infos files

list_of_dfs = [pd.read_csv(file, delimiter='\t') for file in info_file_list]                           # iterate through the infos files and list them

for dataframe, cluster_name in zip(list_of_dfs, info_file_list):                                       # add column 'cluster_name' to the info files
    dataframe['cluster_name'] = path_leaf(cluster_name)


# list all filtered nanopore reads.fastq

filtered_reads_arg = filtered_reads_dir + '*.fastq'
filtered_reads_list = glob.glob(filtered_reads_arg)


# create new directory for the filtered out fastq reads (done below), if it does not exist

filtered_fastq_dir = in_between_samples_clustering_dir  + '/filtered_fastq'      # name the directory for the fastq files with only the content of the cluster

if os.path.exists(filtered_fastq_dir):                                           # check if directory exists
    shutil.rmtree(filtered_fastq_dir)                                            # remove directory with content
if not os.path.exists(filtered_fastq_dir):                                       # if directory does not exist
    os.makedirs(filtered_fastq_dir)                                              # make the directory
    print('Filtered reads directory:', filtered_fastq_dir)


# create a file to log unmapped reads

unmapped_reads_log_dir = in_between_samples_clustering_dir + '/unmapped_reads_log.csv'

if os.path.exists(unmapped_reads_log_dir):                                                          # check if missing reads in the mapping file exists
    os.remove(unmapped_reads_log_dir)                                                               # remove it
if not os.path.exists(unmapped_reads_log_dir):                                                      # if file does not exist
    with open(unmapped_reads_log_dir, 'x') as f:                                                    # make it
        f.write('locus\tcluster\tunmapped reads\ttotal reads in cluster\tPercentage unmapped reads\n')     # and write the header
    print('Unmapped reads are logged to:', unmapped_reads_log_dir)
    

# create a dictionary of all consout consensus sequences generated after the second (among samples) clustering

consout_infile = consout
consout_dict = consout_to_dict(consout_infile)
print("Total count of second (among-samples) clustering consout consensus sequences:", len(consout_dict), "\n")
    

    

loopcounter = 0

for i in range(len(list_of_dfs)):
    loopcounter = loopcounter + 1
    print('infos_cluster number', loopcounter, 'out of', len(info_file_list), '- Currently iterating through', list_of_dfs[i]['cluster_name'][0], '...', end='\r') 
    
    # create seperate directory for the outputs of each info_file
    
    info_output_dir = filtered_fastq_dir + '/output_' + list_of_dfs[i]["cluster_name"][0]
    if os.path.exists(info_output_dir):
        shutil.rmtree(info_output_dir)
    if not os.path.exists(info_output_dir):
        os.makedirs(info_output_dir)
    
    
    # write the corresponding consout consensus sequence of the in-between-samples clustering as a fasta reference to 'info_output_dir'
    
    corresponding_consout_name = info_output_dir + '/' + list_of_dfs[i]['cluster_name'][0] + '_reference.fasta'
    fasta_reference_formatting = '>' + list_of_dfs[i]['cluster_name'][0] + '_reference\n' + consout_dict[list_of_dfs[i]['read_id'][0]]
    
    with open(corresponding_consout_name, 'x') as f:
        f.write(fasta_reference_formatting)

        
        
    # loop through all loci to identify fastq reads of each sample to map against the pseudoreference generated by within-samples clustering vsearch --consout
        
    for j in range(len(list_of_dfs[i])):                                                 # loop through lines (=samples) in the 'list_of_dfs'
        #print(j)
        for k in barcodes:                                                               # loop through the barcodes
            if list_of_dfs[i]['barcode'][j] == k:                                        # if a barcode matches the barcode of current line (=sample):
                #print(list_of_dfs[i]['barcode'][j], 'is', k)
                path_to_fastq_arg = filtered_reads_dir + k + '*.fastq'               
                path_to_fastq = glob.glob(path_to_fastq_arg)                             # choose the bcXX*.fastq Nanopore filtered read 
                if len(path_to_fastq) >= 2:
                    print('More than 1 FASTQ-file was found for sample', k, '. Please make sure to remove similar named .fastq files with the same barcode/identifier from the directory')
                    break                                                                # Stop the loop if somehow 2 files with the same barcode were chosen as 'path_to_fastq'
                if len(path_to_fastq) == 0:
                    print('No filtered reads.fastq could be identified. Make sure to place the exact barcode affilitation into the file name, e.g. "bc01.fastq"')
                    break                                                                # Stop the loop if no reads matching the barcode were identified
                if len(path_to_fastq) == 1:                                              # but continue if there is a single match for Nanopore filtered fastqs in 'path_to_fastq'
                    #print(k, ':', path_leaf(path_to_fastq[0]).replace('_q7_200to1000bp.fastq', ''))
                    wsc_clusters_arg = within_sampe_clustering_dir + '/' + list_of_dfs[i]['barcode'][j] + '/depth_filtered_clusters/*_cluster_*'
                    wsc_clusters = glob.glob(wsc_clusters_arg)                           # list all within-samples clusters of the specified sample
                    #print('Extracting', list_of_dfs[i]['barcode'][j], 'reads:')
                    
                    for cluster in wsc_clusters:                                 # iterates through clusters(.fasta) in path chosen in 'wsc_clusters'
                        with open(cluster, 'r') as f:
                            if list_of_dfs[i]['read_id'][j] in f.read():         # if the sample's centromer read_id matches the read_id in a cluster:
                                #print('Found reads in:', cluster)
                                
                                with open(cluster, 'r') as f:                    # open and work with this cluster
                                    
                                    # parse the cluster.fasta and place read ids in a list
                                    content = f.read().split('>')             # list the content by splitting the fasta by '>'
                                    while '' in content:
                                        content.remove('')                       # removes empty list contents (the first element is always empty)
                                    read_ids = []                                # create list to append read ids of current cluster to
                                    for element in content:                      # loop through each fasta entry in the cluster.fasta
                                        header = element[0:element.find('\n')]   # takes the header = Nanopore read id 
                                        read_ids.append(header)                  # append the header/read id to the read_ids list
                                    #print('Read ids in cluster:', len(read_ids))

                                    
                                    # parse the filtered reads.fastq for the current sample
            
                                    fastq_reads = []                             # an empty list to append fastq reads to (see below)

                                    with open(path_to_fastq[0]) as fq_f:         # open the filtered reads.fastq of the current sample
                                        for entry in take_four(fq_f):
                                            fastq_reads.append(''.join(entry))   # append each fastq entry as a single element to list 'fastq_reads'
                                    
                                    named_filtered_fastq = info_output_dir + '/' + list_of_dfs[i]['barcode'][j]
                                    
                                    with open(named_filtered_fastq, 'x') as f:
                                        for fastq_entry in fastq_reads:
                                            for read_id in read_ids:
                                                if read_id in fastq_entry:
                                                    f.write(fastq_entry)
                                        
                                        #print('Done writing the fastq of', cluster)
                                    
                                    
                                    # mapping reads extracted just now to pseudoreference (within-samples clustering consout)
                                    
                                    corresponding_consout_name   # reference
                                    named_filtered_fastq         # input fastq reads
                                    
                                    stderr_path = named_filtered_fastq + '_stderr'     # stderr output from mapping = minimap2 printout
                                    stdout_path = named_filtered_fastq + '.sam'        # mapping output .sam file
                                    
                                    if os.path.exists(stderr_path):    # remove mapping outputs from previous runs
                                        os.remove(stderr_path)
                                    if os.path.exists(stdout_path):
                                        os.remove(stdout_path)
                                    
                                    with open(stderr_path, 'x') as stderr_outfile:
                                        with open(stdout_path, 'x') as stdout_outfile:
                                            mapping = subprocess.run(['minimap2', 
                                                                      '-ax', 
                                                                      'map-ont', 
                                                                      '-t', threads, 
                                                                      corresponding_consout_name,  # reference path
                                                                      named_filtered_fastq],       # input fastq path
                                                                      stderr = stderr_outfile,     # minimap2 printout
                                                                      stdout = stdout_outfile      # minimap2 .sam output
                                                                    )
                                    
                                    if mapping.returncode == 1:
                                        print('Error. Subprocess for minimap2 mapping returncode = 1. Please check your output.')
                                        break
                                    #if mapping.returncode == 0:
                                    #    print('Successfully mapped', path_leaf(named_filtered_fastq), 'to', path_leaf(corresponding_consout_name))
                                        
                                    
                                    # counting unmapped reads per sample by comparing n-2 lines of the .sam file ('stdout_path') to total reads ('read_ids')
                                    # but first, remove unmapped reads and supplemental alignments (duplicated reads) from the minimap2 .sam output
                                    
                                    filtered_sam = stdout_path.replace('.sam', '_filtered.sam')
                                    filter_sam_process = subprocess.run(['samtools', 'view', 
                                                                                     '-hF', '2052', 
                                                                                     '-o', filtered_sam, 
                                                                                     stdout_path])
                                    
                                    if filter_sam_process.returncode == 1:
                                        print('Error. Subprocess for samtools view to remove supplemental alignments and unmapped reads: Returncode = 1')
                                        break
                                    #if remove_unmapped_reads_from_sam.returncode == 0:    # if samtools view worked, replace original .sam with new .sam
                                    #    os.remove(stdout_path)
                                    #    os.rename(no_unmapped_reads_sam, stdout_path)
                                    
                                    with open(filtered_sam) as sam_f:
                                        num_lines = sum(1 for line in sam_f)
                                        mapped_reads = num_lines - 2
                                        
                                    with open(unmapped_reads_log_dir, 'a') as unmapped_reads_f:
                                        unmapped_reads_f.write(list_of_dfs[i]['cluster_name'][0])
                                        unmapped_reads_f.write('\t')
                                        unmapped_reads_f.write(path_leaf(cluster))
                                        unmapped_reads_f.write('\t')
                                        unmapped_reads_f.write(str(len(read_ids) - mapped_reads))
                                        unmapped_reads_f.write('\t')
                                        unmapped_reads_f.write(str(len(read_ids)))
                                        unmapped_reads_f.write('\t')
                                        unmapped_reads_f.write(str((len(read_ids) - mapped_reads) / len(read_ids)))
                                        unmapped_reads_f.write('\n')        
                    
    #print('\nDone with', list_of_dfs[i]['cluster_name'][0], '\n-------------------------------------------------------------------------------------------------------------------------\n')

print('\n')

# unmapped reads stats

df = pd.read_csv(unmapped_reads_log_dir, delimiter='\t')
print('\nUnmapped reads:\t\t\t\t\t', df['unmapped reads'].sum(), 'out of', df['total reads in cluster'].sum(), 
      '(', round((df['unmapped reads'].sum() / df['total reads in cluster'].sum()) * 100, ndigits=2), '% )')
print('Clusters with unmapped reads:\t\t\t', len(df[(df['unmapped reads'] >= 1)]), 'out of', len(df), 
      '(', round((len(df[(df['unmapped reads'] >= 1)]) / len(df)) * 100, ndigits=2), '% )')
print('Sum of clusters with', minimum_depth, 'or more unmapped reads:\t', len(df[(df['unmapped reads'] >= minimum_depth)]), 
      '(', round((len(df[(df['unmapped reads'] >= minimum_depth)]) / len(df)) * 100, ndigits=2))
print('% clusters with', minimum_depth, 'or more unmapped reads:\t', round(len(df[(df['unmapped reads'] >= minimum_depth)]) / len(df), ndigits=2))
print('A detailed overview can be found in', unmapped_reads_log_dir)
    
print('\n')


# sam to bam conversion, bam sorting and bam indexing

sam_files_list_arg = filtered_fastq_dir + '/*/*_filtered.sam'
sam_files_list = glob.glob(sam_files_list_arg)

print('.sam to .bam conversion, bam sorting and indexing of')

loopcounter = 0

for i in sam_files_list:
    
    loopcounter = loopcounter + 1
    print(loopcounter, 'out of', len(sam_files_list), ':', i, end='\r')
    
    bam_output = i.replace('_filtered.sam', '.bam')
    sam_to_bam = subprocess.run(['samtools', 'view', '-b', '-F', '12', '-o', bam_output, i])
    if sam_to_bam.returncode == 1:
        print('something went wrong: sam to bam')
        break
    
    sorted_bam_output = bam_output.replace('.bam', '_sorted.bam')
    sort_bam = subprocess.run(['samtools', 'sort', '-o', sorted_bam_output, bam_output])                      # run samtools sort -o [out] [in]
    if sort_bam.returncode == 1:
        print('something went wrong: bam sorting')
        break
    
    bai_output = sorted_bam_output.replace('.bam', '.bai')                                                    # the .bai output should habe the same name and directory as the .bam
    index_bam = subprocess.run(['samtools', 'index', '-b', sorted_bam_output, bai_output])                    # run samtools index -b [in] [out]
    if index_bam.returncode == 1:
        print('something went wrong: bam indexing')
        break    
    
    os.remove(bam_output)                                                                                     # removes the unsorted bam file

print('\n')


# collecting SNPs loci per loci with bcftools

# create a log file counting the amount of SNPs called for each locus
snp_log_file_dir = in_between_samples_clustering_dir +  '/snps_log.csv'  # name and path for the log file
if os.path.exists(snp_log_file_dir):                                                                 # check if log file exists
    os.remove(snp_log_file_dir)                                                                      # yes: remove it
if not os.path.exists(snp_log_file_dir):                                                             # if file does not exist:
    with open(snp_log_file_dir, 'x') as snp_log_file:                                                # creates the log file
        snp_log_file.write('locus')                                                                  # write a header for tab-delimited text file
        snp_log_file.write('\t')
        snp_log_file.write('SNPs')
        snp_log_file.write('\n')

loci_list_arg = filtered_fastq_dir + '/*'
list_of_all_loci = glob.glob(loci_list_arg)                      # collect all loci directories in order to loop through them

loopcounter = 0

print('Started calling SNPs for each loci:')

for i in range(len(list_of_all_loci)):                                                          # loop through length of loci directories (= list_of_all_loci)

    loopcounter = loopcounter + 1                                                               # counts + 1 every time the loop is started again
    print(loopcounter, "out of", len(list_of_all_loci), ':', list_of_all_loci[i], end='\r')     # shows progress of this loop
    
    
    '''bcftools mpileup'''
    '''Collects all possible variants found in the mapping against the (pseudo)reference = consout from vsearch 2nd (among-samples) clustering'''
    
    bam_in_dir = glob.glob(list_of_all_loci[i] + '/*_sorted.bam')                               # list filtered and sorted .bam files in current loci directory (= list_of_all_loci[i])
    vcf_calling_ref = glob.glob(list_of_all_loci[i] + '/*_reference.fasta')                     # list the (pseudo)reference (vsearch consout of 2nd among-samples clustering) 
                                                                                                # in current loci directory (= list_of_all_loci[i])
    
    if not len(vcf_calling_ref) == 1:
        print('More than one possible reference.fasta was found in', list_of_all_loci[i], 'make sure there is only one file called "_reference.fasta".')
        break
    
    mpileup_out = vcf_calling_ref[0].replace('_reference.fasta', '_mpileup.vcf')                # correctly names the bcftools mpileup file + path
    
    # run bcftools mpileup
    bcftools_mpileup_args = ['bcftools', 'mpileup', 
                             '--skip-indels', 
                             '--count-orphans', 
                             '--no-BAQ', 
                             '--min-BQ', '0', 
                             '-Ov', 
                             '--annotate', 'DP,AD,DP4', 
                             '-f', vcf_calling_ref[0], 
                             '-o', mpileup_out] + bam_in_dir
    run_bcftools_mpileup = subprocess.run(bcftools_mpileup_args)
    
    if run_bcftools_mpileup.returncode == 1:                                                    # print a warning message in case there is an error creating the mpileup vcf file
        print('Error. bcftools mpileup: returncode =/= 0!\t')
        break

        
    '''Change sample name'''
    '''Samples in the mpileuxp.vcf are now named after the path of the input _filtered_sorted.bam file. This will cause problems in subsequent sample merging. Names will therefore
    be changed to simply the barcode number.'''
    
    mpileup_edit_out = mpileup_out.replace('_mpileup.vcf', '_mpileup_edit.vcf')                    # defines the output path/name of the edited mpileup file with correct sample names
       
    with open(mpileup_out) as f:                                                                   # opens the *_mpileup.vcf file
        correct_sample_names = f.read()                                                            # loads the text in *_mpileup.vcf
        for j, k in enumerate(bam_in_dir):                                                         # loops through .bam files in 'current_locus' list
            if bam_in_dir[j] not in correct_sample_names:                                          # print a warning message that sample in current_locus[j] was not found in mpileup.vcf
                print('COULD NOT FIND', bam_in_dir[j])
            else:
                correct_sample_names = correct_sample_names.replace(bam_in_dir[j], path_leaf(k))   # if sample name was found in *_mpileup.vcf, replace path name for just the file name
        correct_sample_names = correct_sample_names.replace('_sorted.bam', '')            # previous task leaves *_filtered_sorted.bam in sample name, remove this so only
                                                                                                   # the barcode remains as the sample name (e.g. bc01 instead of bc01_filtered_sorted.bam)
        if os.path.exists(mpileup_edit_out):
            os.remove(mpileup_edit_out)                                                            # remove the *_mpileup_edit.vcf file if it already exists (from a previous run)

    with open(mpileup_edit_out, 'x') as samples_named_vcf:                                         # create and open the *mpileup_edit.vcf file
        samples_named_vcf.write(correct_sample_names)                                              # write the corrected sample names in *_mpileup_edit.vcf
    
    
    '''bcftools norm -m -'''
    '''Splits multiallelic sites so every variant gets its own vcf entry. This makes it easier to filter for true variants in the next steps'''
    
    norm_split_out = mpileup_out.replace('_mpileup.vcf', '_norm_split.vcf')                     # correctly names the bcftools norm -m - output file and path
    
    # run bcftools norm -m -
    bcftools_norm_split_args = ['bcftools', 'norm', '-m', '-', '-Ov', '-o', norm_split_out, mpileup_edit_out]
    run_bcftools_norm_split = subprocess.run(bcftools_norm_split_args)
    
    if run_bcftools_norm_split.returncode == 1:                                                 # print a warning message in case there is an error creating the norm -m - vcf file
        print('Error. bcftools norm -m -: returncode =/= 0!\t')
        break
    
    
    '''bcftools view -i (AD[*:1]/FORMAT/DP) > 0.25 & MIN(FORMAT/DP) >= 5'''
    '''Filters for true variants. Variants rarer than 25% and at positions with less than 5 coverage will be filtered out. Therefore SNPs must appear at least twice at minimum coverage. 
    Only applicable for diploid samples'''
    
    view_out = norm_split_out.replace('_norm_split.vcf', '_view.vcf')                           # correctly names the bcftools view output file and path
    
    # run bcftools view -i "(AD[*:1]/FORMAT/DP) > 0.25 & MIN(FORMAT/DP) >= 5"
    bcftools_view_args = ['bcftools', 'view', '-i', '(AD[*:1]/FORMAT/DP) > 0.25 & MIN(FORMAT/DP) >= 5', '-o', view_out, norm_split_out]
    run_bcftools_view = subprocess.run(bcftools_view_args)

    if run_bcftools_view.returncode == 1:                                                       # print a warning message in case there is an error creating the view vcf file
        print('Error. bcftools view: returncode =/= 0!\t')
        break

        
    '''bcftools norm -m +'''
    '''Combines multiallelic sites, so every site/position will be a single entry in the finalized vcf for the locus'''
    
    locus_vcf_out = view_out.replace('_view.vcf', '_locus.vcf')                                 # correctly names the bcftools norm -m + output file and path
    
    # run bcftools norm -m +
    bcftools_norm_combine_args = ['bcftools', 'norm', '-m', '+', '-Ov', '-o', locus_vcf_out, view_out]
    run_bcftools_norm_combine = subprocess.run(bcftools_norm_combine_args)
    
    if run_bcftools_norm_combine.returncode == 1:                                               # print a warning message in case there is an error creating the view vcf file
        print('Error. bcftools norm -m +: returncode =/= 0!\t')
        break

        
    # count amount of SNPs of current locus and append to the log file
    linecount = 0
    with open(locus_vcf_out) as f:
        for line in f:
            if not line.startswith('#'):
                linecount += 1
                
        with open(snp_log_file_dir, 'a') as snp_log_file:
            snp_log_file.write(path_leaf(list_of_all_loci[i]))
            snp_log_file.write('\t')
            snp_log_file.write(str(linecount))
            snp_log_file.write('\n')           

            
# printed outputs

print('\r\nDone with SNP calling.\n')
snps_df = pd.read_csv(snp_log_file_dir, sep='\t')                                      # opens the snp_log file for statistics
print('Total amount of SNPs:', snps_df['SNPs'].sum())                                  # prints total SNPs
print('Mean SNPs per locus:', snps_df['SNPs'].mean(), '+/-', snps_df['SNPs'].std())    # prints mean SNPs per locus + standard deviation
print('Median SNPs per locus:', snps_df['SNPs'].median())                              # prints median SNPs per locus
print('Loci with 0 SNPs:', snps_df['SNPs'].isin([0]).sum(), '\n')                      # prints amount of loci with 0 SNPs


list_of_all_vcf_arg = filtered_fastq_dir + '/*/*_locus.vcf'
list_of_all_vcf = glob.glob(list_of_all_vcf_arg)                                    # lists all locus VCF files
vcfcombine_out = in_between_samples_clustering_dir + '/Analysis_SNPS_dotted.vcf'    # output path


# building the subprocess.run *popenargs list with the command line for vcflib's vcfcombine 

vcfcombine_args = ['vcfcombine']            # start with the vcfcombine module
vcfcombine_args.extend(list_of_all_vcf)     # add the list with the paths for all locus VCF files

if os.path.exists(vcfcombine_out):          # delete the merged VCF file if it already exists from a previous run
    os.remove(vcfcombine_out)

with open(vcfcombine_out, 'x') as vcfcombine_stdout:                                  # create and open the file for the merged VCF output
    run_vcfcombine = subprocess.run(vcfcombine_args, stdout = vcfcombine_stdout)      # runs vcfcombine and outputs in stdout = file made above with path from vcfcombine_out

if run_vcfcombine.returncode == 1:                                                                                                # if the subprocess failed:
    print('Error. vcfcombine returncode =/= 0!\t')                                                       # print a warning message
else:                                                                                                                             # if everything worked:
    print('All', len(list_of_all_vcf), 'loci VCF files have been merged to a single one:', vcfcombine_out)    # print some infos about the run

      
# replace '.' in VCF for 0,0,0:0:0,0,0,0:0,0

no_dot_vcf_outfile = vcfcombine_out.replace('Analysis_SNPS_dotted.vcf', 'analysis_SNPs.vcf')

c = 0

with open(no_dot_vcf_outfile, 'w') as of:
    with open(vcfcombine_out) as f:
        for line in f:
            if line[0] == '#':
                of.write(line)
                #print(line)
            else:
                cols = line.split('\t')
                cols[-1] = cols[-1].rstrip()
                for i in range(9, len(cols)):
                    if cols[i][0] == '.':
                        cols[i] = '0,0,0:0:0,0,0,0:0,0'
                of.write('\t'.join(cols) + '\n')
                
print('\nFinished at:', datetime.now().strftime('%d/%m/%Y %H:%M:%s'))
print('Duration: %s minutes' % round(((time.time() - start_time) / 60), ndigits=2))
