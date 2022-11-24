"""
remove missing data from (SLANG) vcf
A CLI tool to parsimonially remove samples with missing data 
from a (SLANG) vcf file.

Planned for further releases:
- vcf filtering without the need of a distance matrix file (.dst/.dist)
"""

__version__ = '0.1'
__author__ = 'Marco Dorfner'
__email__ = 'marco.dorfner@ur.de'
__date__ = '2022-11-24'


import argparse
import os
import pandas as pd


class rm_missing_data_from_vcf:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            prog = 'rm_missing_data_from_vcf.py',
            description = 'CLI tool to parsimonially remove samples with missing data from a (SLANG) vcf file.'
        )
        self.CLI()
        self.args = self.parser.parse_args()

        self.run()

    def CLI(self):
        self.parser.add_argument(
            '-vcf',
            type = str,
            required = True,
            metavar = '',
            help = 'path to input (SLANG) vcf file.'
        )

        self.parser.add_argument(
            '-dst',
            type = str,
            required = True,
            metavar = '',
            help = 'path to input tab-delimited distance matrix (usually .dst) file without missing data'
        )

        self.parser.add_argument(
            '-o', '--out',
            metavar = '',
            type = str,
            default = os.path.join(os.path.abspath(os.getcwd()), 'output.vcf'),
            help = 'output path and name. Default in current directory.'
        )

        self.parser.add_argument(
            '-msv',
            metavar = '',
            type = int,
            default = 2,
            help = 'Minimum samples per variant needed to be accepted. Default: 2'
        )

        self.parser.add_argument(
            '-v', '--version',
            action = 'version',
            version = __version__
        )

    def get_remaining_samples(self, infile_dst):
        '''Lists samples from "filtered_dst" without missing data.
        Returns the sample list.'''

        samples = []
        with open(infile_dst) as f:
            for line in f.readlines()[1:]:
                samples.append(line.split('\t')[0])
        
        return samples

    def vcf_to_df(self, infile_vcf):
        '''Read the input vcf and transform it into a pandas df.
        Returns the dataframe.'''

        # get rows in vcf file without header
        rows = []
        with open(infile_vcf) as infile:
            for line in infile:
                if line.startswith('##'):
                    continue # skip header lines
                else:       
                    rows.append(line.split('\t'))

        # format columns as dict and make pandas df
        data = {}
        for column_idx, header in enumerate(rows[0]):
            header = header.strip()
            data[header] = [row[column_idx].strip() for row in rows[1:]]
        df = pd.DataFrame(data)

        return df

    def filter_missing_data_samples_from_vcf(self, vcf_df, infile_dst):
        '''Removes the samples with missing data from the vcf.
        Returns filtered vcf as a pandas df.'''
        
        all_samples = vcf_df.columns.tolist()[9:]
        samples_to_remove = self.get_remaining_samples(infile_dst)
        samples_to_keep = [i for i in all_samples if i not in samples_to_remove]
        filtered_vcf_df = vcf_df.drop(samples_to_keep, axis=1)
        
        return filtered_vcf_df

    def filter_redundant_variants_from_vcf(self, filtered_vcf_df, msv):
        '''Removes redundant variant lines in vcf after filtering with
        "filter_missing_data_samples_from_vcf()". A variant line is
        redundant, when < msv.
        Returns filtered vcf as a pandas df.'''
        
        for row in filtered_vcf_df.iterrows():
            variants = row[1].tolist()[9:]
            row_index = row[0]
            samples_with_variants = len([i for i in variants if i != '.'])
            if samples_with_variants < msv:
                filtered_vcf_df.drop(row_index, axis='index', inplace=True)
        
        return filtered_vcf_df

    def run(self):
        infile_vcf = self.args.vcf
        infile_dst = self.args.dst
        msv = self.args.msv

        vcf_df = self.vcf_to_df(infile_vcf)
        filtered_vcf_df = self.filter_missing_data_samples_from_vcf(vcf_df, infile_dst)
        filtered_vcf_df = self.filter_redundant_variants_from_vcf(filtered_vcf_df, msv)

        filtered_vcf_df.to_csv(self.args.out, sep='\t', index=False)

        # print some data after finishing
        samples_before = vcf_df.shape[1]-9
        samples_after = filtered_vcf_df.shape[1]-9
        variants_before = vcf_df.shape[0]-1
        variants_after = filtered_vcf_df.shape[0]-1
        variants_remaining_pct = round(((filtered_vcf_df.shape[0]-1) / (vcf_df.shape[0]-1)) * 100, 2)
        print(f'{samples_after} samples of {samples_before} remaining.')
        print(f'{variants_after} variants of {variants_before} variants ({variants_remaining_pct} %) remaining.')

def main():
    rm_missing_data_from_vcf()

if __name__ == "__main__":  
    main()
