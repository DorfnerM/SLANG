"""
SLANG VCF to SNP matrix.py
This tool parses VCF files outputted by SLANG and computes a SNP matrix.
"""

__version__ = '0.1'
__author__ = 'Marco Dorfner'
__email__ = 'marco.dorfner@ur.de'
__date__ = '2022-10-18'


import argparse
import pandas as pd
import numpy as np
import os
import time
from pathlib import Path


class SLANG_vcf_to_SNP_matrix:
    def __init__(self):
        self.parser = argparse.ArgumentParser(
            prog = 'SLANG VCF to SNP matrix',
            description = 'This tool parses VCF files outputted by SLANG and computes a SNP matrix.'
        )
        self.CLI()
        self.args = self.parser.parse_args()

        self.run()

    def CLI(self): 
        self.parser.add_argument(
            'vcf',
            type = str,
            help = 'path to input SLANG vcf file'
        )

        self.parser.add_argument(
            '-o', '--out',
            metavar = '',
            type = str,
            help = 'output path and name. Default in current directory.'
        )

        self.parser.add_argument(
            '-p', '--ploidy',
            metavar = '',
            type = int,
            default = 2,
            help = 'ploidy level of samples. Default = 2 (diploid)'
        )

        self.parser.add_argument(
            '-v', '--version',
            action = 'version',
            version = __version__
        )

    def vcf_to_pandas_df(self):
        '''
        Parses the vcf from SLANG and returns a pandas dataframe 
        storing relevant data for the SNP matrix generation.
        '''

        # extract relevant information from vcf
        chrom = list()
        pos = list()
        ref = list()
        alt = list()
        samples = list()
        
        with open(self.args.vcf) as infile:
            for line in infile:
                if line.startswith('##'):
                    continue # skip header lines
                else:       
                    columns = line.split('\t')

                    chrom_val = columns[0].strip()
                    pos_val = columns[1].strip()
                    ref_val = columns[3].strip()
                    alt_val = columns[4].strip()
                    sample_val = [__.strip() for __ in columns[9:]]

                    chrom.append(chrom_val)
                    pos.append(pos_val)
                    ref.append(ref_val)
                    alt.append(alt_val)
                    samples.append(sample_val)

        data = dict()
        data['chrom'] = chrom[1:]
        data['pos'] = pos[1:]
        data['ref'] = ref[1:]
        data['alt'] = alt[1:]

        # add 'samples' to 'data' with AD base frequencies for each locus
        for sample in samples[0]:
            sample_format = list()
            for format_fields in samples[1:]:
                idx = samples[0].index(sample)
                format = format_fields[idx].strip()
                if format != '.':
                    format = format.split(':')[1] + ',' + format.split(':')[-1]
                sample_format.append(format)
            data[sample] = sample_format

        # make pandas df from 'data'
        column_names = ['chrom', 'pos', 'ref', 'alt']
        for sample in samples[0]:
            column_names.append(sample)

        df = pd.DataFrame(
            data,
            columns = column_names
        )

        return df
    
    def IUPAC(self, nucleotides):
        '''Finds the correct IUPAC nucleotide code from input nucleotides.'''

        input = np.sort(nucleotides)

        iupac_code_2bp = {
            'R' : np.array(['A', 'G']),
            'Y' : np.array(['C', 'T']),
            'S' : np.array(['C', 'G']),
            'W' : np.array(['A', 'T']),
            'K' : np.array(['G', 'T']),
            'M' : np.array(['A', 'C'])
        }

        iupac_code_3bp = {
            'B' : np.array(['C', 'G', 'T']),
            'D' : np.array(['A', 'G', 'T']),
            'H' : np.array(['A', 'C', 'T']),
            'V' : np.array(['A', 'C', 'G'])
        }

        if len(input) == 2:
            for wobble, code in zip(iupac_code_2bp.keys(), iupac_code_2bp.values()):
                if np.array_equal(input, code) is True:
                    return wobble
                elif np.array_equal(input, code) is False:
                    continue
                else:
                    print('No IUPAC nucleotide could be chosen from {} [2bp]'.format(nucleotides))
        elif len(input) == 3:
            for wobble, code in zip(iupac_code_3bp.keys(), iupac_code_3bp.values()):
                if np.array_equal(input, code) is True:
                    return wobble
                elif np.array_equal(input, code) is False:
                    continue
                else:
                    print('No IUPAC nucleotide could be chosen from {} [3bp]'.format(nucleotides))
        elif len(input) == 4:
            return 'N'
        else:
            print('No IUPAC nucleotide could be chosen from {}'.format(nucleotides))
    
    def get_allele_threshold(self, ploidy):
        '''Calculate the allele acceptance threshold for the given ploidy.
        If the base frequency of a variant > allele threshold, the variant
        is accepted as a SNP.'''
        allele_threshold =  (1 / ploidy) - ((1 / ploidy) * 0.2)
        return allele_threshold
    
    def SNP_chooser(self, allele_threshold, depth, base_frequencies, variants):
        '''Chooses the correct SNP from base frequencies.'''

        base_frequencies = np.asarray(base_frequencies)
        chosen_SNPs = np.empty(0)

        for base_frequency, variant in zip(base_frequencies, variants):
            div = base_frequency / depth
            if div >= allele_threshold:
                chosen_SNPs = np.append(chosen_SNPs, variant)

        if chosen_SNPs.size == 0:
            chosen_SNPs = np.append(chosen_SNPs, 'N')

        return chosen_SNPs
    
    def run(self):
        '''main function'''

        start = time.time()
        df = self.vcf_to_pandas_df()
        columns = df.columns.to_list()

        if os.path.exists(Path(self.args.out)):
            os.remove(Path(self.args.out))
        with open(Path(self.args.out), 'x') as outfile:
            outfile.write('{}\n'.format(len(columns[4:])))

        allele_threshold = self.get_allele_threshold(ploidy = self.args.ploidy)

        for sample in columns[4:]:
            snps = np.array([sample])
            idx = columns.index(sample)

            for row in df.iterrows():
                format = np.asarray(row[1][idx].split(','))
                if str(format[0]) == '.':
                    snps = np.append(snps, 'N')
                else:
                    depth = int(format[0])
                    base_frequencies = np.array([int(i) for i in format[1:]])

                    ref = row[1][2]
                    alt = np.asarray(row[1][3].split(','))
                    variants = np.insert(alt, 0, ref)
                    filtered_variants = self.SNP_chooser(allele_threshold, depth, base_frequencies, variants)

                    if len(filtered_variants) == 1:
                        snps = np.append(snps, filtered_variants)
                    else:
                        snps = np.append(snps, self.IUPAC(filtered_variants))

            __ = snps.tolist()
            snps = ''.join(__[1:])

            with open(Path(self.args.out), 'a') as outfile:
                outfile.write('{}\t{}\n'.format(sample, snps))

        end = time.time()
        print('Finished. Time: {}'.format(end - start))
            

def main():
    SLANG_vcf_to_SNP_matrix()

if __name__ == "__main__":  
    main()