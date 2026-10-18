import os
import pandas as pd
import sys


read_file = sys.argv[1]
output_folder = sys.argv[2]
HTF_gene = sys.argv[3]
chrom_file = sys.argv[4]


### This script will make wig wiles from pyReadCounters GTF files
### showing pileup of deletions


def process_gene_data(read_file):
    df = pd.read_csv(read_file, sep='\t', header=None, skiprows=11, usecols = [6, 8])
    location_list_fwd = []
    location_list_rev = []

    for strand, attribute in zip(df[6].tolist(), df[8].tolist()):
        attribute_list = attribute.split('; ')
        if attribute_list[-1].find('#') != -1:
            deletion_list = attribute_list[-1].split('# ')[-1].split(',')
            for deletion in deletion_list:
                location = int(deletion.replace('D', '').replace(';', ''))
                if strand == '+':
                    location_list_fwd.append(location)
                else:
                    location_list_rev.append(location)

    wig_df_fwd = pd.DataFrame({'location': location_list_fwd})
    wig_df_rev = pd.DataFrame({'location': location_list_rev})

    wig_df_fwd['value'] = wig_df_fwd.groupby('location')['location'].transform('size')
    wig_df_rev['value'] = wig_df_rev.groupby('location')['location'].transform('size')

    wig_df_fwd = wig_df_fwd.drop_duplicates()
    wig_df_fwd = wig_df_fwd.sort_values(['location'])
    new_row = pd.DataFrame({'location': ['variableStep'], 'value': ['chrom=NC_000964.3']})
    wig_df_fwd = pd.concat([new_row, wig_df_fwd])
    wig_df_rev = wig_df_rev.drop_duplicates()
    wig_df_rev = wig_df_rev.sort_values(['location'])
    wig_df_rev = pd.concat([new_row, wig_df_rev])

    return wig_df_fwd, wig_df_rev


wig_df_fwd, wig_df_rev = process_gene_data(read_file)

wig_df_fwd.to_csv(os.path.join(output_folder, f'{HTF_gene}_deletions_fwd.wig'), sep='\t', header=None, index=False)
fwd_file = output_folder + f'/{HTF_gene}_deletions_fwd.wig'
wig_df_rev.to_csv(os.path.join(output_folder, f'{HTF_gene}_deletions_rev.wig'), sep='\t', header=None, index=False)
rev_file = output_folder + f'/{HTF_gene}_deletions_rev.wig'

os.system('wigToBigWig {} {} {}/{}.bw'.format(fwd_file, chrom_file, output_folder, os.path.basename(fwd_file)[:-4]))
os.system('wigToBigWig {} {} {}/{}.bw'.format(rev_file, chrom_file, output_folder, os.path.basename(rev_file)[:-4]))


