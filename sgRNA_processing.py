#! /usr/bin/python3

# System imports
import os
import sys
import json
# Third-party imports
import pandas as pd
# Project imports

__author__ = 'Peeyush Sahu'

"""
This script will perform the sgRNA characterization after alignment.
"""
sample_name = sys.argv[1]
gff_path = sys.argv[2]
sam_path = sys.argv[3]
outpath = sys.argv[4]


def gff_parser(annot):
    all_features = {}
    for i,r in annot.iterrows():
        all_features[i] = {ent.strip().split(' ')[0]: ent.strip().split(' ')[1].strip('"') for ent in r['features'].rstrip(';').split(';')}
    return all_features


def process_gff(path):
    """
    Read and prepare file for analysis
    """
    print("Step 1 of 3: Processing gff file {}".format(path))
    annot = pd.read_csv(path, skiprows=5, header=None, index_col=None, sep='\t')
    #print(annot.shape)
    annot = annot[annot[2].isin(['gene'])]
    annot.reset_index(inplace=True, drop=True)
    annot.columns = ['chrom', 'agen', 'element', 'start', 'stop', 'score', 'strand', 'x', 'features']
    #print(annot.shape)
    #annot.head()

    annot_copy = annot.copy(deep=True)
    mod_annot = gff_parser(annot)
    extracted_features = pd.DataFrame(mod_annot).T
    extracted_features = extracted_features.fillna('-')
    #extracted_features.head()

    annot_copy = pd.concat([annot, extracted_features], axis=1)
    annot_copy = annot_copy.loc[:,['chrom', 'start', 'stop', 'gene_name']]
    annot_copy['start'] = annot_copy['start'].astype('int')

    #annot_copy.to_csv("./filtered_gff.tsv", header=None, index=None, sep='\t')
    return annot_copy


def get_flags(flag):
    """
    Decode sam alignment flags.
    """
    flag_dict = {
        0: 'paired',
        1: 'mapped_in_proper_pair',
        2: 'unmapped',
        3: 'mate_unmapped',
        4: 'reverse_strand',
        5: 'mate_reverse_strand',
        6: 'first_in_pair',
        7: 'second_in_pair',
        8: 'not_primary_alignment',
        9: 'fails_platform/vendor_quality_checks',
        10: 'PCR_or_optical_duplicate',
        11: 'supplementary_alignment',
    }
    flag = int(flag)
    if flag > 0:
        get_inx = [i for i, v in enumerate(list(bin(flag)[::-1])) if v == '1']
        decoded_flag = [flag_dict[flg] for flg in get_inx]
        return ';'.join(decoded_flag)
    else:
        return 'forward_strand'


def check_seq_annotations(sam, annot_grp):
    """
    Find fasta sequences positon on genome and identify the overlapping gene. Next check if the
    given gene name in fasta file matches with found name. Returns a dataframe with sam and annotation data combined.
    """
    print("Step 2 of 3: Assigning annotation to aligned sgRNAs.")
    sam_ext = sam.copy(deep=True)
    sam_ext = sam_ext.reindex(columns=range(0,11), fill_value=None)
    sam_grp = sam.groupby(2)
    sam_col = ['seq_id_sam', 'flag_sam', 'chrom_sam', 'start_sam', 'y', 'cigar_sam', 'gene_name_sam', 'chrom_gff', 'start_gff', 'stop_gff', 'gene_name']
    for sam_ch, sam_df in sam_grp:
        if sam_ch in annot_grp.groups.keys():
            for ind, row in sam_df.iterrows():
                annt_temp = annot_grp.get_group(row[2])
                # Running binary search to find the match between seq aligned location and gene from gff
                df_ind_l = annt_temp['start'].searchsorted(row[3], side='left')
                #df_ind_r = annt_temp['start'].searchsorted(row[3], side='right')
                #if (df_ind_r - df_ind_l) > 1:
                #    print(annt_temp.iloc[df_ind_l:df_ind_r,:])
                #    raise Warning('sgRNA {} is present in more than on gene region'.format(row[0]))
                annot_row = annt_temp.iloc[df_ind_l-1,:]
                sam_ext.iloc[ind, 7] = annot_row['chrom']
                sam_ext.iloc[ind, 8] = int(annot_row['start'])
                sam_ext.iloc[ind, 9] = int(annot_row['stop'])
                sam_ext.iloc[ind, 10] = annot_row['gene_name']
        else:
            print("No annotation entry for chromosome {} found in gff, this chromosome will be ignored.".format(sam_ch))
    sam_ext.iloc[:,1] = sam_ext.iloc[:,1].apply(lambda x: get_flags(x))
    sam_ext.columns = sam_col
    return sam_ext


def add_tcga_data(matched_sam, tcga_path):
    """
    This will add all the gene expression (unstranded fpkm) from TCGA files.
    An api call can be added to download TCGA files but out of scope for this task.
    """
    print("Step 3 of 3: Extracting norm expression from TCGA for matched genes marked my sgRNA.")
    matched_sam_temp = matched_sam.copy(deep=True)
    selected_matched_genes = list(set(matched_sam['gene_name']))
    for path, subdirs, files in os.walk(tcga_path):
        #print(path, subdirs, files)
        for name in files:
            rna_file = os.path.join(path, name)
            rna_data = pd.read_csv(rna_file, header=None, index_col=None, sep='\t', skiprows=6)
            print(rna_data.head())
            tnx_df = rna_data[rna_data[1].isin(selected_matched_genes)].iloc[:,[1,6]]
            tnx_df.columns = ['gene_name', name]
            matched_sam_temp = matched_sam_temp.merge(tnx_df, how='left', on='gene_name')
    print(matched_sam_temp.head())
    return matched_sam_temp


# Start the execution
outinfo = {}

# Read annotation file
print("Starting sgRNA processing for sam file: {}".format(sam_path))
annot_copy = process_gff(path=gff_path)
annot_copy = annot_copy.sort_values(['chrom', 'start'], ascending=[False, True])
annot_copy_grp = annot_copy.groupby('chrom')


# Read SAM file from bowtie2
# A function can be written to get the number of gff header length for skiprows
sam_file = pd.read_csv(sam_path, skiprows=196, header=None, index_col=None, sep='\t', usecols=[0,1,2,3,4,5])

# Some stats
total = len(sam_file)
unaligned = len(sam_file[sam_file[2] == '*'])
sam_file = sam_file[sam_file[2] != '*']
sam_file.index = range(len(sam_file))
aligned = len(sam_file)
sam_file[3] = sam_file[3].astype('int')
sam_file[6] = sam_file[0].str.split('|', expand=True)[2].str.split('_', expand=True)[0]
#print(sam.head())


# Caculate some statistics
# Perform gene annotation check
extended_sam = check_seq_annotations(sam_file, annot_grp=annot_copy_grp)

# Find all matching annotations
matched_sam = len(extended_sam[extended_sam['gene_name_sam'] == extended_sam['gene_name']])
#print("Total number of correct sgRNA asignments: {}".format(len(matched_sam)))

# Find all mismatched annotation
mismatched_sam = len(extended_sam[extended_sam['gene_name_sam'] != extended_sam['gene_name']])
#print("Total number of in-correct sgRNA asignments: {}".format(len(mismatched_sam)))

# Find all without perfect matches
indel_seq = len(extended_sam[extended_sam['cigar_sam'] != '20M'])

# Outinfo
outinfo = {sam_path: {
                        'Number of total gRNA': total,
                        'Number of unaligned gRNA': unaligned,
                        'Number of total aligned gRNA': aligned,
                        'Number of gRNA with correct annotations': matched_sam,
                        'Number of gRNA with in-correct annotations': mismatched_sam,
                        'Number of matches with indels/mismatches': indel_seq
                        }
            }


# All the matched genes from extended_sam to be concatnated with TCGA rna exp data
final_sam_tcga = add_tcga_data(extended_sam, tcga_path=os.path.join("data/TCGA"))
final_sam_tcga.to_csv(os.path.join(outpath, sample_name+'_final_df.tsv'), header=True, sep='\t', index=None)
print("Output of the pipeline is saved in: {}".format(os.path.join(outpath, sample_name+'_final_df.tsv')))

with open(os.path.join(outpath, sample_name+'_summary.json'), 'w') as file_summary:
  json.dump(outinfo, file_summary, indent=2)
  file_summary.close()
