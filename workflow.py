'''
This workflow is for mapping fastq files to a fasta 
file with genes in to calculate coverage.
The reads are mapped with BWA and and filtered.
The coverage for each gene is calculated are compared 
for each IND.


Author: CMK
Date: may. 2018 



 ~ todo ~
* record timeline/core/memusage etc.
'''
import sys, os
from os.path import basename # More paths?
from gwf import Workflow
from workflow_templates import *
import json # For pretty printing trees
import subprocess # For running bash in python
from datetime import datetime # To set dates in titles
import re # to escape symbols in strings

#2019
import pandas as pd


# Define the root gwf object
gwf = Workflow()


def get_individuals(path):
    '''
    Takes a path leading to a sample_info.txt-file
    and returns a list of individuals
    '''
    with open(path, 'r') as file:
        return [line.split('\t')[0] for line in file][1:]
        #                        ^ We are only interested in the names of the individuals (first column)
        #                                             ^ We are only interested in the rows after the labels (labels are on row 0)


def get_filenames(individuals, path):
    '''
    Takes a list of individuals
    and returns a dictionary for pairs of genome-files
    '''
    dico_names = {} # dictionary for unique individuals
    for subject in individuals:  # for each subject in the sample_info.txt-file:
        dico_names[subject] = []     # create an empty list for each subject
        # file = working_dir+'../Names_files/' + subject + '.txt'    # path to the link-file (every subject has its own file)
        with open(path+subject+'.txt') as file:      # close afterwards
            lines = [line.strip() for line in file] # transfer the strippedlines to a list, so we can later pair them together. The stripping removes newlines present in the end of each line.
            for pair in zip(lines[0::2], lines[1::2]):      # zipping the strided lines, puts them in pairs [(1,2), (2,3), ...]
                dico_names[subject].append(pair)     # append the file pairs to the list in the dictionary for each subject
    return dico_names


# This is a list of batches. Each batch is i dictionary with various keys and their values.
# KEY                       VALUE
# title                     The name of the batch
# chromosome                The name of the chromosome being worked on
# rel_ac             The relative path to the reference (artificial chromosome)
# rel_sample_info_file      The relative path to the sample info file, containing a list of individuals
# rel_individual_genomes    The relative path to the directory containing a file for each subject, each file containing a list og filenames which are the haploid genome builds
# abs_genome_dir            The absolute path to the directory containing the files listed in the rel_individual_genomes directory's files.


"""
# chimp
batches = [{"title": "batch_x6_chimp",
    "chromosome": "x",
    "description": "korsel fredag d 30. apr med ar. chrom. for chimp",
    "rel_ac": "ac_chimp_x.fa",
    "rel_sample_info_file": "sample_info_chimp.txt",
    "rel_individual_genomes": "genomes/",
    "abs_genome_dir": "/home/cmkobel/MutationRates/faststorage/NEW_PIPELINE/TrimFASTQ/TrimmedFASTQ/"},
    
    {"title": "batch_x6.3_gorilla",
    "chromosome": "x",
    "description": "gorilla med chimps opsætning ændret",
    "rel_ac": "ac_gorilla_x.fa",
    "rel_sample_info_file": "sample_info_gorilla.txt",
    "rel_individual_genomes": "genomes/",
    "abs_genome_dir": "/home/cmkobel/MutationRates/faststorage/NEW_PIPELINE/TrimFASTQ/TrimmedFASTQ/"},
    
    {"title": "batch_y6.2_chimp",
    "chromosome": "y",
    "description": "ændret opsætning til chimp y",
    "rel_ac": "ac_chimp_y.fa",
    "rel_sample_info_file": "sample_info_chimp.txt",
    "rel_individual_genomes": "genomes/",
    "abs_genome_dir": "/home/cmkobel/MutationRates/faststorage/NEW_PIPELINE/TrimFASTQ/TrimmedFASTQ/"}]
"""    

batch = {
    "prefix": "X_L",
    
    "rel_ac": "pan_tro_3_X.fa",
    "rel_subjects": "subjects/subjects.tsv",
    "description": "pantro3 complete X chromosome. Starting out with Carolina only."
}

subjects = pd.read_csv(batch['rel_subjects'], delimiter = '\t')
print(subjects.columns)




"""
# Protocol:
for batch in batches: # Måske batch bare skal hedde b ?
    pass
    individuals = get_individuals(batch['rel_sample_info_file'])

    print(individuals)
    #dico_names = get_filenames(individuals, batch['rel_individual_genomes'])

"""
# 0: initialize
gwf.target_from_template(batch['prefix']+'_0_initialize', initialize(batch['prefix'],
                                                                     str(batch)))

# 1: index
gwf.target_from_template(batch['prefix'] + '_1_index', index_genome(batch['prefix'], 
                                                                    batch['rel_ac']))


for _row, row in subjects.iterrows(): # for each subject
    subject = row['given_name']
    print('\n' + subject)

    with open(f'subjects/runs/{row["given_name"]}.txt', 'r') as runs_file:
        raw = [str(i).strip() for i in runs_file]
    forward = [i for i in raw[0::2]]
    reverse = [i for i in raw[1::2]]

    BAM_files=[] # for collecting BAM-files for later.
    #for _run, pair in enumerate(dico_names[subject]):
    for _run, (pe_forward, pe_reverse) in enumerate(zip(forward, reverse)):
        # 2: Mapping
        print(' ' + str(pe_forward), pe_reverse, sep = " | ")
        gwf.target_from_template(batch['prefix'] + '_2_map_'+subject.replace('-','_') + str(_run),
                                 bwa_map_pe(batch['prefix'],
                                            batch['rel_ac'],
                                            row['path']+'/'+pe_forward,
                                            row['path']+'/'+pe_reverse,
                                            subject+str(_run)))
        BAM_files.append(subject+str(_run)+'_sort_dedup.bam') # collect the filenames for use in the merge point
    print(' ', subject, 'BAM_files:', BAM_files)

    # 3: Merging the bam files # merge the bam-files for each subject
    gwf.target_from_template(batch['prefix'] + '_3_Merge_BAMS_' + subject.replace('-', '_'),
                             merge_bams_new(title = batch['prefix'],
                                            subject = subject,
                                            infiles = BAM_files,
                                            outfile = subject + '_merged.bam',
                                            input = subject + str(_run) + '_sort_dedup.bam'))

    # 4: Filtering the reads    # Remove bad quality reads?
    gwf.target_from_template(batch['prefix'] + '_4_filter_bam'+subject.replace('-', '_'),
                             filter_bam_file(batch['prefix'],
                                             subject))

    # 5: Get coverage for each subject
    gwf.target_from_template(batch['prefix'] + '_5_get_coverage'+subject.replace('-', '_'),
                             get_coverage(batch['prefix'], 
                                          subject))

    # 6: Calculate CNV
    gwf.target_from_template(batch['prefix'] + '_6_calc_cnv' + subject.replace('-', '_'), get_cnv(batch['prefix'],
                                                                                                  subject,
                                                                                                  'X'))
    #break # debug for only first subject (Carolina)

