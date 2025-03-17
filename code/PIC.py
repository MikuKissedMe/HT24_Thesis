import sys
import os
import time
import math
import argparse
import warnings
import config

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from pathlib import Path
from Bio import Phylo
from Bio import SeqIO
from pyseqlogo.pyseqlogo import draw_logo, setup_axis

#Suppress runtime warnings caused by /pyseqlogo/format_utils.py
#warnings.filterwarnings("ignore", category=RuntimeWarning)


def convert_matrix_to_scores(matrix):
    '''convert the seq_counter matrix to the format used by pyseqlogo '''
    for c in config.charaters:
        if sum(config.seq_counter[c]) == 0:
            config.seq_counter.pop(c)
        elif c == 'X': # removes counts for X in sequnces since it indicates any/unknown
            config.seq_counter.pop(c)
        else:
            config.existing_characters.append(c)

    all_scores = []
    for i in range(0,24):
        scores = []
        for c in config.existing_characters:
            scores.append((str(c),matrix[config.charaters.index(c)][i]))
        all_scores.append(scores)
    return all_scores

def convert_count_to_scores(matrix):
    '''convert the count dictionary of original data (sequnces) to the format used by pyseqlogo '''
    all_scores_count = []
    for i in range(0,config.seq_length):
        position_scores = []
        for key, value in matrix.items():
            position_scores.append((str(key),float(value[i]/config.terminals)))
        all_scores_count.append(position_scores)

    return all_scores_count

def set_seq(leaf_i,leaf_j):
    if (float(leaf_i.branch_length) == 0) and (float(leaf_j.branch_length) == 0):
        return config.matrix
    else:
        seq_matrix = np.add((float(leaf_i.branch_length) / (float(leaf_i.branch_length) + float(leaf_j.branch_length))) * config.seq_dict[leaf_j],
        (float(leaf_j.branch_length) / (float(leaf_i.branch_length) + float(leaf_j.branch_length))) * config.seq_dict[leaf_i])
        return seq_matrix


def add_length(leaf_i,leaf_j):
    if (float(leaf_i.branch_length) == 0) and (float(leaf_j.branch_length) == 0):
        return float(0)
    else:
        return (float(leaf_i.branch_length) * float(leaf_j.branch_length))/(float(leaf_i.branch_length) + float(leaf_j.branch_length))


def unbifrucating(childs):
    ''' For case where a node has more than 2 childs'''
    branch_length_temp = []
    seq_matrix_temp = []
    for child in childs:
        branch_length_temp.append(float(child.branch_length))
        seq_matrix_temp.append(config.seq_dict[child])

    while (len(branch_length_temp)> 1 and len(seq_matrix_temp)> 1):
        leaf_i_branch_length = branch_length_temp.pop(0)
        leaf_j_branch_length = branch_length_temp.pop(0)
        leaf_i_seq_matrix = seq_matrix_temp.pop(0)
        leaf_j_seq_matrix = seq_matrix_temp.pop(0)

        if (leaf_i_branch_length == 0) and (leaf_j_branch_length == 0):
                branch_length_temp.append(0)
                seq_matrix_temp.append(config.matrix)
        else:
            branch_length_temp.append(((leaf_i_branch_length * leaf_j_branch_length)/(leaf_i_branch_length + leaf_j_branch_length)))
            seq_matrix_temp.append(np.add(
                (leaf_i_branch_length / (leaf_i_branch_length + leaf_j_branch_length)) * leaf_i_seq_matrix,
                (leaf_j_branch_length / (leaf_i_branch_length + leaf_j_branch_length)) * leaf_j_seq_matrix))
    return branch_length_temp[0] , seq_matrix_temp[0]


def PIC_postorder (tree):
    for child in tree.clade:
        traverse_postorder(child)
    
    result , seq_matrix = unbifrucating(tree.clade)
    ALL_SCORES1 = convert_matrix_to_scores(seq_matrix)

    '''
    plt.rcParams['figure.dpi'] = 300
    fig, axarr = draw_logo(ALL_SCORES1,seq_type='aa', yaxis='probability', colorscheme='hydrophobicity')
    fig.tight_layout()
    fig.show()
    fig.savefig("with_PIC.png")
    print("with_PIC.png saved")
    '''

    plt.rcParams['figure.dpi'] = 300
    fig, axarr = draw_logo(ALL_SCORES1,data_type='bits',seq_type=config.seq_type, yaxis='bits',colorscheme='hydrophobicity',draw_axis=True)
    fig.tight_layout()
    fig.show()
    fig.savefig("with_PIC_bits.png")
    print("with_PIC_bits.png saved")


def traverse_postorder(clade):
    if len(clade) == 0: #only tips of the tree will have length 0
        clade.seq = str(config.updated_dict[clade.name].seq) #store str(sequnces) as an artribute for the clade object.

        seq_matrix = config.matrix.copy()
        for i in range(0,len(clade.seq)):
            config.seq_counter[clade.seq[i].upper()][i] += 1 #count and store the frequncy of each character

            character_index = config.charaters.index(clade.seq[i].upper())
            seq_matrix[character_index,i] = float(1)

        config.seq_dict[clade] = seq_matrix #stores the matrix to dictionary seq

    if len(clade) > 0:
        for child in clade:
            traverse_postorder(child)
        if len(clade) == 2:
            clade.branch_length = float(clade.branch_length) + add_length(clade[0],clade[1])
            config.seq_dict[clade] = set_seq(clade[0],clade[1])
        if len(clade) > 2:
            result , seq_matrix = unbifrucating(clade)
            clade.branch_length = float(clade.branch_length) + result
            config.seq_dict[clade] = seq_matrix


def parse_dict(filename,filetype = 'fasta'):
    '''removes suffix from sequnce names'''
    record_dict = SeqIO.to_dict(SeqIO.parse(filename, filetype))
    return remove_position_from_key(record_dict)


def remove_position_from_key(record_dict):
    '''removes suffix, also determine length and type of sequnces.  ex. removes /60-74 from SPLA_STAA9/60-74 '''
    updated_dict = {}
    for key, value in record_dict.items():
        if len(value) != config.seq_length:
            if config.seq_length == 0:
                config.seq_length = len(value)
            else:
                raise Exception("Sequnce length is inconsistent")

        if set(value.upper()) not in [{'A','C','T','G'}, {'A','C','T','G','X'}]:
            config.seq_type = 'aa'

        new_key = key.split('/')[0]
        updated_dict [new_key] = value
    return updated_dict


def find_clades(clade, condition):
    """Find clades matching a condition, used in testing."""
    matches = []
    if condition(clade):
        matches.append(clade)
    for sub_clade in clade.clades:
        matches.extend(find_clades(sub_clade, condition))
    return matches


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Apply PIC on tree file and corresbond Name-sequnce file, outputs 2 .png file(sequnce logo with and with out applying PIC)."
    )
    parser.add_argument(
        "input_tree", 
        type=str, 
        help="Path to the .tree file (newick format)"
    )
    parser.add_argument(
        "input_fa", 
        type=str, 
        help="Path to the .fa file (fasta format), assumes sequnce name in the format of: ETA_STAAU/96-110"
    )
    args = parser.parse_args()

    # Validate file paths
    file1 = Path(args.input_tree)
    file2 = Path(args.input_fa)

    if not file1.is_file():
        raise FileNotFoundError(f".tree {file1} does not exist.")
    if not file2.is_file():
        raise FileNotFoundError(f".fa {file2} does not exist.")

    print(f"Processing: {file1} and {file2}")

    tree = Phylo.read(file1, "newick")
    #tree = Phylo.read("PS00673_full.phy.treefile", "newick")
    #print(tree)
    #print(tree.count_terminals())
    config.terminals = tree.count_terminals()
    
    #config.py to store global variables
    config.updated_dict = parse_dict(file2, "fasta")
    print('sequnce type = ' + config.seq_type)
    print('sequnce length = ' + str(config.seq_length))
    #print(config.updated_dict['ABDA_AEDAE'].seq)

    config.seq_dict = {} # dictionary for storing sequnce matrix
    config.seq_counter = {} # dictionary storing counts for seqlogo
    config.matrix = np.zeros((26, config.seq_length), dtype=float) # a default matrix for each individual leaf/node, to store character counts.
    config.charaters = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    config.existing_characters = []

    # initialize the seq_counter dictionary
    for i in range(0,len(config.charaters)):
        config.seq_counter[config.charaters[i]] = np.zeros((config.seq_length,), dtype=int).tolist()

    start_time = time.time()
    PIC_postorder(tree)
    
    #print(tree)
    #print(seq_counter)
    #counts = {'A' : [3,4,5,6], 'C': [2,3,1,1], 'T': [2,1,3,1], 'G': [3,2,1,2]}
    #print(config.seq_counter)
    '''
    fig, axarr = draw_logo(config.seq_counter, data_type='counts', seq_type= config.seq_type , yaxis='probability', colorscheme='hydrophobicity')
    fig.tight_layout()
    fig.show()
    fig.savefig("without_PIC.png")
    print("without_PIC.png saved")'
    '''

    plt.rcParams['figure.dpi'] = 300
    fig, axarr = draw_logo(convert_count_to_scores(config.seq_counter),data_type='bits',seq_type=config.seq_type, yaxis='bits',colorscheme='hydrophobicity',draw_axis=True)
    fig.tight_layout()
    fig.savefig("without_PIC_bits_bits.png")
    print("without_PIC_bits.png saved")

    print("---Done in %s seconds ---" % round(time.time() - start_time))


if __name__ == "__main__":
    try:
        main()
    except argparse.ArgumentError as e:
        print(f"Error: {e}")
    except FileNotFoundError as e:
        print(f"Error: {e}")






