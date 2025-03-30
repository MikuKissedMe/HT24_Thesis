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
def export_scores_to_file(scores_str,file_name):
    try:
        with open(file_name, "x") as f:
            f.write(scores_str)
    except FileExistsError:
        print("Already exists.")

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
    for i in range(0,config.seq_length):
        scores = []
        for c in config.existing_characters:
            scores.append((str(c),matrix[config.charaters.index(c)][i]))
        all_scores.append(scores)
    export_scores_to_file(str(all_scores),'with_PIC.txt')
    return all_scores

def convert_count_to_scores(matrix):
    '''convert the count dictionary of original data (sequnces) to the format used by pyseqlogo '''
    all_scores_count = []
    for i in range(0,config.seq_length):
        position_scores = []
        for key, value in matrix.items():
            position_scores.append((str(key),float(value[i]/config.terminals)))
        all_scores_count.append(position_scores)
    export_scores_to_file(str(all_scores_count),'without_PIC.txt')
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
    fig, axarr = draw_logo(convert_count_to_scores(config.seq_counter),colorscheme='hydrophobicity',draw_axis=True)
    fig.tight_layout()
    fig.savefig("without_PIC_bits.png")
    print("without_PIC_bits.png saved")

    print("---Done in %s seconds ---" % round(time.time() - start_time))

    ALL_SCORES1 = [[('C', 0.02247014831444764),
              ('T', 0.057903843733384308),
              ('A', 0.10370837683591219),
              ('G', 0.24803586793255664)],
             [('T', 0.046608227674354567),
              ('G', 0.048827667087419063),
              ('A', 0.084338697696451109),
              ('C', 0.92994511407402669)],
             [('G', 0.0),
              ('T', 0.011098351287382456),
              ('A', 0.022196702574764911),
              ('C', 1.8164301607015951)],
             [('C', 0.020803153636453006),
              ('T', 0.078011826136698756),
              ('G', 0.11268374886412044),
              ('A', 0.65529933954826969)],
             [('T', 0.017393530660176126),
              ('A', 0.030438678655308221),
              ('G', 0.22611589858228964),
              ('C', 0.45078233627623127)],
             [('G', 0.022364103549245576),
              ('A', 0.043412671595594352),
              ('T', 0.097349627214363091),
              ('C', 0.1657574733649966)],
             [('C', 0.03264675899941203),
              ('T', 0.045203204768416654),
              ('G', 0.082872542075430544),
              ('A', 1.0949220710572034)],
             [('C', 0.0),
              ('T', 0.0076232429756614498),
              ('A', 0.011434864463492175),
              ('G', 1.8867526364762088)],
             [('C', 0.0018955903000026028),
              ('T', 0.0094779515000130137),
              ('A', 0.35637097640048931),
              ('G', 0.58005063180079641)],
             [('A', 0.01594690817903021),
              ('C', 0.017541598996933229),
              ('T', 0.2774762023151256),
              ('G', 0.48638069946042134)],
             [('A', 0.003770051401807444),
              ('C', 0.0075401028036148881),
              ('T', 0.011310154205422331),
              ('G', 1.8624053924928772)],
             [('C', 0.036479877757360731),
              ('A', 0.041691288865555121),
              ('T', 0.072959755514721461),
              ('G', 1.1517218549109602)],
             [('G', 0.011831087684038642),
              ('T', 0.068620308567424126),
              ('A', 0.10174735408273231),
              ('C', 1.0009100180696691)],
             [('C', 0.015871770937774379),
              ('T', 0.018757547471915176),
              ('A', 0.32176408355669878),
              ('G', 0.36505073156881074)],
             [('A', 0.022798100897300954),
              ('T', 0.024064662058262118),
              ('G', 0.24571286522646588),
              ('C', 0.34070495229855319)]]
    
    plt.rcParams['figure.dpi'] = 300
    fig, axarr = draw_logo(ALL_SCORES1,yaxis='bits',draw_axis=True)
    fig.tight_layout()
    fig.savefig('seq_ex.png')
    print('seq_ex.png saved')

if __name__ == "__main__":
    try:
        main()
    except argparse.ArgumentError as e:
        print(f"Error: {e}")
    except FileNotFoundError as e:
        print(f"Error: {e}")






