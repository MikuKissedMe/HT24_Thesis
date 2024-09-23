import sys
import os
import time

from pathlib import Path
from Bio import Phylo



script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
print("Current working directory:", os.getcwd())

def PIC(tree):
    # The PIC algorithm should only collapse leafs for n - 1 times, i.e. (total_leaf - 1) times
    total_leaf = len(tree.get_terminals())

    while total_leaf > 1:
        print(total_leaf)
        # Uppdate the list of leafs since the tree will constantly be modified
        leafs = tree.get_terminals()
        count = len(leafs)

        # Using reversed index seems to make it faster
        for i in reversed(range(1, count)):
            for j in reversed(range(0, i)):
                LCA = tree.common_ancestor(leafs[i],leafs[j])
                collapsed = collapse_leaf(leafs[i],leafs[j],LCA)
                if collapsed:
                    total_leaf -= 1
                    break
                else:
                    pass
            else:
                continue
            break
    print('done')
    print(tree)
    print(tree.get_terminals())

def add_length(leaf_i,leaf_j):
    ''' Currently only computing branch length'''
    if (float(leaf_i.branch_length) == 0) and (float(leaf_j.branch_length) == 0):
        return float(0)
    else:
        return (float(leaf_i.branch_length) * float(leaf_j.branch_length))/(float(leaf_i.branch_length) + float(leaf_j.branch_length))


def collapse_leaf(leaf_i,leaf_j,LCA):
    ''' Returns boolean that shows if leafs are collapsed'''
    if leaf_i in LCA and leaf_j in LCA:
        # If the last common ancestor are also direkt ancestor to leafs

        if LCA.branch_length == None:
            # If the leafs are direkt parents of root, and since root does not have a branch length
            leaf_i.branch_length = str(float(leaf_i.branch_length) + add_length(leaf_i,leaf_j))
            tree.collapse(leaf_j)
            return True

        else:
            # Regular leafs
            LCA.branch_length = str(float(LCA.branch_length) + add_length(leaf_i,leaf_j))
            tree.collapse(leaf_i)
            tree.collapse(leaf_j)
            return True
    else:
        # The last common ancestor are not direkt ancestor for some leaf
        return False



tree = Phylo.read("PS00027.treefile", "newick")
#tree = Phylo.read("PS00673_full.phy.treefile", "newick")
print(tree)


start_time = time.time()

PIC(tree)

print("--- %s seconds ---" % (time.time() - start_time))
