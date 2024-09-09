import sys
from pathlib import Path

path_root = Path(__file__).parents[0]
sys.path.append(str(path_root))
#print(path_root)

from dependency.pyseqlogo_master.pyseqlogo.pyseqlogo import draw_logo, setup_axis 



# Read the tree from a Newick file
tree = Phylo.read("example.nwk", "newick")

# Read the alignment file (assuming it's in FASTA format)
alignment = AlignIO.read("alignment.fasta", "fasta")




