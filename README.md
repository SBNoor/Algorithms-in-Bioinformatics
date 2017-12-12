# Algorithms-in-Bioinformatics

The codes global_affine.py and global_linear.py is the implementation of global pairwise alignment using affine gap cost and linear gap cost respectively. The running time for both codes is O(n^2).

The code sp_exact.py and MSA.py performs alignment of multiple sequences and not just two sequences unlike global_affine.py and global_linear.py. Moreover, sp_exact.py can only work on three sequences. The code will not work on more than three sequences. And MSA.py is a more general implementation and works on as many sequences as the user wants. The downside of MSA.py is that it will give an approximate score since it's a heuristics. But that score will always be 4/3 of the optimal (true) score.

The code rf_dist.py and RF_dynamic.py computes the RF distance between two unrooted evolutionary trees over the same set of species. rf_dist.py is the implementation of Day's algorithm and RF-dynamic.py makes use of quartet distance.

The neighborjoining.py code uses neighbor joining algorithm to construct an unrooted evolutionary tree. It produces a tree in newick format at the end.
