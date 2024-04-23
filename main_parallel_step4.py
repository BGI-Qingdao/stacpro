"""this is an example to run on 27.18.114.42"""
"""note: the protein IDs should not contain '(' or ')'"""
"""this scrip generates the nwk file for tree plot on itol website"""
import clustering.get_clusters
######### this are parameters to be modified
align_all_path = 'D:\\python_toolbox_test\\stacpro\\alignments\\alignment_all.txt'
# this is to define how many node you want to have in one cluster
node_number_upward = 3

######### this lines you do not need to touch
labels, path_tree_folder = clustering.get_clusters.get_tree_file(align_all_path)
df = clustering.get_clusters.clustering_upward(labels, node_number_upward, path_tree_folder)