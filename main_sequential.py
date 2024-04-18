"""this is an example to run on 192.168.66.203"""
"""note: the protein IDs should not contain '(' or ')'"""
import get_lists
import get_alignments
import clustering.get_clusters
pdb_folder_path='/dellfsqd2/ST_OCEAN/USER/wangdantong/python_toolbox_test/stacpro/pdb_files'
usalign_path='/dellfsqd2/ST_OCEAN/USER/wangdantong/toolboxes/usalign/USalign/USalign'
pdb_list, pdb_list_path = get_lists.get_lists(pdb_folder_path)
align_all_path = get_alignments.compute_similarity(pdb_folder_path, usalign_path)
clustering.get_clusters.get_tree_file(align_all_path, plot=1)
