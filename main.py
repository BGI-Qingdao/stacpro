"""this is an example to run on 192.168.66.203"""
"""note: the protein IDs should not contain '(' or ')'"""
import get_lists
import get_alignments
import clustering.get_clusters
pdb_path = '/dellfsqd2/ST_OCEAN/USER/wangdantong/python_toolbox_test/stacpro/pdb_files'
us_path = '/dellfsqd2/ST_OCEAN/USER/wangdantong/toolboxes/usalign/USalign/USalign'
parallel = 1
#get_alignments.run_usalign(pdb_path, us_path, parallel)
get_lists.generate_sub_lists(pdb_path)
align_all_path = get_alignments.compute_similarity(pdb_path, us_path, parallel)
clustering.get_clusters.get_tree_file(align_all_path)
