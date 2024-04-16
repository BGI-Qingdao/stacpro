"""this is an example to run on 192.168.66.203"""
"""note: the protein IDs should not contain '(' or ')'"""
import get_lists
import get_alignments
import clustering.get_clusters
#pdb_folder_path = '/dellfsqd2/ST_OCEAN/USER/wangdantong/python_toolbox_test/stacpro/pdb_files'
usalign_path = '/dellfsqd2/ST_OCEAN/USER/wangdantong/toolboxes/usalign/USalign/USalign'
pdb_folder_path = 'D:\\python_toolbox_test\\stacpro\\pdb_files'
parallel = 1
pdb_list, pdb_list_path = get_lists.get_lists(pdb_folder_path, parallel=1, par_num=3)
align_all_path = get_alignments.compute_similarity(pdb_folder_path, usalign_path, parallel=1,
                                                   par_index=0, par_num=3)
clustering.get_clusters.get_tree_file(align_all_path)

get_alignments.run_usalign(pdb_folder_path, usalign_path, 1, par_index=1,
                par_num=3)


align_folder_path = '/home/share/huadjyin/home/fanguangyi/wangdantong/projects/PETs/predictions_2928/alignments'
get_alignments.cat_align(pdb_folder_path, usalign_path, align_folder_path=align_folder_path,
                         pdb_list=pdb_list, sublist_path=sublist_path)


