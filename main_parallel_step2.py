"""this is an example to run on 27.18.114.42"""
"""note: the protein IDs should not contain '(' or ')'"""
import sys
import get_lists
import get_alignments
import clustering.get_clusters
import pandas as pd
pdb_folder_path='/home/share/huadjyin/home/fanguangyi/wangdantong/projects/PETs/predictions_2928/pdbs'
usalign_path='/home/share/huadjyin/home/fanguangyi/wangdantong/bin/USalign/USalign'
pdb_list_path = '/home/share/huadjyin/home/fanguangyi/wangdantong/projects/PETs/predictions_2928/pdb_list.txt'
sublist_path = '/home/share/huadjyin/home/fanguangyi/wangdantong/projects/PETs/predictions_2928/sublists'
pdb_list = pd.read_csv(pdb_list_path, sep='\t', header=None)
parallel = 1
align_all_path = '/home/share/huadjyin/home/fanguangyi/wangdantong/projects/PETs/predictions_2928/alignments/alignment_all.txt'
pdb_list, pdb_list_path = get_lists.get_lists(pdb_folder_path, parallel=1, par_num=30)
align_all_path = get_alignments.compute_similarity(pdb_folder_path, usalign_path, parallel=parallel,
                                                   par_index=int(sys.argv[1]), par_num=30, pdb_list=pdb_list, sublist_path=sublist_path)
clustering.get_clusters.get_tree_file(align_all_path)
