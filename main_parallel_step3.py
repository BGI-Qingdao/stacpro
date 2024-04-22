"""this is an example to run on 27.18.114.42"""
"""note: the protein IDs should not contain '(' or ')'"""
"""this scrip is to be submitted in parallel to compute the similarity between all proteins"""
import clustering.get_clusters
align_all_path = '/home/share/huadjyin/home/fanguangyi/wangdantong/projects/PETs/predictions_2928/alignments/alignment_all.txt'
clustering.get_clusters.get_tree_file(align_all_path)