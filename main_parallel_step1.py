"""this is an example to run on192.168.66.203"""
"""note: the protein IDs should not contain '(' or ')'"""
import os
import get_lists
pdb_folder_path = '/dellfsqd2/ST_OCEAN/USER/wangdantong/python_toolbox_test/stacpro/pdb_files'
usalign_path = '/dellfsqd2/ST_OCEAN/USER/wangdantong/toolboxes/usalign/USalign/USalign'
parallel = 1
# define how many jobs you want to divide the alignment computation into.
par_num = 30
pdb_list, pdb_list_path = get_lists.get_lists(pdb_folder_path, parallel=parallel, par_num=par_num)
print('Please use this path as the "sublist_path" input for the next step:')
print(pdb_list_path)
pdb_txt_path = os.path.join(pdb_list_path, 'pdb_list.txt')
print('Please use this path as the "pdb_list_path" input for the next step:')
print(pdb_txt_path)
