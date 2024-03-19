import copy
import pandas as pd
import os
def get_pdblist_all(path):
    """get a list of all .pdb files in the provided folder."""
    pdb_list = []
    # loop over all files in the folder, if it is .pdb file, then add it to the list
    for file in os.listdir(path):
        if file.endswith(".pdb"):
            pdb_list.append(file)
    pdb_list = pd.DataFrame(pdb_list)
    # save the list in the same path as the file contains all .pdb files
    prj_path = os.path.dirname(path)
    pdb_list_path = os.path.join(prj_path, 'pdb_list.txt')
    pdb_list.to_csv(pdb_list_path, header=None, index=None)
    return pdb_list, pdb_list_path

def generate_sub_lists(path):
    """generate sub-lists of .pdb files for running structure alignment parallelly."""
    # get path for the sub-lists
    prjfolder = os.path.dirname(path)
    list_folder = 'sublists'
    list_folder_path = os.path.join(prjfolder, list_folder)
    if not os.path.exists(list_folder_path):
        # if the folder for sub-lists is not present
        # then create it.
        os.makedirs(list_folder_path)
        print('Folder for sub-lists of pdb files does not exist, generated!')
    # get the list of all .pdb files, this is one of the outputs
    df_all, _ = get_pdblist_all(path)
    # copy of the overall list
    df4loop = copy.deepcopy(df_all)
    # save the list, then pop one file name out, save again, get all sub-lists in the end
    for i_file in range(df_all.size):
        txt_name = 'list_' + str(df4loop[0].iloc[0][:-4]) + '.txt'
        txt_path = os.path.join(list_folder_path, txt_name)
        df4loop[0].pop(i_file)
        df4loop[0].to_csv(txt_path, header=None, index=None)
    return df_all, list_folder_path
