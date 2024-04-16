import copy
import numpy as np
import pandas as pd
import os


def get_lists(path,
              pdb_list_path=None,
              list_folder_path=None,
              parallel=0, par_num=5):
    if parallel:
        pdb_list, pdb_list_path = generate_sub_lists(path, list_folder_path=list_folder_path, par_num=par_num)
    else:
        pdb_list, pdb_list_path = get_pdblist_all(path, pdb_list_path=pdb_list_path)
    return pdb_list, pdb_list_path

def get_pdblist_all(path, pdb_list_path = None):
    """get a list of all .pdb files in the provided folder.
    Parameters:
    ----------
    path: string
        Path of the folder containing all .pdb files for alignment and clustering.
    pdb_list_path: sting
        Path where to save the generated list of all .pdb files. If provided, it should be end with
        a name of .txt file (e.g. PATH/list_pdb.txt)
    Returns:
    ----------
    pdb_list: pandas dataframe
        A list of all .pdb file names.
    pdb_list_path: string
        Path of the generated .txt file containing all .pdb file names"""
    # initialization of pdb list
    pdb_list = []
    # loop over all files in the folder, if it is .pdb file, then add it to the list
    for file in os.listdir(path):
        if file.endswith(".pdb"):
            pdb_list.append(file)
    pdb_list = pd.DataFrame(pdb_list)
    # save the list in the same path as the file contains all .pdb files
    prj_path = os.path.dirname(path)
    if pdb_list_path is None:
        pdb_list_path = os.path.join(prj_path, 'pdb_list.txt')
    pdb_list.to_csv(pdb_list_path, header=None, index=None)
    return pdb_list, pdb_list_path


def generate_sub_lists(path, list_folder_path=None, par_num=None):
    """generate sub-lists of .pdb files for running structure alignment in parallel.
    Parameters:
    ----------
    path: string
        Path of the folder containing all .pdb files for alignment and clustering.
    list_folder_path: sting
        Folder path where to save the generated sub-lists for computing protein similarity
        in parallel.
    Returns:
    ----------
    df_all: pandas dataframe
         A list of all .pdb file names.
    list_folder_path: string
        Folder path where to save the generated sub-lists for computing protein similarity
        in parallel, if not provided, this is in the same path as the pdb folder."""
    if list_folder_path is None:
        # get path for the sub-lists
        prjfolder = os.path.dirname(path)
        list_folder = 'sublists'
        list_folder_path = os.path.join(prjfolder, list_folder)
        # if the folder for sub-lists is not present
        if not os.path.exists(list_folder_path):
            # then create it.
            os.makedirs(list_folder_path)
            print('Folder for sub-lists of pdb files does not exist, created!')
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
    sub_pdb_list_size = int(np.floor(df_all.size/par_num))
    for ind_par_num in range(par_num):
        if ind_par_num == par_num - 1:
            sub_df_all = df_all[ind_par_num * sub_pdb_list_size:]
        else:
            sub_df_all = df_all[ind_par_num * sub_pdb_list_size:ind_par_num * sub_pdb_list_size + sub_pdb_list_size]
        par_list = 'pdb_list' + str(ind_par_num) + '.txt'
        parlist_file_path = os.path.join(list_folder_path, par_list)
        sub_df_all.to_csv(parlist_file_path, header=None, index=None)
    return df_all, list_folder_path
