import os
import sys
import get_lists
import pandas as pd


def sanitycheck(df, size, alignment_title):
    """Sanity check of the size of alignment files, number of rows should be same as
    the corresponding sub-list file.
    Parameters:
    ----------
    df: pandas dataframe
        The alignment dataframe for sanity check.
    size: int
        The number of rows that the alignment dataframe should have.
    alignment_title: string
        The path of the checked alignment file. (e.g. PATH/align.txt)
    """
    if len(df.axes[0]) != size:
        print(alignment_title, 'did not pass sanitycheck, the size needed is:', size,
              ', but the exist size is:', len(df.axes[0]))


def clean_pro_name_in_align(align_file_path, size_align, pdb_file=None, parallel=1):
    """Rename first two columns.
    The previous first column values are pdb file names with path, while the second column values
    are full pdb file names with '.pdb'. We keep the protein IDs only and remove all the other
    characters.
    Parameters:
    ----------
    align_file_path: string
        The generated .txt file containing the TM-scores of all .pdb files used for alignment.
    size_align: int
        Number of rows the alignment dataframe should contain, used for sanity check.
    pdb_file: string
        ID of the protein used to compare with all the .pdb files in the corresponding list.
    parallel: bool
        If the alignment is computed in parallel."""
    df_align = pd.read_csv(align_file_path, sep='\t')
    # remove '/' in second column
    df_align['PDBchain2'] = df_align['PDBchain2'].str.replace('/', '')
    #df_align['PDBchain2'] = df_align['PDBchain2'].str.replace('.pdb:A', '')
    df_align['PDBchain2'] = df_align['PDBchain2'].str[:-6]
    df_align.rename(columns={'#PDBchain1': 'PDBchain1'}, inplace=True)
    if parallel:
        df_align['PDBchain1'] = pdb_file[:-4]
    else:
        df_align['PDBchain1'] = df_align['PDBchain1'].str.replace('/', '')
        #df_align['PDBchain1'] = df_align['PDBchain1'].str.replace('.pdb:A', '')
        df_align['PDBchain1'] = df_align['PDBchain1'].str[:-6]
    sanitycheck(df_align, size_align, align_file_path)
    df_align.to_csv(align_file_path, index=None, sep='\t')


def run_usalign(pdb_path, usalign_path, parallel, par_index=None,
                par_num=5, align_folder_path=None, align_file=None, pdb_list=None,
                sublist_path=None, pdb_list_path=None):
    """Run us-align to compute the pair-wise similarity of structures.
    Parameters:
    ----------
    pdb_path: string
        Path of the file containg all .pdb files used for alignment.
    usalign_path: string (!!IMPORTANT)
        Path of the US-align tool (e.g. PATH/USalign), which should be installed before computing the similarity matrix.
    parallel: bool
        If the alignment is computed in parallel.
    Returns:
    ----------
    pdb_list: pandas dataframe
        A list of all .pdb file names.
    align_folder_path: string
        Path of the folder containing the alignment files.
    align_file_path: string
        Path of the alignment file."""
    # check if the US-align tool exist, if not, stop here
    if not os.path.exists(usalign_path):
        print('Error: IMPORTANT!! US-align tool does not exist, please provide the correct path, '
              'or install before computing the similarity matrix! Program will stop here.')
        sys.exit()
    prjfolder = os.path.dirname(pdb_path)
    # create folder for alignments
    if align_folder_path is None:
        align_folder = 'alignments'
        align_folder_path = os.path.join(prjfolder, align_folder)
    if not os.path.exists(align_folder_path):
        # if the folder for sub-lists is not present
        # then create it.
        os.makedirs(align_folder_path)
        print('Folder for alignments does not exist, created!')
    if pdb_list is None or sublist_path is None:
        pdb_list, sublist_path = get_lists.get_lists(pdb_path, parallel=parallel, par_num=par_num)
    if parallel:
        sub_par_list_file = 'pdb_list' + str(par_index) + '.txt'
        sub_par_list_path = os.path.join(sublist_path, sub_par_list_file)
        df_sub_par = pd.read_csv(sub_par_list_path, sep='\t', header=None)
        # for every pdb file, run similarity computation according to each pdb files in the pdb folder
        for pdb_file in df_sub_par[0]:
            # sub-list file name
            list_file = 'list_' + pdb_file[:-3] + 'txt'
            sub_list_list_path = os.path.join(sublist_path, list_file)
            df_sub_list = pd.read_csv(sub_list_list_path, sep='\t', header=None)
            size_align = df_sub_list.size
            # sub-alignment file name
            align_file = 'align_' + pdb_file[:-3] + 'txt'
            align_file_path = os.path.join(align_folder_path, align_file)
            # if the alignment file exit already, remove it before write into it. TODO: if it exists, one can check it
            #  first, then decide if it should be removed or skip it to reduce computation time
            if os.path.exists(align_file_path):
                rm_cmd = 'rm ' + align_file_path
                os.system(rm_cmd)
            # do us-align
            usalign_cmd = usalign_path + ' ' + os.path.join(pdb_path, pdb_file) + \
                          ' -dir2 ' + pdb_path + ' ' + \
                          sub_list_list_path + ' -outfmt 2 >> ' + \
                          align_file_path
            os.system(usalign_cmd)
            # clean protein names and do sanity check of the alignment size
            clean_pro_name_in_align(align_file_path, size_align, pdb_file=pdb_file)
            size_align -= 1
    # if not parallel
    else:
        size_align = pdb_list.size * (pdb_list.size - 1) / 2
        if align_file is None:
            align_file = 'alignment_all.txt'
        align_file_path = os.path.join(align_folder_path, align_file)
        usalign_cmd = usalign_path + ' -dir ' + pdb_path + ' ' + sublist_path + \
                      ' -outfmt 2 >> ' + align_file_path
        os.system(usalign_cmd)
        # clean protein names and do sanity check of the alignment size
        clean_pro_name_in_align(align_file_path, size_align, parallel=parallel)
    return pdb_list, align_folder_path, align_file_path


def cat_align(pdb_path, usalign_path, par_index=None,
                par_num=5, align_folder_path=None, align_file=None, pdb_list=None, sublist_path=None):
    """If the alignments are computed parallelly, concatenate them.
    Parameters:
    ----------
    pdb_path: string
        Path of the file containg all .pdb files used for alignment.
    usalign_path: string
        Path of the US-align tool (e.g. PATH/USalign), which should be installed before computing the similarity matrix.
    Returns:
    ----------
    align_all_path: string
        Path of the final pair-wise similarity matrix.
    """
    # for each sub-alignments, concatenate it to the full alignments
    for pdb_file in pdb_list[0]:
        try:
            alignment2 = 'align_' + pdb_file[:-3] + 'txt'
            path2 = os.path.join(align_folder_path, alignment2)
            df_2 = pd.read_csv(path2, sep='\t')
            df_align_all = pd.concat([df_align_all, df_2])
        # initialize the full alignment dataframe
        except NameError:
            alignment1 = 'align_' + pdb_file[:-3] + 'txt'
            path1 = os.path.join(align_folder_path, alignment1)
            df_align_all = pd.read_csv(path1, sep='\t')
    if align_file is None:
        align_file = 'alignment_all.txt'
    align_all_path = os.path.join(align_folder_path, align_file)
    if os.path.exists(align_all_path):
        rm_cmd = 'rm ' + align_all_path
        os.system(rm_cmd)
    # save it
    df_align_all.to_csv(align_all_path, index=None, sep='\t')
    return align_all_path


def compute_similarity(pdb_path, usalign_path, parallel=0, par_index=None,
                par_num=5, align_folder_path=None, align_file=None, pdb_list=None,
                       sublist_path=None, pdb_list_path=None):
    """This is a function to generate the pair-wise similarity matrix.
    Parameters:
    ----------
    pdb_path:
        Path of the file containg all .pdb files used for alignment.
    usalign_path: string (!!IMPORTANT)
        Path of the US-align tool (e.g. PATH/USalign), which should be installed before computing the similarity matrix.
    parallel: bool
        If the alignment is computed in parallel.
    Returns:
    ----------
    align_all_path: string
        Path of the final pair-wise similarity matrix."""
    # if compute similarity in parallel, concatenate the individual alignments
    if parallel:
        print('Computing all alignments in parallel, use the function "cat_align" to concatenate them after all '
              'computations are finished.')
        if par_index is None:
            print('The index of parallel is not specified, computing for par_index=0.')
            par_index = 0
        # get pdb_list and folder path of sub-alignments
        pdb_list, align_folder_path, _ = run_usalign(pdb_path, usalign_path, 1, par_index=par_index,
                                                     par_num=par_num, align_folder_path=align_folder_path,
                                                     pdb_list=pdb_list, sublist_path=sublist_path)
    # in the similarity is not computed in parallel, simplly do us-align
    else:
        _, _, align_all_path = run_usalign(pdb_path, usalign_path, 0, par_index=None,
                                           par_num=None, align_folder_path=align_folder_path,
                                           align_file=align_file, pdb_list=pdb_list, pdb_list_path=pdb_list_path)
    print('Pair-wise similarity computed, saved in: ', align_all_path)
    return align_all_path
