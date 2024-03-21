import os
import get_lists
import pandas as pd


def sanitycheck(df, size, alignment_title):
    """Sanity check of the size of alignment files, number of rows should be same as
    the corresponding sub-list file."""
    if len(df.axes[0]) != size:
        print(alignment_title, 'did not pass sanitycheck, the size needed is:', size,
              ', but the exist size is:', len(df.axes[0]))


def clean_pro_name_in_align(align_file_path, pdb_file, size_align):
    """Rename first two columns. The previous first column values are pdb file names with path. We keep the file names
    only and remove all the other characters"""
    df_align = pd.read_csv(align_file_path, sep='\t')
    # remove '/' in second column
    df_align['PDBchain2'] = df_align['PDBchain2'].str.replace('/', '')
    df_align['PDBchain2'] = df_align['PDBchain2'].str.replace('.pdb:A', '')
    df_align.rename(columns={'#PDBchain1': 'PDBchain1'}, inplace=True)
    df_align['PDBchain1'] = pdb_file[:-4]
    sanitycheck(df_align, size_align, align_file_path)
    df_align.to_csv(align_file_path, index=None, sep='\t')


def run_usalign(pdb_path, usalign_path, parallel):
    """Run us-align to compute the pair-wise similarity of structures"""
    prjfolder = os.path.dirname(pdb_path)
    # create folder for alignments
    align_folder = 'alignments'
    align_folder_path = os.path.join(prjfolder, align_folder)
    if not os.path.exists(align_folder_path):
        # if the folder for sub-lists is not present
        # then create it.
        os.makedirs(align_folder_path)
        print('Folder for alignments does not exist, generated!')
    if parallel:
        # if compute similarity parallelly, generate sub-lists
        pdb_list, sublist_path = get_lists.generate_sub_lists(pdb_path)
        # for every pdb file, run similarity computation according to each pdb files in the pdb folder
        size_align = pdb_list.size - 1
        for pdb_file in pdb_list[0]:
            # sub-list file name
            list_file = 'list_' + pdb_file[:-3] + 'txt'
            # sub-alignment file name
            align_file = 'align_' + pdb_file[:-3] + 'txt'
            align_file_path = os.path.join(align_folder_path, align_file)
            # if the alignment file exit already, remove it before write into it. TODO: if it exists, one can check it
            #  first, then decide if it should be removed or skip it to reduce computation time
            if os.path.exists(align_file_path):
                rm_cmd = 'rm ' + align_file_path
                os.system(rm_cmd)
            # do us-align
            usalign_cmd = usalign_path + ' ' + os.path.join(pdb_path, pdb_file)  + \
                                  ' -dir2 ' + pdb_path + ' ' + \
                                  os.path.join(sublist_path, list_file) + ' -outfmt 2 >> ' + \
                                  align_file_path
            os.system(usalign_cmd)
            # clean protein names and do sanity check of the alignment size
            clean_pro_name_in_align(align_file_path, pdb_file, size_align)
            size_align -= 1
    # if not parallel
    else:
        pdb_list, pdb_list_path = get_lists.get_pdblist_all(pdb_path)
        align_file = 'alignment_all.txt'
        align_file_path = os.path.join(align_folder_path, align_file)
        usalign_cmd = usalign_path + ' -dir ' + pdb_path + ' ' + pdb_list_path  + \
                      ' -outfmt 2 >> ' + align_file_path
        os.system(usalign_cmd)
    return pdb_list, align_folder_path, align_file_path


def cat_align(pdb_path, usalign_path):
    """If the alignments are computed parallelly, concatenate them."""
    # get pdb_list and folder path of sub-alignments
    pdb_list, align_folder_path, _ = run_usalign(pdb_path, usalign_path, 1)
    # for each sub-alignments, concatenate it to the full alignments
    for pdb_file in pdb_list[0]:
        try:
            alignment2 = 'align_' + pdb_file[:-3] + 'txt'
            path2 = os.path.join(align_folder_path, alignment2)
            df_2 = pd.read_csv(path2, sep='\t')
            df_align_all = pd.concat([df_align_all, df_2])
        except:
            alignment1 = 'align_' + pdb_file[:-3] + 'txt'
            path1 = os.path.join(align_folder_path, alignment1)
            df_align_all = pd.read_csv(path1, sep='\t')
    align_all_path = os.path.join(align_folder_path, 'alignment_all.txt')
    df_align_all.to_csv(align_all_path, index=None, sep='\t')
    return align_all_path


def compute_similarity(pdb_path, usalign_path, parallel):
    """This is a function to generate the pair-wise similarity matrix"""
    if parallel:
        align_all_path = cat_align(pdb_path, usalign_path)
    else:
        _, _, align_all_path = run_usalign(pdb_path, usalign_path, parallel)
    print('Pair-wise similarity computed, saved in: ', align_all_path)
