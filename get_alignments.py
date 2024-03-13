import os
import get_lists

def sanitycheck(df, size, alignment_title):
    """Sanity check of the size of alignment files, number of rows should be same as
    the corresponding sub-list file."""
    if len(df.axes[0]) != size:
        print(alignment_title, 'did not pass sanitycheck, the size needed is:', size,
              ', but the exist size is:', len(df.axes[0]))

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
        pdb_list, sublist_path = get_lists.generate_sub_lists(pdb_path)
        for pdb_file in pdb_list[0]:
            list_file = 'list_' + pdb_file[:-3] + 'txt'
            align_file = 'align_' + pdb_file[:-3] + 'txt'
            align_file_path = os.path.join(align_folder_path, align_file)
            if os.path.exists(align_file_path):
                rm_cmd = 'rm ' + align_file_path
                os.system(rm_cmd)
            usalign_cmd = usalign_path + ' ' + os.path.join(pdb_path, pdb_file)  + \
                                  ' -dir2 ' + pdb_path + ' ' + \
                                  os.path.join(sublist_path, list_file) + ' -outfmt 2 >> ' + \
                                  align_file_path
            os.system(usalign_cmd)