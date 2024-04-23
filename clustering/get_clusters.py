import clustering.tree_functions as tree_functions
import pandas as pd
import os


def get_tree_file(path_similarity, method='nj', tm_score='average', path_tree_folder=None, tree_name=None, plot=0,
                  clust_num=3, clust_save_path=None):
    """This is a function to generate a .nwk file, which can be uploaded to "https://itol.embl.de/"
    to plot and edit the tree plot"""
    # reform the pair-wise similarity form to a similarity matrix
    mat, df_pro_uni = tree_functions.align2mat(path_similarity, tm_score=tm_score)
    # get the nearest proteins and update the similarity matrix
    if method == 'upgma':
        mat_new, row, col, min_dis = tree_functions.update_mat_upgma(mat, df_pro_uni.size - 1)
    elif method == 'nj':
        mat_new, row, col, min_dis1, min_dis2 = tree_functions.update_mat_nj(mat, df_pro_uni.size - 1)
    else:
        print('Wrong clustering method defined, or the ', method, ' is still not implemented, using '
                                                                  'NJ method instead')
        mat_new, row, col, min_dis1, min_dis2 = tree_functions.update_mat_nj(mat, df_pro_uni.size - 1)
    # get all protein IDs
    name_all = tree_functions.add_p(df_pro_uni)
    # initialize the list, which will then saved as the .nwk file
    name_distances = list(df_pro_uni)
    # get the nearest two protein IDs
    pro_pre = name_all[row]
    pro_post = name_all[col + 1]
    pro_pre_dis = name_distances[row]
    pro_post_dis = name_distances[col + 1]
    # update the list of all protein IDs by combining the two proteins and delet the later one
    name_all[row] = '(' + pro_pre + ',' + pro_post + ')'
    labels = [[[pro_pre, pro_post]]]
    if method == 'upgma':
        distances = [[min_dis]]
        name_distances[row] = '(' + pro_pre_dis + ':' + str(min_dis) + ',' + pro_post_dis + ':' \
                          + str(min_dis) + ')'
    else:
        distances1 = [[min_dis1]]
        distances2 = [[min_dis2]]
        name_distances[row] = '(' + pro_pre_dis + ':' + str(min_dis1) + ',' + pro_post_dis + ':' \
                              + str(min_dis2) + ')'
    # remove the later protein
    name_all.pop(col + 1)
    name_distances.pop(col + 1)
    # size of matrix, it was (n-1), and then removed one, thus (n-2)
    len_mat = df_pro_uni.size - 2
    # loop over the similarity matrix, remove one dimension in each step
    for iloop in range(len_mat):
        mat = mat_new
        # update matrix
        if method == 'upgma':
            mat_new, row, col, min_dis = tree_functions.update_mat_upgma(mat, len_mat - iloop)
        else:
            mat_new, row, col, min_dis1, min_dis2 = tree_functions.update_mat_nj(mat, len_mat - iloop)
        # first protein ID
        pro_pre = name_all[row]
        pro_pre_dis = name_distances[row]
        # second protein ID
        pro_post = name_all[col + 1]
        pro_post_dis = name_distances[col + 1]
        # if the first ID is a single protein ID
        if pro_pre[-1] != ')':
            # if the second protein ID is a single protein
            if pro_post[-1] != ')':
                if method == 'upgma':
                    name_all, labels, distances, name_distances = \
                        tree_functions.single_single_upgma(row, pro_pre, pro_post, labels, name_all,
                                                           min_dis, distances, name_distances, pro_pre_dis,
                                                           pro_post_dis)
                else:
                    name_all, labels, distances1, distances2, name_distances = \
                        tree_functions.single_single_nj(row, pro_pre, pro_post, labels, name_all,
                                                        min_dis1, min_dis2, distances1, distances2,
                                                        name_distances, pro_pre_dis, pro_post_dis)
            # if the second ID is a combination of multiple protein IDs
            else:
                if method == 'upgma':
                    name_all, labels, distances, name_distances = \
                        tree_functions.single_comp_upgma(row, pro_pre, pro_post, labels, name_all, min_dis,
                                                         distances, name_distances, pro_pre_dis,
                                                         pro_post_dis)
                else:
                    name_all, labels, distances1, distances2, name_distances = \
                        tree_functions.single_comp_nj(row, pro_pre, pro_post, labels, name_all,
                                                      min_dis1, min_dis2, distances1, distances2,
                                                      name_distances, pro_pre_dis, pro_post_dis)
        # if the first ID is a combination of multiple protein IDs
        else:
            # if the second ID is a single protein ID
            if pro_post[-1] != ')':
                if method == 'upgma':
                    name_all, labels, distances, name_distances = \
                        tree_functions.comp_single_upgma(row, pro_pre, pro_post, labels, name_all, min_dis,
                                                         distances, name_distances, pro_pre_dis,
                                                         pro_post_dis)
                else:
                    name_all, labels, distances1, distances2, name_distances = \
                        tree_functions.comp_single_nj(row, pro_pre, pro_post, labels, name_all,
                                                      min_dis1, min_dis2, distances1, distances2,
                                                      name_distances, pro_pre_dis, pro_post_dis)
            # if both IDs are combinations of multiple proteins
            else:
                if method == 'upgma':
                    name_all, labels, distances, name_distances = \
                        tree_functions.comp_comp_upgma(row, pro_pre, pro_post, labels, name_all, min_dis,
                                                       distances, name_distances, pro_pre_dis, pro_post_dis)
                else:
                    name_all, labels, distances1, distances2, name_distances = \
                        tree_functions.comp_comp_nj(row, pro_pre, pro_post, labels, name_all,
                                                    min_dis1, min_dis2, distances1, distances2,
                                                    name_distances, pro_pre_dis, pro_post_dis)
        name_all.pop(col + 1)
        name_distances.pop(col + 1)
    if tree_name is None:
        tree_name = 'tree_' + method + '.nwk'
    if path_tree_folder is None:
        alignfolder = os.path.dirname(path_similarity)
        path_tree_folder = os.path.join(alignfolder, 'trees')
        # if the folder for .nwk file is not present
        if not os.path.exists(path_tree_folder):
            # then create it.
            os.makedirs(path_tree_folder)

    save_tree_path = os.path.join(path_tree_folder, tree_name)
    file_tree = open(save_tree_path, 'w')
    file_tree.write(name_distances[0])
    file_tree.close()
    print('The nwk file for tree plot is saved in: ', save_tree_path)

    if plot:
        if method == 'upgma':
            fig_save_path = os.path.join(path_tree_folder, 'treeplot.png')
            fig, ax = tree_functions.plot_tree(clust_num, name_all, labels, distances)
            fig.savefig(fig_save_path)
        else:
            print('Sorry, we only support the printing of simple tree plot using UPGMA method, update'
                  'will come later, thank you for your understanding!')

    return labels, path_tree_folder


def clustering_upward(labels, node_number_upward, path_tree_folder):
    cluster_list = tree_functions.get_clusters(labels, node_number_upward)
    if cluster_list[-1] == []:
        cluster_list = cluster_list[:-1]
    list_cluster = []
    list_pro = []
    for ind_cluster, i_cluster in enumerate(cluster_list):
        num_cluster = ind_cluster
        for pro_id in i_cluster:
            list_cluster.append(num_cluster + 1)
            list_pro.append(pro_id)
    data_clusters = {'cluster_number': list_cluster, 'protein_ID': list_pro}
    df_clusters = pd.DataFrame(data_clusters)
    clusters_file = 'cluster_info.txt'
    clusters_save_path = os.path.join(path_tree_folder, clusters_file)
    df_clusters.to_csv(clusters_save_path, index=None, sep='\t')
    print('The information of clusters is saved in: ', clusters_save_path)
    return df_clusters
