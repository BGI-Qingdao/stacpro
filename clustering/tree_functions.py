import pandas as pd
import numpy as np
import copy
import random
import matplotlib.pyplot as plt


def align2mat(pairwise_sim_path, tm_score='average'):
    """This is a function to reform the pair-wise similarity to similarity matrix with the size of
    (n-1)*(n-1), with n the number of proteins.
    parameters:
    align_path: path for the pair-wise similarity file
    tm-score: values of TM-score to use for clustering (default = average, means average of TM1 and TM2);
              TM1: TM-score normalized using the first sequence;
              TM2: TM-score normalized using the second sequence;
              long: TM-score normalized using the longer sequence;
              short: TM-score normalized using the shorter sequence."""
    df = pd.read_csv(pairwise_sim_path, index_col=None, sep='\t')
    if tm_score == 'TM1':
        df_ave = df['TM1']
    elif tm_score == 'TM2':
        df_ave = df['TM2']
    elif tm_score == 'long':
        list_long = []
        for ind_tm in range(df['TM2'].size):
            min_tm = min(df['TM1'][ind_tm], df['TM2'][ind_tm])
            list_long.append(min_tm)
        df_ave = pd.Series(list_long)
    elif tm_score == 'short':
        list_short = []
        for ind_tm in range(df['TM2'].size):
            max_tm = max(df['TM1'][ind_tm], df['TM2'][ind_tm])
            list_short.append(max_tm)
        df_ave = pd.Series(list_short)
    else:
        df_ave = (df['TM1'] + df['TM2']) / 2
    df_normal = 1 - (df_ave - df_ave.min()) / (df_ave.max() - df_ave.min())
    df_pro_uni1 = df.PDBchain1.unique()
    df_pro_uni2 = df.PDBchain2.unique()
    df_pro_uni = np.append(df_pro_uni1, df_pro_uni2[-1])
    len_matrix = df_pro_uni1.size
    mat = np.ones([len_matrix, len_matrix]) + 0.0001
    count = 0
    for i in range(len_matrix):
        for j in range(i, len_matrix):
            # print(df_normal[count])
            mat[i, j] = df_normal[count]
            count = count + 1
    return mat, df_pro_uni


def update_mat_upgma(mat, len_matrix):
    min_dis = mat.min() / 2
    ind_all = mat.argmin()
    row = np.floor(ind_all / len_matrix)
    col = ind_all % len_matrix
    row = int(row)
    mat_new = np.ones([len_matrix - 1, len_matrix - 1]) + 0.0001
    for irow in range(len_matrix - 1):
        if irow == row:
            for icol in range(irow, len_matrix - 1):
                if icol >= col:
                    mat_new[irow, icol] = (mat[irow, icol + 1] + mat[col + 1, icol + 1]) / 2
                else:
                    mat_new[irow, icol] = (mat[irow, icol] + mat[icol + 1, col]) / 2
        elif irow > col:
            for icol in range(irow, len_matrix - 1):
                mat_new[irow, icol] = mat[irow + 1, icol + 1]
        else:
            for icol in range(irow, len_matrix - 1):
                if icol == row - 1:
                    mat_new[irow, icol] = (mat[irow, row - 1] + mat[irow, col]) / 2
                elif icol < col:
                    mat_new[irow, icol] = mat[irow, icol]
                else:
                    mat_new[irow, icol] = mat[irow, icol + 1]
    return mat_new, row, col, min_dis


def update_mat_nj(mat, len_matrix):
    if len_matrix == 1:
        ind_all = mat.argmin()
    else:
        r_mat = get_r_mat(mat, len_matrix)
        m_mat = get_m_mat(mat, r_mat, len_matrix)
        ind_all = m_mat.argmin()
    row = np.floor(ind_all / len_matrix)
    col = ind_all % len_matrix
    row = int(row)
    if len_matrix == 1:
        min_dis1 = mat[row, col] / 2
        min_dis2 = mat[row, col] - min_dis1
    else:
        min_dis1 = mat[row, col] / 2 + (r_mat[row] - r_mat[col + 1]) / (2 * (len_matrix - 1))
        min_dis2 = mat[row, col] - min_dis1
    mat_new = np.ones([len_matrix - 1, len_matrix - 1]) + 0.0001
    for irow in range(len_matrix - 1):
        if irow == row:
            for icol in range(irow, len_matrix - 1):
                if icol >= col:
                    mat_new[irow, icol] = (mat[irow, icol + 1] + mat[col + 1, icol + 1] - mat[row, col]) / 2
                else:
                    mat_new[irow, icol] = (mat[irow, icol] + mat[icol + 1, col] - mat[row, col]) / 2
        elif irow > col:
            for icol in range(irow, len_matrix - 1):
                mat_new[irow, icol] = mat[irow + 1, icol + 1]
        else:
            for icol in range(irow, len_matrix - 1):
                if icol == row - 1:
                    mat_new[irow, icol] = (mat[irow, row - 1] + mat[irow, col] - mat[row, col]) / 2
                elif icol < col:
                    mat_new[irow, icol] = mat[irow, icol]
                else:
                    mat_new[irow, icol] = mat[irow, icol + 1]
    if min_dis1 < 0:
        min_dis1 = 0.01
    if min_dis2 < 0:
        min_dis2 = 0.01
    return mat_new, row, col, min_dis1, min_dis2


def get_r_mat(mat, len_matrix):
    r_mat = []
    for ir in range(len_matrix + 1):
        if ir == 0:
            r_pro = sum(mat[ir])
        elif ir == len_matrix:
            r_pro = 0
            for i_list in range(ir):
                r_pro = r_pro + mat[i_list][-1]
        else:
            r_pro = sum(mat[ir][ir:])
            for i_list in range(ir):
                r_pro = r_pro + mat[i_list][ir - 1]
        r_mat.append(r_pro)
    return r_mat


def get_m_mat(mat, r_mat, len_matrix):
    m_mat = np.ones([len_matrix, len_matrix]) + 0.0001
    for irow in range(len_matrix):
        for icol in range(irow, len_matrix):
            m_mat[irow, icol] = mat[irow, icol] - (r_mat[irow] + r_mat[icol + 1]) / (len_matrix - 1)

    return m_mat


def single_single_upgma(row, pro_pre, pro_post, labels, name_all, min_dis, distances,
                        name_distances, pro_pre_dis, pro_post_dis):
    """This is a function using UPGMA method to update the itol file input if the minimum distance is
            between two single proteins"""
    # combine protein IDs
    name_all[row] = '(' + pro_pre + ',' + pro_post + ')'
    # append the new combination to the first layer of labels
    labels[0].append([pro_pre, pro_post])
    # append distance and combination of protein ID
    distances[0].append(min_dis)
    name_distances[row] = '(' + pro_pre_dis + ':' + str(min_dis / 2) + ',' + pro_post_dis + ':' + \
                          str(min_dis / 2) + ')'
    return name_all, labels, distances, name_distances


def single_single_nj(row, pro_pre, pro_post, labels, name_all, min_dis1, min_dis2, distances1,
                     distances2, name_distances, pro_pre_dis, pro_post_dis):
    """This is a function using NJ method to update the itol file input if the minimum distance is
        between two single proteins"""
    # combine protein IDs
    name_all[row] = '(' + pro_pre + ',' + pro_post + ')'
    # append the new combination to the first layer of labels
    labels[0].append([pro_pre, pro_post])
    # append ne distances
    distances1[0].append(min_dis1)
    distances2[0].append(min_dis2)
    name_distances[row] = '(' + pro_pre_dis + ':' + str(min_dis1) + ',' + pro_post_dis + ':' \
                          + str(min_dis2) + ')'
    return name_all, labels, distances1, distances2, name_distances


def get_pro_name_no_num(pro):
    for ind_pro_get, pro_get_str in enumerate(pro):
        if pro_get_str == 'p':
            pro_get_ind = ind_pro_get
            break
    return pro_get_ind


def get_num_bra(pro_post):
    count_max = 0
    count_bra = 0
    # count the number of brackets
    for i_ch in pro_post:
        if i_ch == '(':
            count_bra += 1
        elif i_ch == ')':
            count_max = np.max([count_max, count_bra])
            count_bra -= 1
    return count_max


# 获取复合蛋白质的第一个蛋白质名
def get_1st_pro(pro):
    pro_nobra = pro.replace('(', '')
    pro_nobra = pro_nobra.replace(')', '')
    # 获取逗号前的蛋白质名
    pro_1st = ''
    for letter_str in pro_nobra:
        if letter_str != ',':
            pro_1st = pro_1st + letter_str
        else:
            break
    return pro_1st


# 获得上一层的位置
def get_num_prelayer(label_prelayer, pro_1st):
    for i_sublabel, sublabel in enumerate(label_prelayer):
        pro = sublabel[0]
        # 如果是p
        if pro[0] != 'p':
            pro_get_ind = get_pro_name_no_num(pro)
            pro = pro[pro_get_ind:]
        if pro_1st == pro:
            pro_comb = str(i_sublabel) + pro_1st
    return pro_comb, i_sublabel


def single_comp_upgma(row, pro_pre, pro_post, labels, name_all, min_dis, distances,
                      name_distances, pro_pre_dis, pro_post_dis):
    """This is a function using UPGMA method to update the itol file input if the second ID is a
    combination of multiple protein IDs"""
    # identify the index of layer where the second ID belongs by counting the number of brackets.
    count_bra = get_num_bra(pro_post)
    # identify the position of the ID combination in the previous layer
    pro_1st = get_1st_pro(pro_post)
    pro_post_comb, i_sublabel = get_num_prelayer(labels[count_bra - 1], pro_1st)
    # if the new layer does not exist in the label list
    if len(labels) <= count_bra:
        labels.append([])
    # append the new id
    labels[count_bra].append([pro_pre, pro_post_comb])
    name_all[row] = '(' + pro_pre + ',' + pro_post + ')'
    # if the new layer does not exist in the distance list
    if len(distances) <= count_bra:
        distances.append([])
    # append the new distance
    distances[count_bra].append(min_dis)
    # compute the length if the node by minus the distance by the distance of the former layer
    dis_comp = (min_dis - distances[count_bra - 1][i_sublabel])
    name_distances[row] = '(' + pro_pre_dis + ':' + str(min_dis) + ',' + pro_post_dis + ':' + \
                          str(dis_comp) + ')'
    return name_all, labels, distances, name_distances


def single_comp_nj(row, pro_pre, pro_post, labels, name_all, min_dis1, min_dis2, distances1,
                   distances2, name_distances, pro_pre_dis, pro_post_dis):
    # 确定在哪一层
    count_bra = get_num_bra(pro_post)
    # 去掉括号的复合蛋白质名，对比确定在上一层的第几位
    pro_1st = get_1st_pro(pro_post)
    # 对比蛋白质名
    pro_post_comb, i_sublabel = get_num_prelayer(labels[count_bra - 1], pro_1st)
    if len(labels) <= count_bra:
        labels.append([])
    labels[count_bra].append([pro_pre, pro_post_comb])
    name_all[row] = '(' + pro_pre + ',' + pro_post + ')'
    if len(distances1) <= count_bra:
        distances1.append([])
        distances2.append([])
    distances1[count_bra].append(min_dis1)
    distances2[count_bra].append(min_dis2)
    name_distances[row] = '(' + pro_pre_dis + ':' + str(min_dis1) + ',' + pro_post_dis + ':' + str(min_dis2) + ')'

    return name_all, labels, distances1, distances2, name_distances


def comp_single_upgma(row, pro_pre, pro_post, labels, name_all, min_dis, distances, name_distances,
                      pro_pre_dis, pro_post_dis):
    """This is a function using UPGMA method to update the itol file input if the first ID is a
    combination of multiple protein IDs"""
    # identify the index of layer where the second ID belongs by counting the number of brackets.
    count_bra = get_num_bra(pro_pre)
    # identify the position of the ID combination in the previous layer
    pro_1st = get_1st_pro(pro_pre)
    # compare the protein ID with the previous IDs
    pro_pre_comb, i_sublabel = get_num_prelayer(labels[count_bra - 1], pro_1st)
    # identify the position of the ID combination in the previous layer
    if len(labels) <= count_bra:
        labels.append([])
    # append the new id
    labels[count_bra].append([pro_pre_comb, pro_post])
    name_all[row] = '(' + pro_pre + ',' + pro_post + ')'
    # if the new layer does not exist in the distance list
    if len(distances) <= count_bra:
        distances.append([])
    # append the new distance
    distances[count_bra].append(min_dis)
    # compute the length if the node by minus the distance by the distance of the former layer
    dis_comp = (min_dis - distances[count_bra - 1][i_sublabel])
    name_distances[row] = '(' + pro_pre_dis + ':' + str(dis_comp) + ',' + pro_post_dis + ':' + \
                          str(min_dis) + ')'
    return name_all, labels, distances, name_distances


def comp_single_nj(row, pro_pre, pro_post, labels, name_all, min_dis1, min_dis2, distances1,
                   distances2, name_distances, pro_pre_dis, pro_post_dis):
    # 确定在哪一层
    count_bra = get_num_bra(pro_pre)
    # 去掉括号的复合蛋白质名，对比确定在上一层的第几位
    pro_1st = get_1st_pro(pro_pre)
    # 对比蛋白质名
    pro_pre_comb, i_sublabel = get_num_prelayer(labels[count_bra - 1], pro_1st)
    if len(labels) <= count_bra:
        labels.append([])
    labels[count_bra].append([pro_pre_comb, pro_post])
    name_all[row] = '(' + pro_pre + ',' + pro_post + ')'
    if len(distances1) <= count_bra:
        distances1.append([])
        distances2.append([])
    distances1[count_bra].append(min_dis1)
    distances2[count_bra].append(min_dis2)
    name_distances[row] = '(' + pro_pre_dis + ':' + str(min_dis1) + ',' + pro_post_dis + ':' + str(min_dis2) + ')'
    return name_all, labels, distances1, distances2, name_distances


def comp_comp_upgma(row, pro_pre, pro_post, labels, name_all, min_dis, distances, name_distances,
                    pro_pre_dis, pro_post_dis):
    # identify the index of layer where the second ID belongs by counting the number of brackets.
    count_bra_pre = get_num_bra(pro_pre)
    count_bra_post = get_num_bra(pro_post)
    # identify the position of the ID combination in the previous layer
    pro_1st_pre = get_1st_pro(pro_pre)
    pro_1st_post = get_1st_pro(pro_post)
    # compare the protein ID with the previous IDs
    pro_pre_comb, i_sublabel_pre = get_num_prelayer(labels[count_bra_pre - 1], pro_1st_pre)
    pro_post_comb, i_sublabel_post = get_num_prelayer(labels[count_bra_post - 1], pro_1st_post)
    count_bra_big = np.max([count_bra_pre, count_bra_post])
    # compare the protein ID with the previous IDs
    if len(labels) <= count_bra_big:
        labels.append([])
    # append the new id
    labels[count_bra_big].append([pro_pre_comb, pro_post_comb])
    name_all[row] = '(' + pro_pre + ',' + pro_post + ')'
    # if the new layer does not exist in the distance list
    if len(distances) <= count_bra_big:
        distances.append([])
    # append the new distance
    distances[count_bra_big].append(min_dis)
    # compute the length if the node by minus the distance by the distance of the former layer
    dis_comp_pre = (min_dis - distances[count_bra_pre - 1][i_sublabel_pre])
    dis_comp_post = (min_dis - distances[count_bra_post - 1][i_sublabel_post])
    name_distances[row] = '(' + pro_pre_dis + ':' + str(dis_comp_pre) + ',' + pro_post_dis + ':' + \
                          str(dis_comp_post) + ')'
    return name_all, labels, distances, name_distances


def comp_comp_nj(row, pro_pre, pro_post, labels, name_all, min_dis1, min_dis2, distances1, distances2,
                 name_distances, pro_pre_dis, pro_post_dis):
    # 确定在哪一层
    count_bra_pre = get_num_bra(pro_pre)
    count_bra_post = get_num_bra(pro_post)

    # 去掉括号的复合蛋白质名，对比确定在上一层的第几位
    pro_1st_pre = get_1st_pro(pro_pre)
    pro_1st_post = get_1st_pro(pro_post)
    # 对比蛋白质名
    pro_pre_comb, i_sublabel_pre = get_num_prelayer(labels[count_bra_pre - 1], pro_1st_pre)
    pro_post_comb, i_sublabel_post = get_num_prelayer(labels[count_bra_post - 1], pro_1st_post)
    count_bra_big = np.max([count_bra_pre, count_bra_post])
    if len(labels) <= count_bra_big:
        labels.append([])
    labels[count_bra_big].append([pro_pre_comb, pro_post_comb])
    name_all[row] = '(' + pro_pre + ',' + pro_post + ')'
    if len(distances1) <= count_bra_big:
        distances1.append([])
        distances2.append([])
    distances1[count_bra_big].append(min_dis1)
    distances2[count_bra_big].append(min_dis2)

    name_distances[row] = '(' + pro_pre_dis + ':' + str(min_dis1) + ',' + pro_post_dis + ':' + str(min_dis2) + ')'
    return name_all, labels, distances1, distances2, name_distances


# In[ ]:


def align_y2x(labels, distances):
    distances_labels = copy.deepcopy(labels)
    for i_dis in range(len(distances_labels[0])):
        distances_labels[0][i_dis][0] = distances[0][i_dis]
        distances_labels[0][i_dis][1] = distances[0][i_dis]


def not_in_last_layer(k_labels, labels, ind_i, ind_j, ind_k):
    labels_pop = copy.deepcopy(labels)
    for ind_pop in range(ind_i - 1, len(labels)):
        labels_pop.pop(-1)
    for ind_i_pop, i_labels_pop in enumerate(labels_pop):
        for ind_j_pop, j_labels_pop in enumerate(i_labels_pop):
            for ind_k_pop, k_labels_pop in enumerate(j_labels_pop):
                if k_labels_pop[0] != 'p':
                    k_labels_ind_get = get_pro_name_no_num(k_labels_pop)
                    k_labels_pop = k_labels_pop[k_labels_ind_get:]
                if k_labels == k_labels_pop:
                    i_out = ind_i_pop
                    j_out = ind_j_pop
    return i_out, j_out


def align_name_num(name_all, labels, distances):
    labels_num = copy.deepcopy(labels)
    labels_y1 = copy.deepcopy(labels)
    labels_y2 = copy.deepcopy(labels)
    pro_all = name_all[0]
    pro_all_nobra = pro_all.replace('(', '')
    pro_all_nobra = pro_all_nobra.replace(')', '')
    pro_all_nobra = pro_all_nobra + ','
    count_p = 0
    pro_single = ''
    for pro_str in pro_all_nobra:
        if pro_str == ',':
            count_p += 1
            for ind_i, i_labels in enumerate(labels_num):
                for ind_j, j_labels in enumerate(i_labels):
                    for ind_k, k_labels in enumerate(j_labels):
                        if pro_single == k_labels:
                            labels_num[ind_i][ind_j][ind_k] = count_p - 1
                            labels_y1[ind_i][ind_j][ind_k] = 0
                            labels_y2[ind_i][ind_j][ind_k] = distances[ind_i][ind_j]
            pro_single = ''
        else:
            pro_single = pro_single + pro_str

    for ind_i, i_labels in enumerate(labels_num):
        for ind_j, j_labels in enumerate(i_labels):
            for ind_k, k_labels in enumerate(j_labels):
                if isinstance(k_labels, str):
                    label_num_get_ind = get_pro_name_no_num(k_labels)
                    ind_j_num = int(k_labels[:label_num_get_ind])
                    if ind_j_num >= len(labels[ind_i - 1]):
                        ind_i_ref, ind_j_ref = not_in_last_layer(k_labels[label_num_get_ind:], labels, ind_i, ind_j,
                                                                 ind_k)
                    else:
                        labels_ind_comp = labels[ind_i - 1][ind_j_num][0]
                        if labels_ind_comp[0] != 'p':
                            label_comp_get_ind = get_pro_name_no_num(labels_ind_comp)
                            labels_ind_comp = labels_ind_comp = labels[ind_i - 1][ind_j_num][0][label_comp_get_ind:]
                        if k_labels[label_num_get_ind:] == labels_ind_comp:
                            ind_i_ref = ind_i - 1
                            ind_j_ref = ind_j_num
                        else:
                            ind_i_ref, ind_j_ref = not_in_last_layer(k_labels[label_num_get_ind:], labels, ind_i, ind_j,
                                                                     ind_k)

                    num_pre = labels_num[ind_i_ref][ind_j_ref][0]
                    num_post = labels_num[ind_i_ref][ind_j_ref][1]
                    labels_num[ind_i][ind_j][ind_k] = (num_pre + num_post) / 2
                    labels_y1[ind_i][ind_j][ind_k] = labels_y2[ind_i_ref][ind_j_ref][0]
                    labels_y2[ind_i][ind_j][0] = distances[ind_i][ind_j]
                    labels_y2[ind_i][ind_j][1] = distances[ind_i][ind_j]
    return labels_num, labels_y1, labels_y2


def get_similarity_clust(cluster_num, distances):
    dis_squeez = []
    for i_dis in distances:
        for j_dis in i_dis:
            dis_squeez.append(j_dis)
    dis_squeez.sort()
    dis_sublim = dis_squeez[- cluster_num]
    dis_uplim = dis_squeez[- cluster_num + 1]
    similarity = (dis_sublim + dis_uplim) / 2
    return similarity


def get_xcut(labels_num, labels_y1, labels_y2, cluster_num, distances):
    similarity = get_similarity_clust(cluster_num, distances)
    count_clust = 0
    labels_num_squeez = []
    labels_y1_squeez = []
    labels_y2_squeez = []
    x_cut = []
    for i_labels_num in range(len(labels_num)):
        for j_labels_num in range(len(labels_num[i_labels_num])):
            labels_num_squeez.append(labels_num[i_labels_num][j_labels_num])
            labels_y1_squeez.append(labels_y1[i_labels_num][j_labels_num])
            labels_y2_squeez.append(labels_y2[i_labels_num][j_labels_num])
    for i_all_num in range(len(labels_num_squeez)):
        for j_all_num in range(len(labels_num_squeez[i_all_num])):
            if labels_y1_squeez[i_all_num][j_all_num] < similarity < \
                    labels_y2_squeez[i_all_num][j_all_num]:
                count_clust += 1
                if labels_y1_squeez[i_all_num][j_all_num] == 0:
                    x_mid = labels_num_squeez[i_all_num][j_all_num]
                else:
                    mid_ind = labels_y2_squeez.index([labels_y1_squeez[i_all_num][j_all_num],
                                                      labels_y1_squeez[i_all_num][j_all_num]])
                    y1_pair = labels_y1_squeez[mid_ind]
                    x_mid_pair = labels_num_squeez[mid_ind]
                    x_mid = x_mid_pair[1]
                    count_while = 0
                    while y1_pair[1] != 0:
                        count_while += 1
                        mid_ind_new = labels_y2_squeez.index([y1_pair[1], y1_pair[1]])
                        y1_pair_new = labels_y1_squeez[mid_ind_new]
                        x_mid_pair_new = labels_num_squeez[mid_ind_new]
                        x_mid = x_mid_pair_new[1]
                        y1_pair = y1_pair_new
                        x_mid_pair = x_mid_pair_new
                        count_while = 0
                x_cut.append(x_mid)
                x_cut.sort()
    return x_cut, similarity


def get_random_colors(cluster_num):
    color_list = []
    for i_color in range(cluster_num):
        color_list.append([random.random(), random.random(), random.random()])
    return color_list


def get_plot_color(x_plot3, y_plot1, y_plot2, y_plot3, similarity, x_cut, color_list):
    if y_plot3[0] < similarity:
        for ind_x_cut, value_x_cut in enumerate(x_cut):
            if x_plot3[1] <= value_x_cut:
                color_plot3 = color_list[ind_x_cut]
                color_plot1 = color_list[ind_x_cut]
                color_plot2 = color_list[ind_x_cut]
                break
    else:
        color_plot3 = 'k'
        if y_plot2[0] > similarity:
            color_plot2 = 'k'
        elif y_plot2[0] < similarity and y_plot2[1] > similarity:
            color_plot2 = 'k'
        else:
            for ind_x_cut, value_x_cut in enumerate(x_cut):
                if x_plot3[1] <= value_x_cut:
                    color_plot2 = color_list[ind_x_cut]
                    break
        if y_plot1[0] > similarity:
            color_plot1 = 'k'
        elif y_plot1[0] < similarity and y_plot1[1] > similarity:
            color_plot1 = 'k'
        else:
            for ind_x_cut, value_x_cut in enumerate(x_cut):
                if x_plot3[0] <= value_x_cut:
                    color_plot1 = color_list[ind_x_cut]
                    break
    return color_plot1, color_plot2, color_plot3


def get_xticklabel(name_all):
    name_nobra = name_all.replace('(', '')
    name_nobra = name_nobra.replace(')', '')
    xticklabels_list = []
    name_single = ''
    for name_single_str in name_nobra:
        if name_single_str == ',':
            xticklabels_list.append(name_single)
            name_single = ''
        else:
            name_single = name_single + name_single_str
    xticklabels_list.append(name_single)
    return xticklabels_list


def add_p(pro_list):
    ppro_list = []
    for single_pro in pro_list:
        ppro = 'p' + single_pro
        ppro_list.append(ppro)
    return ppro_list


def plot_tree(cluster_num, name_all, labels, distances, figsize=None):
    if figsize is None:
        figsize = [85, 10]
    labels_num, labels_y1, labels_y2 = align_name_num(name_all, labels, distances)
    xticklabels = get_xticklabel(name_all[0])
    x_cut, similarity = get_xcut(labels_num, labels_y1, labels_y2, cluster_num, distances)
    color_list = get_random_colors(cluster_num)
    fig, ax = plt.subplots(figsize=figsize)
    for i_plot in range(len(labels_y1)):
        for j_plot in range(len(labels_y1[i_plot])):
            x_plot1 = [labels_num[i_plot][j_plot][0], labels_num[i_plot][j_plot][0]]
            x_plot2 = [labels_num[i_plot][j_plot][1], labels_num[i_plot][j_plot][1]]
            x_plot3 = labels_num[i_plot][j_plot]
            y_plot1 = [labels_y1[i_plot][j_plot][0], labels_y2[i_plot][j_plot][0]]
            y_plot2 = [labels_y1[i_plot][j_plot][1], labels_y2[i_plot][j_plot][1]]
            y_plot3 = [labels_y2[i_plot][j_plot][0], labels_y2[i_plot][j_plot][1]]
            color_plot1, color_plot2, color_plot3 = get_plot_color(x_plot3, y_plot1, y_plot2, y_plot3,
                                                                   similarity, x_cut, color_list)
            ax.plot(x_plot1, y_plot1, color=color_plot1, linewidth=1)
            ax.plot(x_plot2, y_plot2, color=color_plot2, linewidth=1)
            ax.plot(x_plot3, y_plot3, color=color_plot3, linewidth=1)

    ax.set_xticks(range(len(xticklabels)))
    ax.set_xticklabels(xticklabels, rotation=60, fontsize=12)
    return fig, ax
