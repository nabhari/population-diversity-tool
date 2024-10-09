#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#author: Niloufar Abhari 
#Licensed under the Apache License 2.0

import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
import time
import os
import matplotlib.cm as cm


################################################# READING DATA #######################################################

def load_frequencies(d_path='test_example.csv' ,pops_col_name = 'pops'):
    data = pd.read_csv(d_path)
    data2 = data.set_index(pops_col_name)
    #print(data2)
    return data2

########################################### CORRELATION PLOTS ##################################################

    
def scatter_plots_all_pairs(k=3, path_df='test_example__k2_output.csv', col_index='Unnamed: 0'):
    data = pd.read_csv(os.path.expanduser(path_df))
    data2 = data.set_index(col_index)

    # create dataframe without 'Set_of_pops' column
    s = pd.DataFrame(data2.drop(columns=['Set_of_pops']), columns=data2.columns[1:])

    # Get all pairs of column names
    column_pairs = list(itertools.combinations(s.columns, 2))

    # Generate a big figure with subplots
    num_plots = len(column_pairs)
    num_cols = 4  # Adjust the number of columns in the big figure as per your preference
    num_rows = num_plots // num_cols + 1

    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(20, 20))
    fig.suptitle(f"Scatter plots showing the correlation of all pairs of measures for population sets of size {k}", fontsize=16)
    plt.subplots_adjust(hspace=0.85,wspace=0.45)
    colors = cm.Greens(np.linspace(0.2, 0.9, len(s)))
    for idx, (col1, col2) in enumerate(column_pairs):
        ax = axes[idx // num_cols, idx % num_cols]
        ax.scatter(s[col1], s[col2], alpha=0.4,c=colors)
        ax.set_xlabel(col1)
        ax.set_ylabel(col2)
        # Calculate Pearson's correlation coefficient
        #pearsons_coefficient = np.corrcoef(s[col1], s[col2])[1, 0]
        pearsons_coefficient =round(np.corrcoef(s[col1], s[col2])[1, 0], 15)
        #print((col1,col2),pearsons_coefficient)
        #slope, residuals = np.polyfit(s[col1], s[col2], 1) # uncomment to calculate the slope and residu of the regression line
        #print(slope, residuals)
        text_pos_x = np.min(s[col1]) + 0.1 * (np.max(s[col1]) - np.min(s[col1]))
        text_pos_y = np.max(s[col2]) - 0.1 * (np.max(s[col2]) - np.min(s[col2]))
        ax.text(text_pos_x, text_pos_y, f"Pearson's correlation:\n{pearsons_coefficient:.2f}", fontsize=8,
                bbox=dict(facecolor='white', alpha=0.5))

    # Remove empty subplots
    for idx in range(num_plots, num_cols * num_rows):
        ax = axes[idx // num_cols, idx % num_cols]
        ax.axis('off')

    #plt.show() #uncomment to show the plot automatically
    #  Save the big figure with subplots
    plt.savefig("scatter.png", dpi=300, format="png")
    plt.close()


##############################################DIVERSITY_MEASURES############################################

############################################## Pairwise Differencing #########################################
def pi_minus_pj_squared_all_efficient(df):
    # Calculate pairwise distances for each column
    distances = {}

    # Iterate over the columns
    for col in df.columns:
        column_values = df[col].values
        # Calculate pairwise distances within the column using combinations
        pairwise_distances = np.subtract.outer(column_values, column_values) ** 2
        column_distances = 2*pairwise_distances[np.triu_indices(len(column_values), k=1)]
        
        # Store the pairwise distances for the current column
        distances[col] = column_distances

    # Create a DataFrame from the distances dictionary
    distances_df = pd.DataFrame(distances)
    #print(distances)
    return distances_df
def Het_Interpop_square(freqs):
    #freqs : an array with #subpopulations rows and #loci column, each entry is the frequency of 1 in pop_i for locus_j
    T = pi_minus_pj_squared_all_efficient(freqs)
    T = np.array(T)
    # Compute the score for the column by summing the pairwise distances and dividing by 2*(number of elements in the freqs)^2
    score = (np.sum(T)) / (len(freqs) ** 2)
    return score
def SSD_Interpop_square(freqs):
    T = pi_minus_pj_squared_all_efficient(freqs)
    # Take the maximum value of each column
    max_values = T.max()
    # Sum over the maximum values
    return max_values.sum()
############################################## POOLING ########################################################
def Het_Pooling(freqs):
    #freqs : an array with #subpopulations rows and #loci column, each entry is the frequency of 1 in pop_i for locus_j
    # Initialize the score variable
    score = 0
    # Iterate over the columns
    for col in freqs.columns:
        # Calculate the average for the current column
        average = freqs[col].mean()
        score += average * (1 - average)

    return 2*score
def SSD_Pooling(freqs):
    # Initialize the score variable
    score = 0
    # Iterate over the columns
    for col in freqs.columns:
        # Calculate the average for the current column
        average = freqs[col].mean()
        
        # Check if the average is above the threshold
        if average > 0 and average < 1:
            # Increment the score by 1
            score += 1
    return score

############################################## AVERAGING #############################################
def Het_Avg(freqs):
    # Perform element-wise multiplication by (1-p) and calculate column sums
    column_sums = freqs.apply(lambda col: (2 * col * (1 - col)).sum())

    # Calculate the final score
    score = column_sums.sum()

    return score/len(np.array(freqs))
def SSD_Avg(freqs):
    # Count the number of elements in each column that satisfy the condition
    score = freqs.apply(lambda col: ((col > 0) & (col < 1)).sum())

    # Sum over the counts per column
    return score.sum()/len(np.array(freqs))

############################################## FIXING #################################################
def Het_Fixing(freqs):
    # Calculate the score for each column
    column_scores = {}

    for col in freqs.columns:
        column_values = freqs[col].values
        #The outer product of a 1D array with itself results in a 2D matrix where each element at position [i, j] 
        # is the product of the i-th element of the first array and the j-th element of the second array.
        outer_product1 = np.outer(column_values, column_values)
        outer_product2 = np.outer(1-column_values, 1-column_values)
        res_i = 1 - outer_product1 - outer_product2
        score = 2 * np.sum(np.triu(res_i, k=1))
        column_scores[col] = score
    # Calculate the total score
    total_score = sum(column_scores.values()) / (len(freqs)**2)
    return total_score
def SSD_Fixing2(freqs):
    #the other version is more efficient
     # Calculate the score for each column
    column_scores = {}

    for col in freqs.columns:
        column_values = freqs[col].values

        prod_ones = np.prod(column_values)
        prod_zeros = np.prod(1-column_values)
        column_scores[col] = 1 - prod_ones - prod_zeros
    # Calculate the total score
    total_score = sum(column_scores.values())
    return total_score
def SSD_Fixing(freqs):
    '''
        freqs is a dataframe with one pop per row and the proportion of each feature per pop in columns
    '''
    #all feature weights are set to be one
    mu_f=np.ones(np.array(freqs).shape[1])
    multi_ones = np.prod(np.array(freqs),axis=0)
    multi_zeros = np.prod(np.array(1-freqs),axis=0)
    score_per_feature = 1 - (multi_zeros+multi_ones)
    weighted_sub_maxmin = np.multiply(mu_f,score_per_feature)
    score_fixing = weighted_sub_maxmin.sum()
    return score_fixing
##############################################BRUTE_FORCE############################################
def findsubsets(s, k):
    return list(itertools.combinations(s, k))

def brute_force_all_HETandSSD(pop_freqs, k=2, save_file_as='brute_force_output.csv'):
    #assuming that the pop names are se as the index column in pop_freqs,
    #we get all possible k-size subset of size k of pops 
    list_of_sets = findsubsets(list(pop_freqs.index),k)
    scores=[]
    i = 0
    for item in list_of_sets:
        i+=1
        print("This is the",i,"-th iteration")
        temp_max = pop_freqs.loc[list(item)]
        scores.append([item,Het_Pooling(temp_max),Het_Avg(temp_max),Het_Interpop_square(temp_max),Het_Fixing(temp_max),
                       SSD_Pooling(temp_max),SSD_Avg(temp_max),SSD_Interpop_square(temp_max),SSD_Fixing(temp_max)])


    df = pd.DataFrame(scores,columns=['Set_of_pops','Het Pooling','Het Averaging','Het Pairwise Differencing','Het Fixing',
                       'SSD Pooling','SSD Averaging','SSD Pairwise Differencing','SSD Fixing'])
    df.to_csv(save_file_as)
    return df

def main_():
    d = load_frequencies()
    # d = load_frequencies('Atlantic-salmon-ONEfreqs_perPop.csv') #uncomment if you want read the Atlantic salmon data    
    brute_force_all_HETandSSD(pop_freqs=d, k=2, save_file_as = 'test_example__k2_output.csv')
    scatter_plots_all_pairs(k=2,path_df='test_example__k2_output.csv', col_index='Unnamed: 0')

################################################# RUN ##################

#Measure the runtime of the code
start_time = time.time()

main_()
###########################################RunTime######################
# Print the runtime of the code
end_time = time.time()
print(f"Runtime: {end_time - start_time} seconds")
 
