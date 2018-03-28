
# coding: utf-8

# In[1]:


import scipy.io as spio
from MyMartinIndex3 import MyMartinIndex3


# In[2]:


def MyClustEvalHyper3(ClusterIm):
    print("Calculating ... ")
    res = 0
    #For Pavia Image
    #load ground truth mask image
    mat_ground_truth_mask = spio.loadmat('PaviaGrTruthMask.mat', squeeze_me=True)
    #load ground truth image
    mat_ground_truth = spio.loadmat('PaviaGrTruth.mat', squeeze_me=True)
    #ground truth mask array
    mat_ground_truth_mask_array = mat_ground_truth_mask.get('PaviaGrTruthMask')
    #ground truth array
    mat_ground_truth_array = mat_ground_truth.get('PaviaGrTruth')
    #ground truth mask size
    mat_ground_truth_height = len(mat_ground_truth_array)
    mat_ground_truth_width = len(mat_ground_truth_array[1])
    #ground truth mask size
    mat_ground_truth_mask_height = len(mat_ground_truth_mask_array)
    mat_ground_truth_mask_width = len(mat_ground_truth_mask_array[1])
    #production result array
    ClusterImAfterFilter = []
    #The evaluation should be performed by multiplying your segmentation by the ground truth Mask. 
    #Don’t use connected components; just use cluster labels. 
    #You can set the number of clusters to 9 since there are 9 ground truthed regions.
    for i in range(mat_ground_truth_mask_height):
        ClusterImAfterFilter.append([])
        for j in range(mat_ground_truth_mask_width):
            ClusterImAfterFilter[i].append(ClusterIm[i][j] * mat_ground_truth_mask_array[i][j])
    res = MyMartinIndex3(mat_ground_truth_array, ClusterImAfterFilter)
    return res

