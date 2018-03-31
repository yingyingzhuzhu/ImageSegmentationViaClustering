
# coding: utf-8

# In[ ]:


import scipy.io as spio
from MyClustEvalHyper3 import MyClustEvalHyper3
from MyGMM3 import MyGMM3

Im = 'PaviaHyperIm.mat';
pavia_img_array = spio.loadmat(Im, squeeze_me=True).get('PaviaHyperIm')
height = len(pavia_img_array);
width = len(pavia_img_array[0]);
nCluster = height * width / 4;
OCE = 1.0;
for c in range(2, nCluster):
    [ClusterIm, CCIm] = MyGMM3(Im, 'ImType', 'Hyper', 'NumClusts', c); #algorithm
    temp_score = MyClustEvalHyper3(ClusterIm) #evaluate
    OCE = min(OCE, temp_score) #find minimum OCE
print(OCE)

