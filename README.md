# Image Segmentation Via Clustering
Machine Learning Project 1

## Overview

Image segmentation refers to a class of algorithms that use unsupervised learning to find “interesting” regions in digital images. In this project, you should use clustering algorithms to cluster two types of images: RGB (color) images and Hyperspectral images. You will use the following clustering algorithms:

• K-means

• Self-Organizing Map

• Fuzzy C-means (FCM)

• Spectral Clustering

• Gaussian Mixture Models (GMM).


## Evaluation Code. 

You will write 2 functions that evaluate clustering algorithms called “EvalClustRGB” and “EvalClustHyper”. 

EvalClustRGB that takes CCIm and a ground-truth segmentation as input and produces a score using the Martin index. 

EvalClustHyper takes ClusterIm as an input.
