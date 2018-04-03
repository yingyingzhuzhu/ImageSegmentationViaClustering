import scipy.io as io
import os
import numpy
from all_code.MyClustEvalRGB3 import MyClustEvalRGB3

dir_list = os.listdir("./data/GMM")

value = []

# using clusterIm and gt, mean is 0.804 and stdv is 0.125

for file in dir_list:
    tmp_data = io.loadmat("./data/GMM/" + file)
    clusters = ["ClusterIm1", "ClusterIm2", "ClusterIm3"]
    ccims = ["CCIm1", "CCIm2", "CCIm3"]
    segs = ["Seg1", "Seg2", "Seg3"]
    print("start to calculate")
    tmp_value = []
    for i in range(3):
        ccim = tmp_data[ccims[i]]
        seg = tmp_data[segs[i]]
        cluster = tmp_data[clusters[i]]
        val = MyClustEvalRGB3(ccim, seg)
        tmp_value.append(val)
        print("value is %.3f" % val)
    print("get the smallest value")
    value.append(min(tmp_value))

mean = numpy.mean(value)
stdv = numpy.std(value)
print("mean is %.3f and stdv is %.3f" % (mean, stdv))