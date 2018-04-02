import scipy.io as io
import os
import numpy
from all_code.MyClustEvalRGB3 import MyClustEvalRGB3

dir_list = os.listdir("./data/GMM")

value = []

for file in dir_list:
    tmp_data = io.loadmat("./data/GMM/" + file)
    ccims = ["CCIm1", "CCIm2", "CCIm3"]
    segs = ["Seg1", "Seg2", "Seg3"]
    for i in range(3):
        ccim = tmp_data[ccims[i]]
        seg = tmp_data[segs[i]]
        val = MyClustEvalRGB3(ccim, seg)
        value.append(val)
        print("value is %.3f" % val)

mean = numpy.mean(value)
stdv = numpy.std(value)
print("mean is %.23f and stdv is %.3f" % (mean, stdv))