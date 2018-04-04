import matplotlib.pyplot as plt
import scipy.io as io

data = io.loadmat("./data/GMM/GMM183.mat")
origin = data["Im"]
seg = data["Seg1"]
ccim = data["CCIm1"]
cluster = data["ClusterIm1"]
plt.subplot(221)
plt.imshow(origin)
plt.subplot(222)
plt.imshow(seg)
plt.subplot(223)
plt.imshow(ccim)
plt.subplot(224)
plt.imshow(cluster)
plt.show()