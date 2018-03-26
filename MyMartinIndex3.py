
# coding: utf-8

# In[169]:


def getClusterAndPixelSet(seg):
    height_seg = len(seg)
    width_seg = len(seg[1])
    clusterAndPixelSet = [] #key: cluster num, value: pixel index
    for i in range(height_seg):
        for j in range(width_seg):
            if seg[i][j] == 0:
                continue
            while len(clusterAndPixelSet) < seg[i][j]:
                clusterAndPixelSet.append(set())
            clusterAndPixelSet[seg[i][j] - 1].add(str(i) + " " + str(j))
    return clusterAndPixelSet


# In[170]:


def getW(clusterAndPixelSet):
    Aj = [] #key: cluster num, value: how many pixels
    for i in range(len(clusterAndPixelSet)):
        Aj.append(len(clusterAndPixelSet[i]))
    sum = getPixelSum(Aj)
    Wj = [] #key: cluster num, value: percentage of pixels in this cluster to total pixels in this image
    for i in range(len(Aj)):
        Wj.append(Aj[i] * 1.0 / sum)
    return Wj


# In[175]:


def getPixelSum(Aj):
    sum = 0
    for i in range(len(Aj)):
        sum = sum + Aj[i]
    return sum


# In[176]:


def getWji(clusterAndPixelSet1, clusterAndPixelSet2):
    Wji = [] # 2D
    for j in range(len(clusterAndPixelSet1)):
        Wji.append([])
        set_Ajj = clusterAndPixelSet1[j]
        sum_Bii = 0
        for i in range(len(clusterAndPixelSet2)):
            set_Bii = clusterAndPixelSet2[i]
            set_intersect = set_Ajj.intersection(set_Bii)
            if len(set_intersect) != 0:
                sum_Bii = sum_Bii + len(set_Bii)
                Wji[j].append(len(set_Bii))
            else:
                Wji[j].append(0)
        for i in range(len(clusterAndPixelSet2)):
                Wji[j][i] = Wji[j][i] * 1.0 / sum_Bii
    return Wji


# In[194]:


def getPartialError(clusterAndPixelSet1, clusterAndPixelSet2, Wji, Wj):
    E = 0
    for j in range(len(clusterAndPixelSet1)):
        sum = 0
        for i in range(len(clusterAndPixelSet2)):
            len1 = len(clusterAndPixelSet1[j].intersection(clusterAndPixelSet2[i]))
            len2 = len(clusterAndPixelSet1[j].union(clusterAndPixelSet2[i]))
            sum += len1 * 1.0 / len2 * Wji[j][i]
        E = E + (1 - sum) * Wj[j]
    return E


# In[201]:


def MyMartinIndex3(seg1, seg2):
    clusterAndPixelSet1 = getClusterAndPixelSet(seg1)
    clusterAndPixelSet2 = getClusterAndPixelSet(seg2)
    Wj = getW(clusterAndPixelSet1)
    Wi = getW(clusterAndPixelSet2)
    Wji = getWji(clusterAndPixelSet1, clusterAndPixelSet2)
    Wij = getWji(clusterAndPixelSet2, clusterAndPixelSet1)
    E1 = getPartialError(clusterAndPixelSet1, clusterAndPixelSet2, Wji, Wj)
    E2 = getPartialError(clusterAndPixelSet2, clusterAndPixelSet1, Wij, Wi)
    res = E1
    if E1 > E2:
        res = E2
    return res

