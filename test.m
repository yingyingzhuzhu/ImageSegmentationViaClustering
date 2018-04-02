
clear
clc
close all;
tic

load('CCIm_kMeans_Seg.mat')
load('ClusterIm_kMeans_Seg.mat')
load('ImagesForTest.mat')
load('numCluster.mat')
load('Seg.mat')

for i = 1: 198
    Im = Image{i};
    Seg1 = Seg{1}{i};
    Seg2 = Seg{2}{i};
    Seg3 = Seg{3}{i};
    CCIm1 = CCIm_kMeans_Seg{1}{i};
    CCIm2 = CCIm_kMeans_Seg{2}{i};
    CCIm3 = CCIm_kMeans_Seg{3}{i};
    numCluster1 = list_numCluster(i,1);
    numCluster2 = list_numCluster(i,2);
    numCluster3 = list_numCluster(i,3);
    ClusterIm1 = ClusterIm_kMeans_Seg{1}{i};
    ClusterIm2 = ClusterIm_kMeans_Seg{2}{i};
    ClusterIm3 = ClusterIm_kMeans_Seg{3}{i};
    save(sprintf('kMeans%d.mat',i),'Im','Seg1','Seg2','Seg3','CCIm1','CCIm2','CCIm3','numCluster1','numCluster2','numCluster3','ClusterIm1','ClusterIm2','ClusterIm3');
    clear Im;
    clear Seg1; clear Seg2; clear Seg3;
    clear CCIm1; clear CCIm2; clear CCIm3;
    clear numCluster1; clear numCluster2; clear numCluster3;
    clear ClusterIm1; clear ClusterIm2; clear ClusterIm3;
    i
end

toc