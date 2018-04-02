function [ClusterIm, CCIm] = MyClust3(Im, Algorithm, algorithm, ImType, type, NumClusts, clusts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if strcmp(algorithm, 'Kmeans')== 1
    [ClusterIm, CCIm] = MyKmeans3(Im, type, clusts);
end

if strcmp(algorithm, 'SOM')== 1
   [ClusterIm, CCIm] = MySOM3(Im, type, clusts);
end

if strcmp(algorithm, 'Spectral')
	[ClusterIm, CCIm] = MySpectral3()
end

if strcmp(algorithm, 'GMM')
	[ClusterIm, CCIm] = MyGMM3()
end

if strcmp(algorithm, 'FCM')
	[ClusterIm, CCIm] = MyFCM()
end

end

