function [ClusterIm, CCIm] = MyClust3(Im, Algorithm, algorithm, ImType, type, NumClusts, clusts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if strcmp(algorithm, 'Kmeans')== 1
    ClusterIm = MyKmeans3(Im, clusts);
end

if strcmp(algorithm, 'SOM')== 1
   ClusterIm = MySOM3(Im, clusts);
end

if strcmp(algorithm, 'Spectral')
	ClusterIm = MySpectral3(Im, clusts);
end

if strcmp(algorithm, 'GMM')
	ClusterIm = MyGMM3(Im, clusts);
end

if strcmp(algorithm, 'FCM')
	ClusterIm = MyFCM3(Im, clusts);
end

switch type
    case ('RGB')
        CCIm = getCCIm(ClusterIm);%CCIm_GMM_temp;
    case ('Hyper')
        CCIm = [];
end

end

function CCIm = getCCIm(ClusterIm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT: ClusterIm
%%% OUTPUT: CCIm

%%% NOTE: minCluIdx >= 0, or it will crash, MAR 27 730AM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minCluIdx = min(min(ClusterIm));
if( minCluIdx == 0)
    ClusterIm = ClusterIm + 1;
    minCluIdx = minCluIdx + 1;
end
CCIm = ClusterIm;
maxCluIdx = max(max(ClusterIm));
numClusters = maxCluIdx - minCluIdx + 1;
CC = cell(numClusters,1);
for n = minCluIdx:maxCluIdx
    tempIm = ClusterIm;
    tempIm(find(tempIm == n)) = -1;
    tempIm(find(tempIm ~= -1)) = -2;
    tempIm(find(tempIm == -1)) = 1;
    tempIm(find(tempIm == -2)) = 0;
    CC{n,1} = bwconncomp(tempIm);
end

ctr = numClusters;
for n = minCluIdx:maxCluIdx
    temp = CC{n,1}.PixelIdxList;
    if(size(temp,2)==1)
        continue;
    else
        endIdx = size(temp,2);
        for m = 2:endIdx
            ctr = ctr + 1;
            CCIm(temp{1,m}) = ctr;
        end
    end
end



