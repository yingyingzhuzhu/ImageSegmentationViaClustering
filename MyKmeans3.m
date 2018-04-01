function [ClusterIm_kMeans,CCIm_kMeans] = MyKmeans3(Im,ImType,Type,NumClusts,k)

% tic
warning off

ImType = inputParser;
paramName = 'ImType';
defaultVal = 'RGB';
addParameter(ImType,paramName,defaultVal)

NumClusts = inputParser;
paramName = 'NumClusts';
defaultVal = 2;
addParameter(NumClusts,paramName,defaultVal)

[m,n,p] = size(Im);

Im_Reshape = reshape(Im,m*n,p);

X = double(Im_Reshape);

ClusterIm_kMeans = reshape(kmeans(X, k, 'Distance','city'),m,n);

switch Type
    case ('RGB')
        CCIm_kMeans = getCCIm(ClusterIm_kMeans);%CCIm_GMM_temp;
    case ('Hyper')
        CCIm_kMeans = [];
end
