function [ClusterIm_GMM,CCIm_GMM] = MyGMM3(Im,ImType,Type,NumClusts,k)

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

gm = fitgmdist(X,k,'RegularizationValue',0.1);
idx = cluster(gm,X);

ClusterIm_GMM = reshape(idx,m,n);

% figure,imagesc(ClusterIm_GMM)
% axis image

% CCIm_GMM_temp = zeros(m*n,k);
% 
% for i = 1 : k
%     CCIm_GMM_temp(find(ClusterIm_GMM == i),i) = 1;
% end
% 
% CCIm_GMM_temp = reshape(CCIm_GMM_temp,m,n,k);

switch Type
    case ('RGB')
        CCIm_GMM = getCCIm(ClusterIm_GMM);%CCIm_GMM_temp;
    case ('Hyper')
        CCIm_GMM = [];
end

% toc
end

