function ClusterIm_GMM = MyGMM3(Im,k)

% tic
warning off

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

% toc
end

