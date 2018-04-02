function ClusterIm_kMeans = MyKmeans3(Im, k)

% tic
warning off

[m,n,p] = size(Im);

Im_Reshape = reshape(Im,m*n,p);

X = double(Im_Reshape);

ClusterIm_kMeans = reshape(kmeans(X, k, 'Distance','city'),m,n);

end
    