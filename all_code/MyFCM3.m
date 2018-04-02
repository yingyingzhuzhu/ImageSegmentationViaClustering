function ClusterIm = MyFCM3(Im, c)

[m,n,p] = size(Im);

Im_Reshape = reshape(Im,m*n,p);

X = double(Im_Reshape);

[~, U] = fcm(X, c);

maxU = max(U);

indexs = {};
for i = 1:c
    indexs{i} = {find(U(i, :) == maxU)};
end 

fcmImage(1:length(Im))=0; 
for i = 1:c
    fcmImage(indexs{i}{1})= i;
end

ClusterIm = reshape(fcmImage,m,n);

% figure;imagesc(ClusterIm),axis image,colorbar;


end
