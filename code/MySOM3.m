function ClusterIm = MySOM3(Im, c)

[m,n,p] = size(Im);

Im_Reshape = reshape(Im, m*n, p);
Im_Reshape = Im_Reshape';

X = double(Im_Reshape);

net = selforgmap([c 1]);

net = train(net, X);

y = net(X);

classes = vec2ind(y);

ClusterIm = reshape(classes,m,n);

figure;imagesc(ClusterIm),axis image,colorbar;

end
