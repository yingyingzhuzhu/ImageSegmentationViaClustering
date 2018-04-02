pavia_img = load('PaviaHyperIm.mat');
pavia_img_array = pavia_img.PaviaHyperIm;
[height,width,p] = size(pavia_img_array);
nCluster = height * width / 4;
OCE = 1.0;
for c = 2 : nCluster
    [ClusterIm, CCIm] = MyKmeans3(pavia_img_array, 'ImType', 'Hyper', 'NumClusts', c);
    npClusterIm = py.numpy.array(ClusterIm(:).');
    npClusterIm = npClusterIm.astype('int');
    ClusterIm = npClusterIm.reshape(size(ClusterIm,1), size(ClusterIm,2));
    temp_score = py.MyClustEvalHyper3.MyClustEvalHyper3(ClusterIm);
    disp(c)
    disp(temp_score)
    if temp_score < OCE
        OCE = temp_score;
        cluster = c;
    end
end
fprintf("OCE = ");
disp(OCE);
fprintf("cluster = ");
disp(cluster);