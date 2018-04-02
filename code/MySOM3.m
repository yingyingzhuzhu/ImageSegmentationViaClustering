function [ClusterIm, CCIm] = MySOM3(Im, Type, c)

[m,n,p] = size(Im);

Im_Reshape = reshape(Im, m*n, p);
Im_Reshape = Im_Reshape';

X = double(Im_Reshape);

net = selforgmap([c 1]);

net = train(net, X);

y = net(X);

classes = vec2ind(y);

ClusterIm = reshape(classes,m,n);

figure;imagesc(ClusterIm),axis image,colorbar

switch Type
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
    
end

