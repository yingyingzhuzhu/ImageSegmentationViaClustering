%RGB_Kmeans_Experiment.m
FileList = dir('./ImsAndSegs/*.mat');
addpath('./ImsAndSegs');
NFiles   = length(FileList);
% OCEs is to store MartinIndex OCE error of each RGB image
OCEs = double.empty;
% segs is to store seg number(1,2 or 3) with which groundtruth each image get best result
segs = int16.empty;
for n = 1:NFiles
    ImAndSegs = load(FileList(n).name);
    %output image file name
    disp(FileList(n).name);
    Im = ImAndSegs.Im;
    seg123 = {ImAndSegs.Seg1,ImAndSegs.Seg2,ImAndSegs.Seg3};
%     seg1 = ImAndSegs.Seg1;
%     seg2 = ImAndSegs.Seg2;
%     seg3 = ImAndSegs.Seg3;
%     [height,width,p] = size(Im);
%     nCluster = height * width / 4;
    OCE = 1.0;
    %for c = 2 : nCluster
    score = 1.0;
    seg = 1;
    for segIndex = 1:3
        %in each loop, evaluate one of seg1,seg2,seg3
        temp_seg = seg123{segIndex};
        NumClust = getNumClusts(temp_seg);
        disp('NumClust:');
        disp(NumClust);
        %CCIM type change: from matlab to python
        nptemp_seg = py.numpy.array(temp_seg(:).');
        nptemp_seg = nptemp_seg.astype('int');
        temp_seg = nptemp_seg.reshape(size(temp_seg,1), size(temp_seg,2));
        %call algorithm, modify function name to test other four algorithms
        [~, CCIm] = MySOM3(Im, 'ImType', 'RGB', 'NumClusts', NumClust);
        %CCIM type change: from matlab to python
        npCCIm = py.numpy.array(CCIm(:).');
        npCCIm = npCCIm.astype('int');
        CCIm = npCCIm.reshape(size(CCIm,1), size(CCIm,2));
        temp_score = py.MyClustEvalRGB3.MyClustEvalRGB3(CCIm, temp_seg);
        if temp_score < score
            score = temp_score;
            seg = segIndex;
        end
        
        %output temp_score, it is MartinIndex OCE error
        disp('MartinIndex OCE error:');
        disp(score);
    end
    %output cluster number
    %disp(c)
    
    %output temp_score, it is MartinIndex OCE error
    disp('Minimum MartinIndex OCE error:');
    disp(score);
    %get minimum OCE and corresponding cluster number c
%     if temp_score < OCE
%         OCE = temp_score;
%         cluster = c;
%     end
    %end
    %sum up MartinIndex OCE error
    OCEs(n) = score;
    segs(n) = seg;
end

fprintf("OCE means = ");
disp(mean(OCEs));
fprintf("OCE std. dev. = ");
disp(std(OCEs));
% fprintf("cluster = ");
% disp(cluster);


% function getNumClusts is used to compute the number of clusters in a
% segmentation
function numClusts = getNumClusts(segment)
    numClusts = 2;
    [height, width, ~] = size(segment);
    for h = 1:height
        for w = 1 : width
            numClusts_temp = segment(h:h,w:w);
            if numClusts_temp > numClusts
                numClusts = numClusts_temp;
            end
        end
    end

end