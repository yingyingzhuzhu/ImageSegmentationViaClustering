function [ClusterIm, CCIm] = MyKmeans3_2(Im,ImType,NumClusts)
%MYKMEANS Summary of this function goes here
%   Detailed explanation goes here
    if strcmp(ImType,'RGB')== 1
        FileName = Im;
        ImsAndSegs = load(FileName);
        [m,n,d] = size(ImsAndSegs.Im);
        ClusterIm = randi([1,NumClusts],m,n,'uint16');
        CCIm = zeros(m,n,'uint16');
        fin = 0;
        while fin == 0
            fin = 1;
            total = zeros(NumClusts,6,'double');
            mean = zeros(NumClusts,5);
            for i = 1:m
                for j = 1:n
                    total(ClusterIm(i,j),1)=total(ClusterIm(i,j),1)+1;
                    total(ClusterIm(i,j),2)=total(ClusterIm(i,j),2)+i;
                    total(ClusterIm(i,j),3)=total(ClusterIm(i,j),3)+j;
                    total(ClusterIm(i,j),4)=total(ClusterIm(i,j),4)+double(ImsAndSegs.Im(i,j,1));
                    total(ClusterIm(i,j),5)=total(ClusterIm(i,j),5)+double(ImsAndSegs.Im(i,j,2));
                    total(ClusterIm(i,j),6)=total(ClusterIm(i,j),6)+double(ImsAndSegs.Im(i,j,3));
                end
            end
            for i = 1:NumClusts
                num = total(i,1);
                for j = 1:5
                    mean(i,j) = total(i,j+1)/num;
                end
            end
             for i = 1:m
                for j = 1:n
                    nowclust = ClusterIm(i,j);
                    for l = 1:NumClusts
                        nowdist = 0;
                        for k = 1:3
                            nowdist = nowdist + (double(ImsAndSegs.Im(i,j,k))-mean(nowclust,k+2)).^2;
                        end
%                          nowdist = nowdist + (i-mean(nowclust,1)).^2;
%                         nowdist = nowdist + (j-mean(nowclust,2)).^2;
                         tdist = 0;
                        for k = 1:3
                            tdist = tdist + (double(ImsAndSegs.Im(i,j,k))-mean(l,k+2)).^2;
                        end
%                         tdist = tdist + (i-mean(l,1)).^2;
%                         tdist = tdist + (j-mean(l,2)).^2;
                        if tdist < nowdist
                            if nowclust ~= l
                                fin = 0;
                                nowclust = l;
                            end
                        end
                    end
                    ClusterIm(i,j) = nowclust;
                end
             end
        end
    comfin = 0;
    comnum = 1;
    while comfin == 0
        comfin = 1;
        nowclust = 0;
        for i = 1:m
            if comfin == 0
                    break;
            end
            for j = 1:n
                if comfin == 0
                    break;
                end
                if CCIm(i,j) == 0
                    CCIm(i,j) = comnum;
                    nowclust = ClusterIm(i,j);
                    comfin = 0;
                    break;
                end
            end
        end
        scanfin = 0;
        while scanfin == 0
            scanfin = 1;
            for i = 1:m
                for j = 1:n
                    if CCIm(i,j) == comnum
                        if (j+1<=n)&&(CCIm(i,j+1)==0)&&(ClusterIm(i,j+1)==nowclust)
                            CCIm(i,j+1) = comnum;
                            scanfin =0;
                        end
                        if (i+1<=m)&&(CCIm(i+1,j)==0)&&(ClusterIm(i+1,j)==nowclust)
                            CCIm(i+1,j) = comnum;
                            scanfin =0;
                        end
                        if (i-1>=1)&&(CCIm(i-1,j)==0)&&(ClusterIm(i-1,j)==nowclust)
                            CCIm(i-1,j) = comnum;
                            scanfin =0;
                        end
                        if (j-1>=1)&&(CCIm(i,j-1)==0)&&(ClusterIm(i,j-1)==nowclust)
                            CCIm(i,j-1) = comnum;
                            scanfin =0;
                        end
                        if (i-1>=1)&&(j-1>=1)&&(CCIm(i-1,j-1)==0)&&(ClusterIm(i-1,j-1)==nowclust)
                            CCIm(i-1,j-1) = comnum;
                            scanfin =0;
                        end
                        if (i+1<=m)&&(j-1>=1)&&(CCIm(i+1,j-1)==0)&&(ClusterIm(i+1,j-1)==nowclust)
                            CCIm(i+1,j-1) = comnum;
                            scanfin =0;
                        end
                        if (i-1>=1)&&(j+1<=n)&&(CCIm(i-1,j+1)==0)&&(ClusterIm(i-1,j+1)==nowclust)
                            CCIm(i-1,j+1) = comnum;
                            scanfin =0;
                        end
                        if (i+1<=m)&&(j+1<=n)&&(CCIm(i+1,j+1)==0)&&(ClusterIm(i+1,j+1)==nowclust)
                            CCIm(i+1,j+1) = comnum;
                            scanfin =0;
                        end
                    end
                end
            end
        end
        comnum = comnum +1;
    end
    
    figure(1010);
    subplot(2,2,1); imagesc(ImsAndSegs.Im);
    subplot(2,2,2); imagesc(ClusterIm);
    subplot(2,2,3); imagesc(CCIm);
    end
    
    if strcmp(ImType,'Hyper')== 1
        FileName = Im;
        matObj  = matfile(FileName);
        details = whos(matObj);
        [maxBytes, index] = max([details.bytes]);
        maxName = details(index).name;
        tIm =load(FileName,'-mat' ,maxName);
        nIm = getfield(tIm, maxName);
        [m,n,d] = size(nIm);
        ClusterIm = randi([1,NumClusts],m,n,'uint16');
        CCIm = zeros(m,n,'uint16');
        fin = 0;
        while fin == 0
            fin = 1;
            total = zeros(NumClusts,d+3,'double');
            mean = zeros(NumClusts,d+2);
            for i = 1:m
                for j = 1:n
                    total(ClusterIm(i,j),1)=total(ClusterIm(i,j),1)+1;
                    total(ClusterIm(i,j),2)=total(ClusterIm(i,j),2)+i;
                    total(ClusterIm(i,j),3)=total(ClusterIm(i,j),3)+j;
                    for k = 1:d
                        total(ClusterIm(i,j),k+3)=total(ClusterIm(i,j),k+3)+double(nIm(i,j,k));
                    end
                end
            end
            for i = 1:NumClusts
                num = total(i,1);
                for j = 1:d+2
                    mean(i,j) = total(i,j+1)/num;
                end
            end
             for i = 1:m
                for j = 1:n
                    nowclust = ClusterIm(i,j);
                    for l = 1:NumClusts
                        nowdist = 0;
                        for k = 1:d
                            nowdist = nowdist + (double(nIm(i,j,k))-mean(nowclust,k+2)).^2;
                        end
%                          nowdist = nowdist + (i-mean(nowclust,1)).^2;
%                         nowdist = nowdist + (j-mean(nowclust,2)).^2;
                         tdist = 0;
                        for k = 1:d
                            tdist = tdist + (double(nIm(i,j,k))-mean(l,k+2)).^2;
                        end
%                         tdist = tdist + (i-mean(l,1)).^2;
%                         tdist = tdist + (j-mean(l,2)).^2;
                        if tdist < nowdist
                            if nowclust ~= l
                                fin = 0;
                                nowclust = l;
                            end
                        end
                    end
                    ClusterIm(i,j) = nowclust;
                end
             end
        end 
    end
end

