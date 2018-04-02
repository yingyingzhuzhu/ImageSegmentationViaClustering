function [ClusterIm, CCIm] = MySOM3(Im,NImType,ImType,NNumClusts,NumClusts)
%MYSOM Summary of this function goes here
%   Detailed explanation goes here
    if strcmp(ImType,'RGB')== 1
        [m,n,d] = size(Im);
        ClusterIm = zeros(m,n,'uint16');
        CCIm = zeros(m,n,'uint16');
        wide = ceil(sqrt(double(NumClusts)));
        Grid = randi([0,255],wide,wide,d,'double');
        for ltime = 1: 5000000
            ti = randi(m);
            tj = randi(n);
            nowdist = 0;
            ni = randi(wide);
            nj = randi(wide);
            for k = 1:d
                nowdist = double(nowdist) +(double(Im(ti,tj,k))-double(Grid(ni,nj,k))).^2;
            end
            for i = 1:wide
                for j = 1:wide
                    tdist = 0;
                    for k = 1:d
                        tdist = double(tdist) +(double(Im(ti,tj,k))-double(Grid(i,j,k))).^2;
                    end
                    if tdist < nowdist
                        ni = i;
                        nj = j;
                        nowdist = tdist;
                    end
                end
            end
            for k = 1:d
                Grid(ni,nj,k)= Grid(ni,nj,k)+(double(Im(ti,tj,k))-double(Grid(ni,nj,k)))/3;
                if (nj+1<=wide)
                    Grid(ni,nj+1,k)= Grid(ni,nj+1,k)+(double(Im(ti,tj,k))-double(Grid(ni,nj+1,k)))/9;
                end
                if (ni+1<=wide)
                    Grid(ni+1,nj,k)= Grid(ni+1,nj,k)+(double(Im(ti,tj,k))-double(Grid(ni+1,nj,k)))/9;
                end
                if (nj-1>=1)
                    Grid(ni,nj-1,k)= Grid(ni,nj-1,k)+(double(Im(ti,tj,k))-double(Grid(ni,nj-1,k)))/9;
                end
                if (ni-1>=1)
                    Grid(ni-1,nj,k)= Grid(ni-1,nj,k)+(double(Im(ti,tj,k))-double(Grid(ni-1,nj,k)))/9;
                end
                if (nj+1<=wide) && (ni+1<=wide)
                    Grid(ni+1,nj+1,k)= Grid(ni+1,nj+1,k)+(double(Im(ti,tj,k))-double(Grid(ni+1,nj+1,k)))/27;
                end
                if (nj-1>=1) && (ni+1<=wide)
                    Grid(ni+1,nj-1,k)= Grid(ni+1,nj-1,k)+(double(Im(ti,tj,k))-double(Grid(ni+1,nj-1,k)))/27;
                end
                if (nj+1<=wide) && (ni-1>=1)
                    Grid(ni-1,nj+1,k)= Grid(ni-1,nj+1,k)+(double(Im(ti,tj,k))-double(Grid(ni-1,nj+1,k)))/27;
                end
                if (nj-1>=1) && (ni-1>=1)
                    Grid(ni-1,nj-1,k)= Grid(ni-1,nj-1,k)+(double(Im(ti,tj,k))-double(Grid(ni-1,nj-1,k)))/27;
                end
            end
        end
        
        for ti =1:m
            for tj = 1:n
                 nowdist = 0;
                ni = randi(wide);
                nj = randi(wide);
                ClusterIm(ti,tj)=(ni-1)*wide+nj;
                for k = 1:d
                    nowdist = nowdist +(double(Im(ti,tj,k))-double(Grid(ni,nj,k))).^2;
                end
                for i = 1:wide
                    for j = 1:wide
                        tdist = 0;
                        for k = 1:d
                            tdist = tdist +(double(Im(ti,tj,k))-double(Grid(i,j,k))).^2;
                        end
                        if (tdist < nowdist) && ((i-1)*wide+j<=NumClusts)
                            ClusterIm(ti,tj)=(i-1)*wide+j;
                            nowdist = tdist;
                        end
                    end
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
    subplot(2,2,1); imagesc(Im);
    subplot(2,2,2); imagesc(ClusterIm);
    subplot(2,2,3); imagesc(CCIm);
    end
    
    if strcmp(ImType,'Hyper')== 1
        nIm = Im;
        [m,n,d] = size(nIm);
        ClusterIm = zeros(m,n,'uint16');
        CCIm = zeros(m,n,'uint16');
        wide = ceil(sqrt(double(NumClusts)));
        Grid = randi([0,255],wide,wide,d,'double');
        for ltime = 1: 5000000
            ti = randi(m);
            tj = randi(n);
            nowdist = 0;
            ni = randi(wide);
            nj = randi(wide);
            for k = 1:d
                nowdist = double(nowdist) +(double(nIm(ti,tj,k))-double(Grid(ni,nj,k))).^2;
            end
            for i = 1:wide
                for j = 1:wide
                    tdist = 0;
                    for k = 1:d
                        tdist = double(tdist) +(double(nIm(ti,tj,k))-double(Grid(i,j,k))).^2;
                    end
                    if tdist < nowdist
                        ni = i;
                        nj = j;
                        nowdist = tdist;
                    end
                end
            end
            for k = 1:d
                Grid(ni,nj,k)= Grid(ni,nj,k)+(double(nIm(ti,tj,k))-double(Grid(ni,nj,k)))/3;
                if (nj+1<=wide)
                    Grid(ni,nj+1,k)= Grid(ni,nj+1,k)+(double(nIm(ti,tj,k))-double(Grid(ni,nj+1,k)))/9;
                end
                if (ni+1<=wide)
                    Grid(ni+1,nj,k)= Grid(ni+1,nj,k)+(double(nIm(ti,tj,k))-double(Grid(ni+1,nj,k)))/9;
                end
                if (nj-1>=1)
                    Grid(ni,nj-1,k)= Grid(ni,nj-1,k)+(double(nIm(ti,tj,k))-double(Grid(ni,nj-1,k)))/9;
                end
                if (ni-1>=1)
                    Grid(ni-1,nj,k)= Grid(ni-1,nj,k)+(double(nIm(ti,tj,k))-double(Grid(ni-1,nj,k)))/9;
                end
                if (nj+1<=wide) && (ni+1<=wide)
                    Grid(ni+1,nj+1,k)= Grid(ni+1,nj+1,k)+(double(nIm(ti,tj,k))-double(Grid(ni+1,nj+1,k)))/27;
                end
                if (nj-1>=1) && (ni+1<=wide)
                    Grid(ni+1,nj-1,k)= Grid(ni+1,nj-1,k)+(double(nIm(ti,tj,k))-double(Grid(ni+1,nj-1,k)))/27;
                end
                if (nj+1<=wide) && (ni-1>=1)
                    Grid(ni-1,nj+1,k)= Grid(ni-1,nj+1,k)+(double(nIm(ti,tj,k))-double(Grid(ni-1,nj+1,k)))/27;
                end
                if (nj-1>=1) && (ni-1>=1)
                    Grid(ni-1,nj-1,k)= Grid(ni-1,nj-1,k)+(double(nIm(ti,tj,k))-double(Grid(ni-1,nj-1,k)))/27;
                end
            end
        end
        
        for ti =1:m
            for tj = 1:n
                 nowdist = 0;
                ni = randi(wide);
                nj = randi(wide);
                ClusterIm(ti,tj)=(ni-1)*wide+nj;
                for k = 1:d
                    nowdist = nowdist +(double(nIm(ti,tj,k))-double(Grid(ni,nj,k))).^2;
                end
                for i = 1:wide
                    for j = 1:wide
                        tdist = 0;
                        for k = 1:d
                            tdist = tdist +(double(nIm(ti,tj,k))-double(Grid(i,j,k))).^2;
                        end
                        if (tdist < nowdist) && ((i-1)*wide+j<=NumClusts)
                            ClusterIm(ti,tj)=(i-1)*wide+j;
                            nowdist = tdist;
                        end
                    end
                end
            end
        end
    end
end

