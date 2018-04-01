function [ClusterIm_Spectral,CCIm_Spectral] = MySpectral3(Im,ImType,Type,NumClusts,k)

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

num = 2;
[m,n,p] = size(Im);
Im_subsample = Im(1:num:end,1:num:end,:);
[m_num,n_num,p_num] = size(Im_subsample);

Im_reshape = reshape(Im_subsample,m_num*n_num,p_num);


M = [ones(m_num*n_num,1),Im_reshape];

M = double(M)';

k1 = 4;
Type1 = 2;
sigma = 1;

W = SimGraph_NearestNeighbors(M, k1, Type1, sigma);

fprintf('Similarity Grapht Done\n');

k2 = k;
Type2 = 1;
[C, L, U] = SpectralClustering(W, k2, Type2);

A = full(C);
ClusterIm_Spectral = reshape(A,m_num,n_num);


% CCIm_Spectral_temp = zeros(m*n/num/num,k);
% 
% for i = 1 : k
%     CCIm_Spectral_temp(find(ClusterIm_Spectral == i),i) = 1;
% end
% 
% CCIm_Spectral_temp = reshape(CCIm_Spectral_temp,m/num,n/num,k);

switch Type
    case ('RGB')
        CCIm_Spectral = getCCIm(ClusterIm_Spectral);%CCIm_Spectral_temp;
    case ('Hyper')
        CCIm_Spectral = [];
end
% toc
end


function W = SimGraph_NearestNeighbors(M, k, Type, sigma)
% SIMGRAPH_NEARESTNEIGHBORS Returns kNN similarity graph
%   Returns adjacency matrix for an k-Nearest Neighbors 
%   similarity graph
%
%   'M' - A d-by-n matrix containing n d-dimensional data points
%   'k' - Number of neighbors
%   'Type' - Type if kNN Graph
%      1 - Normal
%      2 - Mutual
%   'sigma' - Parameter for Gaussian similarity function. Set
%      this to 0 for an unweighted graph. Default is 1.
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

if nargin < 3
   ME = MException('InvalidCall:NotEnoughArguments', ...
       'Function called with too few arguments');
   throw(ME);
end
% 
if ~any(Type == (1:2))
   ME = MException('InvalidCall:UnknownType', ...
       'Unknown similarity graph type');
   throw(ME);
end

n = size(M, 2);

% Preallocate memory
indi = zeros(1, k * n);
indj = zeros(1, k * n);
inds = zeros(1, k * n);

for ii = 1:n
    % Compute i-th column of distance matrix
    dist = distEuclidean(repmat(M(:, ii), 1, n), M);
    
    % Sort row by distance
    [s, O] = sort(dist, 'ascend');
    
    % Save indices and value of the k 
    indi(1, (ii-1)*k+1:ii*k) = ii;
    indj(1, (ii-1)*k+1:ii*k) = O(1:k);
    inds(1, (ii-1)*k+1:ii*k) = s(1:k);
%     ii
end

% Create sparse matrix
W = sparse(indi, indj, inds, n, n);
% 
clear indi indj inds dist s O;

% Construct either normal or mutual graph
if Type == 1
    % Normal
    W = max(W, W');
else
    % Mutual
    W = min(W, W');
end

if nargin < 4 || isempty(sigma)
    sigma = 1;
end

% Unweighted graph
if sigma == 0
    W = (W ~= 0);
    
% Gaussian similarity function
elseif isnumeric(sigma)
    W = spfun(@(W) (simGaussian(W, sigma)), W);
    
else
    ME = MException('InvalidArgument:NotANumber', ...
        'Parameter epsilon is not numeric');
    throw(ME);
end

end

function [C, L, U] = SpectralClustering(W, k, Type)
%SPECTRALCLUSTERING Executes spectral clustering algorithm
%   Executes the spectral clustering algorithm defined by
%   Type on the adjacency matrix W and returns the k cluster
%   indicator vectors as columns in C.
%   If L and U are also called, the (normalized) Laplacian and
%   eigenvectors will also be returned.
%
%   'W' - Adjacency matrix, needs to be square
%   'k' - Number of clusters to look for
%   'Type' - Defines the type of spectral clustering algorithm
%            that should be used. Choices are:
%      1 - Unnormalized
%      2 - Normalized according to Shi and Malik (2000)
%      3 - Normalized according to Jordan and Weiss (2002)
%
%   References:
%   - Ulrike von Luxburg, "A Tutorial on Spectral Clustering", 
%     Statistics and Computing 17 (4), 2007
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

% calculate degree matrix
degs = sum(W, 2);
D    = sparse(1:size(W, 1), 1:size(W, 2), degs);

% compute unnormalized Laplacian
L = D - W;

% compute normalized Laplacian if needed
switch Type
    case 2
        % avoid dividing by zero
        degs(degs == 0) = eps;
        % calculate inverse of D
        D = spdiags(1./degs, 0, size(D, 1), size(D, 2));
        
        % calculate normalized Laplacian
        L = D * L;
    case 3
        % avoid dividing by zero
        degs(degs == 0) = eps;
        % calculate D^(-1/2)
        D = spdiags(1./(degs.^0.5), 0, size(D, 1), size(D, 2));
        
        % calculate normalized Laplacian
        L = D * L * D;
end

% compute the eigenvectors corresponding to the k smallest
% eigenvalues
diff   = eps;
[U, ~] = eigs(L, k, diff);

% in case of the Jordan-Weiss algorithm, we need to normalize
% the eigenvectors row-wise
if Type == 3
    U = bsxfun(@rdivide, U, sqrt(sum(U.^2, 2)));
end

% now use the k-means algorithm to cluster U row-wise
% C will be a n-by-1 matrix containing the cluster number for
% each data point
C = kmeans(U, k, 'start', 'cluster', ...
                 'EmptyAction', 'singleton');
             
% now convert C to a n-by-k matrix containing the k indicator
% vectors as columns
% C = sparse(1:size(D, 1), C, 1);

end