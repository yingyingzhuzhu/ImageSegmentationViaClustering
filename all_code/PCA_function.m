function HyperImPCA = PCA_function(Im, num)
%%% num means reduced dimension
    [m,n,p] = size(Im);
    Im_reshape = reshape(Im,m*n,p);
    Covar = Im_reshape' * Im_reshape;
    [V,D] = eig(Covar);
    vect = V(:,end - num + 1: end);
    HyperImPCA_temp = Im_reshape*vect;
    HyperImPCA = reshape(HyperImPCA_temp,m,n,num);
end