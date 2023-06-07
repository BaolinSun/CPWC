function out_var = weigthed_var(D,W,dim)
%weigthed_var Multidimensional weigthed variance
% 
%   Usage: out_var = weigthed_var(D,W,dim)
%
%   Computes the weigthed variance of a multidimensional matrix
%   D for a multidimensional weigth W along dimension dim
 
    s=tools.weigthed_mean(D,W,dim);
    s=bsxfun(@plus,D,-s).^2;
    s=bsxfun(@times,W,s);
    s=sum(s,dim);
    s0=sum(W,dim);
    out_var=bsxfun(@times,s,1./s0);
    
end

