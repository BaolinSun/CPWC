function out_mean = weigthed_mean(D,W,dim)
%weigthed_mean Multidimensional weigthed mean
% 
%   Usage: out_std = weigthed_mean(D,W,dim)
%
%   Computes the weigthed mean of a multidimensional matrix
%   D for a multidimensional weigth W along dimension dim
 
    s=bsxfun(@times,W,D);
    s=sum(s,dim);
    s0=sum(W,dim);
    out_mean=bsxfun(@times,s,1./s0);
    
end

