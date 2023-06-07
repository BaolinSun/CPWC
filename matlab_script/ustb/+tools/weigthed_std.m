function out_std = weigthed_std(D,W,dim)
%weigthed_std Multidimensional weigthed std
% 
%   Usage: out_std = weigthed_std(D,W,dim)
%
%   Computes the weigthed standard deviation of a multidimensional matrix
%   D for a multidimensional weigth W along dimension dim
 
    out_std=sqrt(tools.weigthed_var(D,W,dim));
  
end

