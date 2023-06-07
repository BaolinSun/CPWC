classdef dimension < int32
%DIMENSION   Enumeration for dimension types. To see the options available write "dimension." and press <TAB>.
%
%   See also PROCESS

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Date: 2017/09/22 $
    
   enumeration
      none(0)
      receive(1)
      transmit(2)
      both(3)
   end
end
