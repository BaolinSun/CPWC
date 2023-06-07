classdef wavefront < int32
%wavefront   Enumeration for wave types. To see the options available write "wavefront." and press <TAB>.
%
%   See also WAVE

%   authors: Alfonso Rodriguez-Molares (alfonso.r.molares@ntnu.no)
%   $Lat updated: 2017/10/13 $
    
   enumeration
      plane(0)
      spherical(1)
      photoacoustic(2)
   end
end
