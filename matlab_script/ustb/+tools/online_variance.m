classdef online_variance < handle
%ONLINE_VARIANCE Welford algorithm for computing online variance
%
% B. P. Welford (1962), "Note on a method for calculating corrected 
% sums of squares and products". Technometrics 4(3):419–420.
    
    %% public properties 
    properties (SetAccess = public)
        n
        mean
        m2
    end
    
    %% dependent properties
    properties  (Dependent)
        var          % variance of the samples
        std          % std of the samples
    end
    
    %% constructor
    methods
        function h=online_variance()
            h.n=0;
        end
    end
    
    %% add samples
    methods
        function add(h,value)
            if h.n==0
                h.mean=zeros(size(value));
                h.m2=zeros(size(value));
            end
            h.n=h.n+1;
            
            delta = value - h.mean;
            h.mean = h.mean + delta/h.n;
            
            delta2 = value - h.mean;
            h.m2 = h.m2 + delta.*delta2;
        end
    end
    
    %% get methods
    methods
        function value=get.var(h)
            value=h.m2/(h.n-1);
        end
        function value=get.std(h)
            value=sqrt(h.var);
        end
        
    end
end

