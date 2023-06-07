classdef weighted_incremental_variance < handle
%WEIGHTED_INCREMENTAL_VARIANCE Welford algorithm improved to compute a
%weighted variance incrementally as described by
%
% D. H. D. West (1979). Communications of the ACM, 22, 9, 532-535: "Updating 
% Mean and Variance Estimates: An Improved Method"

    %% public properties 
    properties (SetAccess = public)
        n
        wSum
        wSum2
        mean
        S
    end
    
    %% dependent properties
    properties  (Dependent)
        var          % variance of the samples
        std          % std of the samples
    end
    
    %% constructor
    methods
        function h=weighted_incremental_variance()
            h.n=0;
        end
    end
    
    %% add samples
    methods
        function add(h,weight,value)
            if h.n==0
                h.wSum=zeros(size(weight));
                h.wSum2=zeros(size(weight));
                h.mean=zeros(size(weight));
                h.S=zeros(size(weight));
            end
            h.n=h.n+1;
            
            % update
            h.wSum = h.wSum + weight;
            h.wSum2 = h.wSum2 + weight.*weight;
            meanOld = h.mean;
            if h.wSum>0
                h.mean = meanOld + (weight ./ h.wSum) .* (value - meanOld);
            else
                h.mean = meanOld;
            end
            h.S = h.S + weight.* (value - meanOld) .* (value - h.mean);
        end
    end
    
    %% get methods
    methods
        function value=get.var(h)
            value=h.S./(h.wSum-1);
        end
        function value=get.std(h)
            value=sqrt(h.var);
        end
        
    end
end

