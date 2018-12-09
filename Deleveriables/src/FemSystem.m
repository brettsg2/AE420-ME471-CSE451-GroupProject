classdef FemSystem < handle
    %FEMSYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        K
        R
        D
        Elements
    end
    
    methods
        function obj = FemSystem(k, r, d, ele)
            obj.K = k;
            obj.R = r;
            obj.D = d;
            obj.Elements = ele;
        end
    end
    
end

