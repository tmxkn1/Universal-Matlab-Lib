classdef point < handle
    %POINT defines a point
    
    properties
        Value (1,:) double
    end
    
    properties (SetAccess = immutable)
        Type geometry.type
    end
    
    properties (Dependent)
        Dimension        
    end
    
    methods
        function obj = point(val)
            obj.Type = geometry.type.Point;
            if nargin == 0
                return
            end
            
            obj.Value = val;
        end
        
        function val = get.Dimension(obj)
            val = numel(obj.Value);
        end            
    end
end

