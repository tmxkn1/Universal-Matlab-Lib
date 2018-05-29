classdef line < matlab.mixin.Copyable
    properties
        Point geometry.point
        DirectionVector (1,:) double
    end
    properties (SetAccess = immutable)
        Type geometry.type
    end
    properties (Dependent)
        Dimension
    end
    
    methods
        function obj = line(point, directionvector)
            
            if isa(point, 'geometry.point')
                obj.Point = copy(point);
            else
                obj.Point = geometry.point(point);
            end
            
            obj.DirectionVector = directionvector/norm(directionvector);
            
            obj.Type = geometry.type.Line;
        end
        
        function val = get.Dimension(obj)
            val = numel(obj.DirectionVector);
        end            
    end
end