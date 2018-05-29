classdef plane
    properties
        Point geometry.point
        NormalVector (1,:) double
    end
    properties (SetAccess = immutable)
        Type geometry.type
    end
    properties (Dependent)
        Dimension
    end
    
    methods
        function obj = plane(point, normalvector)
            
            if isa(point, 'geometry.point')
                obj.Point = copy(point);
            else
                obj.Point = geometry.point(point);
            end
            
            obj.NormalVector = normalvector/norm(normalvector);
            obj.Type = geometry.type.Plane;
        end
        
        function val = get.Dimension(obj)
            val = numel(obj.NormalVector);
        end            
    end
end