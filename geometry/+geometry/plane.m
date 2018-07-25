classdef plane
    properties
        PointObj geometry.point
        NormalVector (1,:) double
    end
    properties (SetAccess = immutable)
        Type geometry.type
    end
    properties (Dependent)
        Dimension
        Point
    end
    
    methods
        function obj = plane(point, normalvector)
            obj.Type = geometry.type.Plane;
            if nargin == 0
                return
            end
            
            if isa(point, 'geometry.point')
                obj.PointObj = copy(point);
            else
                obj.PointObj = geometry.point(point);
            end
            
            obj.NormalVector = normalvector/norm(normalvector);
        end
        
        function val = get.Dimension(obj)
            val = numel(obj.NormalVector);
        end      
        
        function val = get.Point(obj)
            val = obj.PointObj.Value;
        end       
    end
end