classdef line < matlab.mixin.Copyable
    properties
        PointObj geometry.point
        DirectionVector (1,:) double
    end
    properties (SetAccess = immutable)
        Type geometry.type
    end
    properties (Dependent)
        Dimension
        Point
    end
    
    methods
        function obj = line(point, directionvector)
            obj.Type = geometry.type.Line;            
            if nargin == 0
                return
            end
            
            if isa(point, 'geometry.point')
                obj.PointObj = copy(point);
            else
                obj.PointObj = geometry.point(point);
            end
            
            obj.DirectionVector = directionvector/norm(directionvector);
            
        end
        
        function val = get.Dimension(obj)
            val = numel(obj.DirectionVector);
        end 
        
        function val = get.Point(obj)
            val = obj.PointObj.Value;
        end 
        
        function re = isOnLine(obj, p, tol)
            if nargin < 3
                tol = 1e-10;
            end
            g = p - obj.PointObj.Value;
            if obj.Dimension == 3
                re = all(cross(g, obj.DirectionVector)<tol);
            else
                idx = find(obj.DirectionVector~=0, 1);
                t = g(idx) / obj.DirectionVector(idx);
                re = all(obj.PointObj.Value + obj.DirectionVector * t - g < tol);
            end
        end
    end
end