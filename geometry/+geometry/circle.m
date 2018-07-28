classdef circle
    properties
        CentrePointObj geometry.point
        Radius (1,1) double
        NormalVector (1,:) double
    end
    properties (SetAccess = immutable)
        Type geometry.type
    end
    properties (Dependent)
        Dimension
        Centre
    end
    
    methods
        function obj = circle(centre, radius, normalvector)            
            obj.Type = geometry.type.Circle;
            if nargin == 0
                return
            end
            
            if isa(centre, 'geometry.point')
                obj.CentrePointObj = copy(centre);
            else
                obj.CentrePointObj = geometry.point(centre);
            end
            
            obj.Radius = radius;
            
            if nargin >= 3
                obj.NormalVector = normalvector/norm(normalvector);
            else
                obj.NormalVector = zeros(size(centre));
                obj.NormalVector(end) = 1;
            end
        end
        
        function val = get.Dimension(obj)
            val = numel(obj.NormalVector);
        end
        
        function val = get.Centre(obj)
            val = obj.CentrePointObj.Value;
        end 
    end
end