classdef vector
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        norm
        data
    end
    
    methods
        function obj = vector(x,y)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.norm = (x^2+y^2)^0.5;
            obj.data = [x,y];
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

