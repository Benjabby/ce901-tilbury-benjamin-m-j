classdef LiDehazer < ExternalPythonDehazer
   % Courtesy wrapper 
    
    methods
        function self = LiDehazer(showWindow, postponeStart)
            if nargin<1, showWindow = false; end
            if nargin<2, postponeStart = false; end
            self = self@ExternalPythonDehazer("AOD",showWindow, postponeStart);
        end
    end
end

