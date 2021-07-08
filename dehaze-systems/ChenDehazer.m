classdef ChenDehazer < ExternalPythonDehazer
   % Courtesy wrapper 

    methods
        function self = ChenDehazer(showWindow, postponeStart)
            if nargin<1, showWindow = false; end
            if nargin<2, postponeStart = false; end
            self = self@ExternalPythonDehazer("PMHLD", showWindow, postponeStart);
        end
    end
end

