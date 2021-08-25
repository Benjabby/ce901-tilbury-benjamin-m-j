classdef (Sealed) HeDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        % Defaults from paper
        omega       = 0.95;
        winSize    = 15;
        rGF           = 20;
        eps         = 0.001;
        t0          = 0.1;
    end
    
    methods
        function self = HeDehazer(omega, winSize, rGF, eps, t0)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(omega), self.omega = omega; end
            if nargin>1 && ~isempty(winSize), self.winSize = winSize; end
            if nargin>2 && ~isempty(rGF), self.rGF = rGF; end
            if nargin>3 && ~isempty(eps), self.eps = eps; end
            if nargin>4 && ~isempty(t0), self.t0 = t0; end
            
        end
        
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img, ~)
            [m, n, ~] = size(img);

            ATic = tic;

            darkChannel = min(img,[],3);
            se = strel('square', self.winSize);
            darkChannel = imerode(darkChannel,se);

            nPixels = m * n;

            nSearchPixels = floor(nPixels * 0.01);
            darkVec = reshape(darkChannel, nPixels, 1);
            imageVec = reshape(img, nPixels, 3);
            [~,ind] = maxk(darkVec, nSearchPixels);
            predA = mean(imageVec(ind,:),1);
            predA = reshape(predA, [1, 1, 3]);
            timeA = toc(ATic);

            predTic = tic;
            repAtmosphere = repmat(predA, m, n);
            normed = img ./ repAtmosphere;
            darkChannel = min(normed,[],3);
            darkChannel = imerode(darkChannel,se);

            transEst = 1 - self.omega * darkChannel;

            x = BaseDehazer.guidedFilter(rgb2gray(img), transEst, self.rGF, self.eps);

            predT = reshape(x, m, n);
            clampedTransmission = max(predT, self.t0);
            clampedTransmission = repmat(clampedTransmission, [1, 1, 3]);

            predImage = ((img - repAtmosphere) ./ clampedTransmission) + repAtmosphere;
            timeImage = toc(predTic);
            
            predImage = BaseDehazer.clip(predImage);
            predT = BaseDehazer.clip(predT);
        end
        
    end
end

