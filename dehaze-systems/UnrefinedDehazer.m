classdef (Sealed) UnrefinedDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        % Values equivalent to unoptimised TilburyDehazer params.
        omega       = 0.95;
        winSize     = 31; % Note that r=15 == winSize=31, or generally r=x == winSize=x*2+1.
        winSizeA    = 30;
        t0          = 0.01;
    end
    
    methods
        function self = UnrefinedDehazer(omega, winSize,  t0)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(omega), self.omega = omega; end
            if nargin>1 && ~isempty(winSize), self.winSize = winSize; end
            if nargin>2 && ~isempty(t0), self.t0 = t0; end
            
        end
        
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img, ~)
            [m, n, ~] = size(img);

            ATic = tic;

            se = strel('square', self.winSizeA);
            darkChannel = min(img,[],3);
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
            
            se = strel('square', self.winSize);
            darkChannel = min(normed,[],3);
            darkChannel = imerode(darkChannel,se);

            predT = 1 - self.omega * darkChannel;

            clampedTransmission = max(predT, self.t0);
            clampedTransmission = repmat(clampedTransmission, [1, 1, 3]);

            predImage = ((img - repAtmosphere) ./ clampedTransmission) + repAtmosphere;
            timeImage = toc(predTic);
            
            predImage = BaseDehazer.clip(predImage);
            predT = BaseDehazer.clip(predT);
        end
        
    end
end

