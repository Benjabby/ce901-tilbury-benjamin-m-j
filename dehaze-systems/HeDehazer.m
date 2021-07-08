classdef HeDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        % Defaults from paper
        omega       = 0.95;
        win_size    = 15;
        r           = 20;
        eps         = 0.001;
        t0          = 0.1;
    end
    
    methods
        function self = HeDehazer(omega, win_size, r, eps, t0)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(omega), self.omega = omega; end
            if nargin>1 && ~isempty(win_size), self.win_size = win_size; end
            if nargin>2 && ~isempty(r), self.r = r; end
            if nargin>3 && ~isempty(eps), self.eps = eps; end
            if nargin>4 && ~isempty(t0), self.t0 = t0; end
            
        end
        
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img)
            [m, n, ~] = size(img);

            ATic = tic;

            darkChannel = min(img,[],3);
            se = strel('square', self.win_size);
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

            x = BaseDehazer.guidedFilter(rgb2gray(img), transEst, self.r, self.eps);

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

