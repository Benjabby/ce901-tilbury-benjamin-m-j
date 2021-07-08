classdef ZhuDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        % Defaults from paper & code
        beta        = 0.95;
        r           = 15;
        eps         = 0.001;
        t0          = 0.1;
        t1          = 0.9;
        theta       = [0.1893 1.0267 -1.2966]; 
        sigma       = 0.041337;
    end
    
    methods
        function self = ZhuDehazer(r, theta, beta, sigma, eps, t0, t1)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(r), self.r = r; end
            if nargin>1 && ~isempty(theta), self.theta = theta; end
            if nargin>2 && ~isempty(beta), self.beta = beta; end
            if nargin>3 && ~isempty(sigma), self.sigma = sigma; end
            if nargin>4 && ~isempty(eps), self.eps = eps; end
            if nargin>5 && ~isempty(t0), self.t0 = t0; end
            if nargin>6 && ~isempty(t1), self.t1 = t1; end
            
        end
        
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img)

            [m, n, ~] = size(img);

            initialTic = tic;
            [dR, dP] = self.calVSMap(img);
            refineDR = BaseDehazer.fastGuidedFilterColor(img, dP, self.r, self.eps, self.r/4);
            %tP = exp(-beta*dP);

            initialTime = toc(initialTic);

            ATic = tic;
            predA = self.estA(img, dR);
            predA = reshape(predA, [1,1,3]);
            timeA = toc(ATic) + initialTime;

            predTic = tic;
            repAtmosphere = repmat(predA, m, n);
            tR = exp(-self.beta*refineDR);
            predT = tR;
            clampedTransmission = min(max(tR, self.t0),self.t1);
            clampedTransmission = repmat(clampedTransmission, [1, 1, 3]);
            predImage = ((img - repAtmosphere) ./ clampedTransmission) + repAtmosphere;
            timeImage = toc(predTic) + initialTime;
            
            predImage = BaseDehazer.clip(predImage);
            predT = BaseDehazer.clip(predT);
        end 
    end
    
    methods (Access = private)
        function [outputRegion, outputPixel] = calVSMap(self, img)

            hsvI = rgb2hsv(img);
            s = hsvI(:,:,2);
            v = hsvI(:,:,3);
            sigmaMat = normrnd(0, self.sigma, size(img, 1), size(img, 2));

            outputPixel = self.theta(1) + v*self.theta(2) + self.theta(3)*s + sigmaMat;
            se = strel('square',self.r);
            outputRegion = imerode(outputPixel,se);
        end
        
        function A = estA(~, img, Jdark)
            [m, n, ~] = size(img);
            nPixels = m * n;
            nSearchPixels = ceil(nPixels * 0.001);
            darkVec = reshape(Jdark, nPixels, 1);
            imageVec = reshape(img, nPixels, 3);
            [~,ind] = maxk(darkVec, nSearchPixels);
            topA = imageVec(ind,:);
            norms = sqrt(sum(topA.^2,2));
            [~, top] = max(norms);
            A = topA(top,:);
        end
    end
end

