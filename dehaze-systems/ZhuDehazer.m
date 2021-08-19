classdef (Sealed) ZhuDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        % Defaults from paper & code
        beta        = 1;
        r           = 15;
        eps         = 0.001;
        t0          = 0.1;
        t1          = 0.9;
        theta       = [0.121779 0.959710 -0.780245]; 
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
        
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img, ~)

            [m, n, ~] = size(img);

            initialTic = tic;
            [dR, dP] = self.calVSMap(img);
			% There's potentially an error in the authors original code. The paper specified refining the depth map after performing the patch minimum ('dR' in code, 'd_r' in paper)
			% however the author's code used the pixel depth map ('dP' in code, 'd' in paper).
			% In any case, following the paper, dR is used here.
			% Interestingly, this might have already been spotted by the authors of the DehazeNet paper, without them explicitly saying so.
			% In that paper, unlike the other methods they compare against, they use two values for the results of this method.
			% The first set of results they say is from the author's code (presumably the same code this was adapted from)
			% The second set of results they say is from their implementation of this method, and their implementation has much better results.
			% I can't help but speculate that they also spotted the error and corrected but decided to keep the original results in their comparison as well for some reason.
            refineDR = BaseDehazer.fastGuidedFilterColor(img, dR, self.r, self.eps, self.r/4); 

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

