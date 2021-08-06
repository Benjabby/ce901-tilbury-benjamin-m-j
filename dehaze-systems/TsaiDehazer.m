classdef (Sealed) TsaiDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        % Defaults from paper
        AThresh = 0.85;
        skyThresh = 0.6;
        skyLim = 0.3;
        ADiffThresh = 5.0/255;
        AUpdate = 0.95;
        r_min = 7;
        r_mean = 28;
        eps = 0.001;
        t0 = 0.1;
        omega = 0.95;
    end
    
    methods
        function self = TsaiDehazer(omega, r_min, eps, t0, AThresh, skyThresh, skyLim, ADiffThresh, AUpdate)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(omega), self.omega = omega; end
            if nargin>1 && ~isempty(r_min), self.r_min = r_min; self.r_mean = 4*r_min; end
            if nargin>2 && ~isempty(eps), self.eps = eps; end
            if nargin>3 && ~isempty(t0), self.t0 = t0; end
            if nargin>4 && ~isempty(AThresh), self.AThresh = AThresh; end
            if nargin>5 && ~isempty(skyThresh), self.skyThresh = skyThresh; end
            if nargin>6 && ~isempty(skyLim), self.skyLim = skyLim; end
            if nargin>7 && ~isempty(ADiffThresh), self.ADiffThresh = ADiffThresh; end
            if nargin>8 && ~isempty(AUpdate), self.AUpdate = AUpdate; end
            
        end
        
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img, ~)
            ATic = tic;
            gray = rgb2gray(img);
            mask = gray>=self.AThresh;
            At = ones(1,1,3);
            red = img(:,:,1);
            blue = img(:,:,2);
            green = img(:,:,3);
            At(:,:,1) = mean(red(mask));
            At(:,:,2) = mean(blue(mask));
            At(:,:,3) = mean(green(mask));
            
            if ~isfield(self.SequenceState,'Ap')
                predA = At;
            else
                Ap = self.SequenceState.Ap;

                if norm(squeeze(Ap-At))>self.ADiffThresh
                    predA = Ap;
                else
                    predA = Ap.*self.AUpdate + (1-self.AUpdate).*At;
                end
            end
            
            timeA = toc(ATic);
            
    
            self.SequenceState.Ap = predA;
            
            [m, n, ~] = size(img);
    
            predTic = tic;
            se = strel('square',self.r_min);
            repAtmosphere = repmat(predA, m, n);
            normed = img ./ repAtmosphere;
            darkChannel = min(normed,[],3);
            darkChannel = imerode(darkChannel,se);
            transEst = 1 - self.omega * darkChannel;


            avgDenom = BaseDehazer.windowSumFilter(ones(m, n), self.r_mean);
            mean_g = BaseDehazer.windowSumFilter(gray, self.r_mean) ./ avgDenom;
            corr_gg = BaseDehazer.windowSumFilter(gray .* gray, self.r_mean) ./ avgDenom;
            var_g = corr_gg - mean_g .* mean_g;
            mask = (mean_g-var_g)>self.skyThresh; % create a skymask from the grayscale guide image
            transEst(mask & transEst<self.skyLim) = self.skyLim;

            mean_t = BaseDehazer.windowSumFilter(transEst, self.r_mean) ./ avgDenom;
            corr_gt = BaseDehazer.windowSumFilter(gray .* transEst, self.r_mean) ./ avgDenom;
            cov_gt = corr_gt - mean_g .* mean_t;
            a = cov_gt ./ (var_g + self.eps);
            b = mean_t - a .* mean_g;
            mean_a = BaseDehazer.windowSumFilter(a, self.r_mean) ./ avgDenom;
            mean_b = BaseDehazer.windowSumFilter(b, self.r_mean) ./ avgDenom;

            x = mean_a .* gray + mean_b;
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

