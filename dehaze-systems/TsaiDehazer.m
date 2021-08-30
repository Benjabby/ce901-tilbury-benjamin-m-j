classdef (Sealed) TsaiDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        % Defaults from paper
        AThresh     = 0.85;     % Brightness threshold for considering a pixel in atmospheric light estimation
        skyThresh   = 0.6;      % Threshold for considering a region as sky
        skyLim      = 0.3;      % Minimum transmission for regions identified as sky.
        ADiffThresh = 5.0/255;  % Threshold for change in atmospheric light brightness to not update
        AUpdate     = 0.95;     % Atmospheric light update factor
        rMin        = 7;        % Erosion filter radius
        rGF         = 28;       % Guided filter radius
        eps         = 0.001;    % Epsilon (error preventing small value) for the guided filter
        t0          = 0.1;      % Lower bound for transmission map
        omega       = 0.95;     % Proportion of haze to remove
    end
    
    methods
        function self = TsaiDehazer(omega, rMin, eps, t0, AThresh, skyThresh, skyLim, ADiffThresh, AUpdate)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(omega), self.omega = omega; end
            if nargin>1 && ~isempty(rMin), self.rMin = rMin; self.rGF = 4*rMin; end
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
            se = strel('square',self.rMin);
            repAtmosphere = repmat(predA, m, n);
            normed = img ./ repAtmosphere;
            darkChannel = min(normed,[],3);
            darkChannel = imerode(darkChannel,se);
            transEst = 1 - self.omega * darkChannel;


            avgDenom = BaseDehazer.windowSumFilter(ones(m, n), self.rGF);
            mean_g = BaseDehazer.windowSumFilter(gray, self.rGF) ./ avgDenom;
            corr_gg = BaseDehazer.windowSumFilter(gray .* gray, self.rGF) ./ avgDenom;
            var_g = corr_gg - mean_g .* mean_g;
            mask = (mean_g-var_g)>self.skyThresh; % create a skymask from the grayscale guide image
            transEst(mask & transEst<self.skyLim) = self.skyLim;

            mean_t = BaseDehazer.windowSumFilter(transEst, self.rGF) ./ avgDenom;
            corr_gt = BaseDehazer.windowSumFilter(gray .* transEst, self.rGF) ./ avgDenom;
            cov_gt = corr_gt - mean_g .* mean_t;
            a = cov_gt ./ (var_g + self.eps);
            b = mean_t - a .* mean_g;
            mean_a = BaseDehazer.windowSumFilter(a, self.rGF) ./ avgDenom;
            mean_b = BaseDehazer.windowSumFilter(b, self.rGF) ./ avgDenom;

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

