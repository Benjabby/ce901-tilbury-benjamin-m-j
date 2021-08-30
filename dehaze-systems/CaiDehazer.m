classdef (Sealed) CaiDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 1;
        PredictsA   = true;
        PredictsT   = true;
    end
    properties (SetAccess = private)
        % Defaults from paper & original code (paper params supersede code params)
        omega       = 0.7;      % Proportion of haze to remove
        t0          = 0.1;      % Lower bound for transmission map
        r           = 20;       % Spatial radius of the spatio-temporal MRF
        eps         = 0.0001;   % Epsilon (error preventing small value) for spatio-temporal MRF
        f           = 1;        % Number of adjacent frames to use (temporal radius of the spatio-temporal MRF)
        rho         = 0.1;      % Atmospheric light update factor
        rA          = 15;       % Radius for the erosion filter used in estimating atmospheric light
        downsample  = 1;        % Downsample ratio.
        
        lambda;                 % Temporal window weights. Derived from f.
    end
    
    methods
        function self = CaiDehazer(omega, r, eps, rA, t0, f, rho, downsample)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(omega), self.omega = omega; end
            if nargin>1 && ~isempty(r), self.r = r; end
            if nargin>2 && ~isempty(eps), self.eps = eps; end
            if nargin>3 && ~isempty(rA), self.rA = rA; end
            if nargin>4 && ~isempty(t0), self.t0 = t0; end
            if nargin>5 && ~isempty(f), self.f = f; end
            if nargin>6 && ~isempty(rho), self.rho = rho; end
            if nargin>7 && ~isempty(downsample), self.downsample = downsample; end
            
            self.lambda=hann(self.f*2+3);
            self.lambda(self.lambda==0)=[];
            self.lambda=self.lambda./sum(self.lambda);
            self.lambda=reshape(self.lambda,1,1,numel(self.lambda));
            
        end
        
        function self = newSequence(self, knowns)
            if nargin>1
                newSequence@BaseDehazer(self, knowns);
            else
                newSequence@BaseDehazer(self);
            end
            self.SequenceState.frame = 0;
        end
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img, ~)
            
            minFunc = @(block_struct) min(block_struct.data,[],[1 2]);
            
            if self.SequenceState.frame == 0
                minValue = blockproc(img,[self.rA self.rA], minFunc);
                minValue = reshape(minValue, [], 3);
                L = mean(minValue,2);
                predA = minValue(L==max(L),:);
                predA = mean(predA,1);
                predA = reshape(predA,1,1,3);
                self.SequenceState.A = predA;
                
                gray = rgb2gray(img);
                D = min(img,[],3);
                D = imresize(D, 1/self.downsample, 'nearest');

                % Note that box filter is used here instead of the standardized sumImage filter divided by a denominator, 
                % as the original function used here (A compiled C++ function) gave different results. However the 
                % the results for matlab's imboxfilt were consistent with the original. Using this removes the optimized
                % C++ code to make it slightly more of a fair speed comparison
                uD = imboxfilt(D,self.r*2+1); 
                
                self.SequenceState.dBuff = cat(3,D,D,D);
                self.SequenceState.udBuff = cat(3,uD,uD,uD);
                self.SequenceState.prevGray = gray;
                self.SequenceState.prevFrame = img;

                % Cannot do anything first frame
                timeA = [];
                timeImage = [];
                predT = [];
                predImage = [];
            else
                predTic = tic;
                pGray  = self.SequenceState.prevGray;
                pFrame = self.SequenceState.prevFrame;

                if isempty(img) % Because of the frame delay, the evaluator will run this again with an empty image after all images have been processed
                    D = self.SequenceState.dBuff(:,:,3);
                    uD = self.SequenceState.udBuff(:,:,3);
                else
                    D = min(img,[],3);
                    D = imresize(D, 1/self.downsample, 'nearest');        
                    uD = imboxfilt(D,self.r*2+1); % NOTE: As mentioned before, this function is used instead of using BaseDehazer.windowSumFunction / denominator because this handles edges differently. The difference in speed is very small.
                    self.SequenceState.prevGray = rgb2gray(img);
                    self.SequenceState.prevFrame = img;
                end
                
                self.SequenceState.dBuff(:,:,1)=[];       
                self.SequenceState.dBuff = cat(3,self.SequenceState.dBuff,D);
                self.SequenceState.udBuff(:,:,1)=[];     
                self.SequenceState.udBuff = cat(3,self.SequenceState.udBuff,uD);

                refined = self.STMRF(pGray);

                predT = 1 - self.omega.*(refined./self.SequenceState.A);
                clampedTransmission = max(self.t0,predT);
                predImage = (pFrame-(self.SequenceState.A.*(1-clampedTransmission)))./clampedTransmission;

                timeImage = toc(predTic);
                
                ATic = tic;
                if isempty(img)
                    predA = self.SequenceState.A;
                else
                    minValue = blockproc(img,[self.rA self.rA], minFunc);
                    minValue = reshape(minValue, [], 3);
                    L = mean(minValue,2);
                    A = minValue(L==max(L),:);
                    A = mean(A,1);
                    A = reshape(A,1,1,3);
                    predA = self.rho*A+(1-self.rho)*self.SequenceState.A;
                end
                timeA = toc(ATic);

                self.SequenceState.A = predA;
            end
            
            self.SequenceState.frame = self.SequenceState.frame + 1;
            
            predImage = BaseDehazer.clip(predImage);
            predT = BaseDehazer.clip(predT);
        end
    end
    
    methods (Access = private)
        function [ refined, w, b ] = STMRF(self, img)
            winSize = self.r*2+1;
            V = imresize(img, 1/self.downsample, 'nearest');
            uV = imboxfilt(V,winSize);
            uVV = imboxfilt(V.*V, winSize);
            uVDT = imboxfilt(self.SequenceState.dBuff.*V, winSize);

            numerator = sum((uVDT-(uV.*self.SequenceState.udBuff)).*self.lambda,3);
            denominator = uVV - uV.^2;

            w = numerator ./ (denominator + self.eps);
            b = sum(self.SequenceState.udBuff.*self.lambda,3)-w .* uV;

            w = imresize(w, [size(img, 1), size(img, 2)], 'bilinear');
            b = imresize(b, [size(img, 1), size(img, 2)], 'bilinear');

            refined = w .* img + b;

        end
        
    end
end

