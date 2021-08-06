classdef (Sealed) TilburyDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
		
    end
    
    properties (SetAccess = private)
        omega = 0.95;
        alpha = 20;
        lambda = 1.25;
        r = 15;
        t0 = 0.01;
        method = 'default';
		
		minvd = 50.0;
		rc = 1.65; 	  % camera height
        vh          = 175;  % Assumption of horizon line for KITTI data
    end
    
    methods
        function self = TilburyDehazer(method, alpha, lambda, omega, r, t0, minvd, rc, vh)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(method), self.method = method; end
            if nargin>1 && ~isempty(alpha), self.alpha = alpha; end
            if nargin>2 && ~isempty(lambda), self.lambda = lambda; end
            if nargin>3 && ~isempty(omega), self.omega = omega; end
            if nargin>4 && ~isempty(r), self.r = r; end
            if nargin>5 && ~isempty(t0), self.t0 = t0; end
            if nargin>6 && ~isempty(minvd), self.minvd = minvd; end
            if nargin>7 && ~isempty(rc), self.rc = rc; end
            if nargin>8 && ~isempty(vh), self.vh = vh; end
            
        end
        
		function self = newSequence(self, knowns)
            newSequence@BaseDehazer(self, knowns);
            self.SequenceState.rcalib = self.rc*knowns.K(1);
        end
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img, ~)
            
            [m, n, ~] = size(img);
            
            ATic = tic;
        %     predA = smoothAtmLight(img);

            darkChannel = min(img,[],3);
            se = strel('square',self.r);
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
%             if(all(size(predA)~=size(img)))
            repeatedA = repmat(predA, m, n);
%             else
%                 repeatedA = predA;
%             end

            normed = img ./ repeatedA;

            if self.method=="sat"
                predT = self.transSat(normed);
            elseif self.method=="global"
                predT = self.transGlobal(normed);
            else
                predT = self.transDefault(normed);
            end
			
            %% Horizon line estimation
            
%             bigSE = strel('square',self.r*2+1);
%             h = (imdilate(darkChannel,bigSE)-imerode(darkChannel,bigSE))==0;
%             VH = find(max(h,[],2), 1, 'last');
%             
%             if isempty(VH)
%                 if isfield(self.SequenceState,'prevVH')
%                     VH = self.SequenceState.prevVH
%                 else
%                     VH = floor(m/2);
%                 end
%             end
%             
%             self.SequenceState.prevVH = VH;
            %%
            
            plane = repmat(exp((log(0.05)*self.SequenceState.rcalib)./(self.minvd*max([1:m]-self.vh,0)))',1,n);

            predT = max(predT,plane);

            maxTransmission = max(predT, self.t0);

            predImage = ((img - repeatedA) ./ maxTransmission) + repeatedA;
            timeImage = toc(predTic);
            
            predImage = BaseDehazer.clip(predImage);
            predT = BaseDehazer.clip(predT);
        end
        
    end
    
    methods (Access = private)
        
        function t = transGlobal(self, normed)
            darkChannel = min(normed,[],3);

            expD = exp(-self.alpha.*darkChannel);

            denom = BaseDehazer.windowSumFilter(expD,self.r);
            darkChannel = BaseDehazer.windowSumFilter(expD.*darkChannel,self.r)./denom;
            darkChannel = darkChannel.^self.lambda;
            t = 1 - self.omega * darkChannel; 
        end
        
        function t = transDefault(self, normed)
            darkChannel = min(normed,[],3);

            se = strel('square',self.r);
            ero = imerode(darkChannel,se);
            dil = imdilate(darkChannel,se);

            V = self.alpha*(dil-ero);
            expD = exp(-V.*darkChannel);

            denom = BaseDehazer.windowSumFilter(expD,self.r);
            darkChannel = BaseDehazer.windowSumFilter(expD.*darkChannel,self.r)./denom;
            darkChannel = darkChannel.^self.lambda;
            t = 1 - self.omega * darkChannel; 
        end

        function t = transSat(self, normed)
            darkChannel = min(normed,[],3);
            lightChannel = max(normed,[],3);

            se = strel('square',self.r);
            ero = imerode(darkChannel,se);
            dil = imdilate(darkChannel,se);

            sat = (lightChannel-darkChannel)./lightChannel;
            sat = imdilate(sat,se);
            sat = 1-sat;
            darkChannel = min(darkChannel,sat.^self.lambda);

            V = self.alpha*(dil-ero);
            expD = exp(-V.*darkChannel);

            denom = BaseDehazer.windowSumFilter(expD,self.r);
            darkChannel = BaseDehazer.windowSumFilter(expD.*darkChannel,self.r)./denom;
            t = 1 - self.omega * darkChannel;
        end

    end
end

