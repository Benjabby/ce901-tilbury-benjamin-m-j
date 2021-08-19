classdef (Sealed) TilburyDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
		
    end
    
    properties (SetAccess = private)
        
        %% Original values
%         method      = 'masked';
%         alpha       = 20;
%         gamma      = 1.25;
%         omega       = 0.95;
%         t0          =  0.1;
%         r           = 15;
%         rA          = 30;
%         rM          = 15;
        %% Best optimised result
        method      = 'local';
        alpha       = 18.1283727963125;
        gamma      = 1.28331033279787;
        omega       = 0.944081856347016;
        t0          = 0.122535149860325;
        r           = 37;
        rA          = 61;
        rM          = 81;
        alphaM      = 0; % Only used for 'opening' method
		%%
		minvd       = 50.0;
		rc          = 1.65; 	  % camera height
        vh          = 175;  % Assumption of horizon line for KITTI data
    end
    
    methods (Static)
        
        function system = initial
            %% Returns a TilburyDehazer with the initial (unoptimised) parameters.
            system = TilburyDehazer('local',20,1.25,0.95,0.1,15,30,15);
        end
        
        function system = bestLocal
            %% Returns a Local type TilburyDehazer with the best found parameters
            system = TilburyDehazer('local',10.1994823198406,1.28476110308647,0.851134429960658,5.41478542049196e-05,42,39,40);
        end
        
        function system = bestOpening
            %% Returns an Opening type TilburyDehazer with the best found parameters
            system = TilburyDehazer('opening',23.8610032629807,1.27706489799881,0.957621936853399,0.183191485781434,37,46,41,5.59288358085337);
        end
        
        function system = bestGlobal
            %% Returns a Global type TilburyDehazer with the best found parameters
            system = TilburyDehazer('global',7.21237201480562,1.28137736777996,0.909942155725948,0.160068556079595,40,43);
        end
    end
    
    methods
        function self = TilburyDehazer(method, alpha, gamma, omega, t0, r, rA, rM, alphaM, minvd, rc, vh)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(method)
                self.method = method;
            end
            
            if nargin>1 && ~isempty(alpha), self.alpha = alpha; end
            if nargin>2 && ~isempty(gamma), self.gamma = gamma; end
            if nargin>3 && ~isempty(omega), self.omega = omega; end
            if nargin>4 && ~isempty(t0), self.t0 = t0; end
            if nargin>5 && ~isempty(r), self.r = r; end
            if nargin>6 && ~isempty(rA), self.rA = rA; end
            if nargin>7 && ~isempty(rM), self.rM = rM; end
            if nargin>8 && ~isempty(alphaM), self.alphaM = alphaM; end
            if nargin>9 && ~isempty(minvd), self.minvd = minvd; end
            if nargin>10 && ~isempty(rc), self.rc = rc; end
            if nargin>11 && ~isempty(vh), self.vh = vh; end
            
        end
        
		function self = newSequence(self, knowns)
            newSequence@BaseDehazer(self, knowns);
            self.SequenceState.rcalib = self.rc*knowns.K(1);
        end
        
        function [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self, img, ~)
            
            if isempty(self.Knowns) || ~isfield(self.Knowns,'K')
                fprintf("This dehazer requires projection matrix from camera calibration before dehazing.\nPlease use newSequence to pass a structure containing the field 'K' which has the 3x4 projection matrix\n");
                predImage = [];
                predT = [];
                predA = [];
                timeImage = [];
                timeA = [];
                return
            end
            
            [m, n, ~] = size(img);
            
            ATic = tic;
        %     predA = smoothAtmLight(img);

            darkChannel = min(img,[],3);
            w = self.rA*2 + 1; % radius to window width
            se = strel('square', w); 
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
            
            minChannel = min(normed,[],3);

            if self.method=="opening"
                darkChannel = self.seOpen(minChannel);
            elseif self.method=="global"
                darkChannel = self.seGlobal(minChannel);
            else
                darkChannel = self.seLocal(minChannel);
            end
			
            darkChannel = darkChannel.^self.gamma;
            
            predT = 1 - self.omega * darkChannel; 

            
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
            
            plane = repmat(exp((log(0.05)*self.SequenceState.rcalib)./(self.minvd*max((1:m)-self.vh,0)))',1,n);

            predT = max(predT,plane);

            maxTransmission = max(predT, self.t0);

            predImage = ((img - repeatedA) ./ maxTransmission) + repeatedA;
            timeImage = toc(predTic);
            
            predImage = BaseDehazer.clip(predImage);
            predT = BaseDehazer.clip(predT);
        end
        
    end
    
    methods (Access = private)
        
        function darkChannel = seGlobal(self, minChannel)
            expD = exp(-self.alpha.*minChannel);

            denom = BaseDehazer.windowSumFilter(expD,self.r);
            darkChannel = BaseDehazer.windowSumFilter(expD.*minChannel,self.r)./denom;
        end
        
        function darkChannel = seLocal(self, minChannel)
            
            w = self.rM*2 + 1; % radius to window witdh;
            se = strel('square',w);
            ero = imerode(minChannel,se);
            dil = imdilate(minChannel,se);

            V = self.alpha*(dil-ero);
            expD = exp(-V.*minChannel);

            denom = BaseDehazer.windowSumFilter(expD,self.r);
            darkChannel = BaseDehazer.windowSumFilter(expD.*minChannel,self.r)./denom;
        end

        function darkChannel = seOpen(self, minChannel)
            
            % Smooth minimum
            expD = exp(-self.alpha.*minChannel);
            denom = BaseDehazer.windowSumFilter(expD,self.r);
            darkChannel = BaseDehazer.windowSumFilter(expD.*minChannel,self.r)./denom;
            
            % Smooth maximum
            expD = exp(self.alphaM.*darkChannel);
            denom = BaseDehazer.windowSumFilter(expD,self.rM);
            darkChannel = BaseDehazer.windowSumFilter(expD.*darkChannel,self.rM)./denom;
        end

    end
end

