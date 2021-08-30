classdef (Sealed) TilburyDehazer < BaseDehazer
% TILBURYDEHAZER
% method:       The three different TilburyDehazer methods.
%               Either 'local', 'global' or 'opening'. Defaults to 'local'.
%               If any of the following parameters are not specified or empty
%               the default will be the optimised value for the 'method'.
%
% alpha:        Peak factor for the SEF or ASELF. A higher value results in the filters acting more like an erosion/dilation filter.
%  
% gamma:        Transmission gamma correction factor. (Note; NOT image gamma correction)
%
% omega:        Proportion of haze to remove from the image, leaving some in for visual cues
%
% t0:           Lower bound for transmission map
%
% r:            Radius of the SEF or ASELF filter
%
% rA:           Radius used in calculating the atmospheric light.
%
% rM:           With method=='local':   Radius of morphological gradient
%               With method=='opening': Radius of smooth maximum
%
% alphaM:       With method=='opening': Peak factor for the smooth-dilation.
%
% minvd:        Assumed minimum possible visibility.
%
% rc:           Height of the camera from the ground
%
% vh:           Assumption of horizon line location.

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = true;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        method = 'local';               % The three different TilburyDehazer methods. Either 'local', 'global' or 'opening'. Defaults to 'local'. If any of the following parameters are not specified or empty the default will be the optimised value for the 'method'.
        alpha  = 10.1994823198406;      % Peak factor for the SEF or ASELF. A higher value results in the filters acting more like an erosion/dilation filter.
        gamma  = 1.28476110308647;      % Transmission gamma correction factor.
        omega  = 0.851134429960658;     % Proportion of haze to remove from the image, leaving some in for visual cues
        t0     = 5.41478542049196e-05;  % Lower bound for final transmission map
        r      = 42;                    % Radius of the SEF or ASELF filter
        rA     = 39;                    % Radius used in calculating the atmospheric light.
        rM     = 40;                    % With method=='local':   Radius of morphological gradient. With method=='opening': Radius of smooth maximum
        alphaM = 0;                     % With method=='opening': Peak factor for the smooth-dilation.

		minvd       = 50.0;             % Assumed minimum possible visibility.
		rc          = 1.65;             % Height of the camera from the ground
        vh          = 175;              % Assumption of horizon line location.
    end
    
    methods (Static)
        
        function system = initial
            %% Returns a TilburyDehazer with the initial (unoptimised) parameters.
            system = TilburyDehazer('local',20,1.25,0.95,0.1,15,30,15);
        end
        
        function system = bestLocal
            %% Returns a Local type TilburyDehazer with the best found parameters
            system = TilburyDehazer('local');
        end
        
        function system = bestGlobal
            %% Returns a Global type TilburyDehazer with the best found parameters
            system = TilburyDehazer('global');
        end
        
        function system = bestOpening
            %% Returns an Opening type TilburyDehazer with the best found parameters
            system = TilburyDehazer('opening');
        end
        

        function params = bestLocalParams
            params.method = 'local';
            params.alpha  = 10.1994823198406;
            params.gamma  = 1.28476110308647;
            params.omega  = 0.851134429960658;
            params.t0     = 5.41478542049196e-05;
            params.r      = 42;
            params.rA     = 39;
            params.rM     = 40;
            params.alphaM = 0; % not used
            params.Name   = "TilburyLocal";
        end
           
        function params = bestGlobalParams
            params.method = 'global';
            params.alpha  = 7.21237201480562;
            params.gamma  = 1.28137736777996;
            params.omega  = 0.909942155725948;
            params.t0     = 0.160068556079595;
            params.r      = 40;
            params.rA     = 43;
            params.rM     = 1; % not used
            params.alphaM = 0; % not used
            params.Name = "TilburyGlobal";
        end
        
        function params = bestOpeningParams
            params.method = 'opening';
            params.alpha  = 23.8610032629807;
            params.gamma  = 1.27706489799881;
            params.omega  = 0.957621936853399;
            params.t0     = 0.183191485781434;
            params.r      = 37;
            params.rA     = 46;
            params.rM     = 41;
            params.alphaM = 5.59288358085337;
            params.Name   = "TilburyOpening";
        end
     
    end
    
    methods
        function self = TilburyDehazer(method, alpha, gamma, omega, t0, r, rA, rM, alphaM, minvd, rc, vh)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(method)
                self.method = method;
            end
            
            if self.method=="global"
                proto = TilburyDehazer.bestGlobalParams;
            elseif self.method=="opening"
                proto = TilburyDehazer.bestOpeningParams;
            else
                proto = TilburyDehazer.bestLocalParams;
            end
            
            self.rename(proto.Name);
            
            if nargin>1 && ~isempty(alpha), self.alpha = alpha;     else, self.alpha = proto.alpha; end
            if nargin>2 && ~isempty(gamma), self.gamma = gamma;     else, self.gamma = proto.gamma; end
            if nargin>3 && ~isempty(omega), self.omega = omega;     else, self.omega = proto.omega; end
            if nargin>4 && ~isempty(t0), self.t0 = t0;              else, self.t0 = proto.t0; end
            if nargin>5 && ~isempty(r), self.r = r;                 else, self.r = proto.r;  end
            if nargin>6 && ~isempty(rA), self.rA = rA;              else, self.rA = proto.rA;  end
            if nargin>7 && ~isempty(rM), self.rM = rM;              else, self.rM = proto.rM;  end
            if nargin>8 && ~isempty(alphaM), self.alphaM = alphaM;  else, self.alphaM = proto.alphaM;  end
            if nargin>9 && ~isempty(minvd), self.minvd = minvd; end
            if nargin>10 && ~isempty(rc), self.rc = rc; end
            if nargin>11 && ~isempty(vh), self.vh = vh;  end
            
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
        %% Global Method (Fig.17)
            darkChannel = TilburyDehazer.SEF(minChannel,self.r,-self.alpha);
        end
        
        
        function darkChannel = seLocal(self, minChannel)
        %% Local Method (Fig.16)
        
            % Calculate morphological gradient as smoothing mask
            w = self.rM*2 + 1; % radius to window witdh;
            se = strel('square',w);
            ero = imerode(minChannel,se);
            dil = imdilate(minChannel,se);
            grad = dil-ero;
            
            darkChannel = TilburyDehazer.ASELF(minChannel,self.r,-self.alpha,grad);
        end

        
        function darkChannel = seOpen(self, minChannel)
        %% Opening Method (Fig.18)
        
            % Smooth minimum
            darkChannel = TilburyDehazer.SEF(minChannel,self.r,-self.alpha);
            
            % Smooth maximum
            darkChannel = TilburyDehazer.SEF(darkChannel,self.rM,self.alphaM);
        end

    end
    
    methods (Static, Access=private)
        
        %% SEF (Algorithm 1)
        function out = SEF(img, radius, alpha)
            expD = exp(alpha.*img);
            denom = BaseDehazer.windowSumFilter(expD,radius);
            out = BaseDehazer.windowSumFilter(expD.*img,radius)./denom;
        end
        
        %% ASELF (Algorithm 2)
        function out = ASELF(img, radius, alpha, map)
            expD = exp(map.*alpha.*img);
            denom = BaseDehazer.windowSumFilter(expD,radius);
            out = BaseDehazer.windowSumFilter(expD.*img,radius)./denom;
        end
    end
end