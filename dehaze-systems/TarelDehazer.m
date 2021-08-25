classdef (Sealed) TarelDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = false;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        % Defaults from paper except vh & rc
        winSize     = 11;
        omega       = 0.95;
        balance     = -1;
        gamma       = 1;
        rc          = 1.65; % Multiplier based on height of KITTI cameras from the ground
        minvd       = 50.0;
        vh          = 175;  % Assumption of horizon line for KITTI data
         
    end
    
    methods
        function self = TarelDehazer(winSize, omega, balance, gamma, rc, minvd, vh)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(winSize), self.winSize = winSize; end
            if nargin>1 && ~isempty(omega), self.omega = omega; end
            if nargin>2 && ~isempty(balance), self.balance = balance; end
            if nargin>3 && ~isempty(gamma), self.gamma = gamma; end
            if nargin>4 && ~isempty(rc), self.rc = rc; end
            if nargin>5 && ~isempty(mindvd), self.minvd = minvd; end
            if nargin>6 && ~isempty(vh), self.vh = vh; end
            
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
            
            [m, n, ncol]=size(img);

            predTic = tic;
            
            if (self.balance==0.0) % global white balance on clear pixels
                w=min(img,[],3); 
                ival = quantile(w(:),.99);
                mask = w>=ival;
                mask = repmat(mask, 1,1,3);
                sel = reshape(img(mask),[],3);
                white = mean(sel,1);
                white = white./max(white);
                white = reshape(white, [1 1 3]);
                img = img./white;
            elseif (self.balance>0.0) % local white balance
                fo(:,:,1)=medfilt2(img(:,:,1), [self.winSize, self.winSize], 'symmetric');
                fo(:,:,2)=medfilt2(img(:,:,2), [self.winSize, self.winSize], 'symmetric');
                fo(:,:,3)=medfilt2(img(:,:,3), [self.winSize, self.winSize], 'symmetric');
                nbfo=mean(fo,3);
                fo = (fo./nbfo).^self.balance;
                nbfo=mean(fo,3);
                fo = fo./nbfo;
                img = img./fo;
            end
            % compute photometric bound
            w=min(img,[],3); 
            %nbo=mean(img,3);

            % compute saturation bound
            wm=medfilt2(w, [self.winSize, self.winSize], 'symmetric');
            sw=abs(w-wm);
            swm=medfilt2(sw, [self.winSize, self.winSize], 'symmetric');
            b=wm-swm;
            %compute planar assumption bound	
			c = 1-repmat(exp((log(0.05)*self.SequenceState.rcalib)./(self.minvd*max((1:m)-self.vh,0)))',1,n);
            % combining bounds
            b=min(b,c);
            % infered athmospheric veil respecting w and b bounds
            v=self.omega*max(min(w,b),0);

            predT = 1.0-v;

            % restoration with inverse Koschmieder's law
            factor=1.0./(1.0-v);
            recovered=zeros(size(img));
            if (ncol==1) 
                recovered=(img-v).*factor; 
                %nbr=r;
            end
            if (ncol==3) 
                recovered = (img-v).*factor; 
                % restore original light colors
                if (self.balance==0.0) 
                    recovered = recovered.*white;
                end
                if (self.balance>0.0) 
                    recovered = recovered.*fo;
                end
                %nbr=mean(r,3);
            end
            
            % final gamma correction 
            u=recovered.^(1.0/self.gamma);

            % final tone mapping for a gray level between O and 1
            mnbu=max(u(:));
            predImage=u./(1.0+(1.0-1.0/mnbu)*u);

            timeImage = toc(predTic);
            timeA = [];
            predA = [];
            
            predImage = BaseDehazer.clip(predImage);
            predT = BaseDehazer.clip(predT);
        end
        
    end
end

