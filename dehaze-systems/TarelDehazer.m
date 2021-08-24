classdef (Sealed) TarelDehazer < BaseDehazer

    properties (Constant)
        FrameDelay  = 0;
        PredictsA   = false;
        PredictsT   = true;
    end
    
    properties (SetAccess = private)
        % Defaults from paper except vh & rc
        sv          = 11;
        p           = 0.95;
        balance     = -1;
        gfactor     = 1;
        rc          = 1.65; % Multiplier based on height of KITTI cameras from the ground
        minvd       = 50.0;
        vh          = 175;  % Assumption of horizon line for KITTI data
         
    end
    
    methods
        function self = TarelDehazer(sv, p, balance, gfactor, rc, minvd, vh)
            self = self@BaseDehazer;
            
            if nargin>0 && ~isempty(sv), self.sv = sv; end
            if nargin>1 && ~isempty(p), self.p = p; end
            if nargin>2 && ~isempty(balance), self.balance = balance; end
            if nargin>3 && ~isempty(gfactor), self.gfactor = gfactor; end
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
                [rind,cind]=find(w>=ival);
                sel(:,1)=img(sub2ind(size(img),rind,cind,ones(size(rind))));
                sel(:,2)=img(sub2ind(size(img),rind,cind,2*ones(size(rind))));
                sel(:,3)=img(sub2ind(size(img),rind,cind,3*ones(size(rind))));
                white=mean(sel,1);
                white=white./max(white);
                img(:,:,1)=img(:,:,1)./white(1);
                img(:,:,2)=img(:,:,2)./white(2);
                img(:,:,3)=img(:,:,3)./white(3);
            elseif (self.balance>0.0) % local white balance
                fo(:,:,1)=medfilt2(img(:,:,1), [self.sv, self.sv], 'symmetric');
                fo(:,:,2)=medfilt2(img(:,:,2), [self.sv, self.sv], 'symmetric');
                fo(:,:,3)=medfilt2(img(:,:,3), [self.sv, self.sv], 'symmetric');
                nbfo=mean(fo,3);
                fo(:,:,1)=(fo(:,:,1)./nbfo).^self.balance;
                fo(:,:,2)=(fo(:,:,2)./nbfo).^self.balance;
                fo(:,:,3)=(fo(:,:,3)./nbfo).^self.balance;
                nbfo=mean(fo,3);
                fo(:,:,1)=fo(:,:,1)./nbfo;
                fo(:,:,2)=fo(:,:,2)./nbfo;
                fo(:,:,3)=fo(:,:,3)./nbfo;
                img(:,:,1)=img(:,:,1)./fo(:,:,1);
                img(:,:,2)=img(:,:,2)./fo(:,:,2);
                img(:,:,3)=img(:,:,3)./fo(:,:,3);
            end
            % compute photometric bound
            w=min(img,[],3); 
            %nbo=mean(img,3);

            % compute saturation bound
            wm=medfilt2(w, [self.sv, self.sv], 'symmetric');
            sw=abs(w-wm);
            swm=medfilt2(sw, [self.sv, self.sv], 'symmetric');
            b=wm-swm;
            %compute planar assumption bound	
			c = 1-repmat(exp((log(0.05)*self.SequenceState.rcalib)./(self.minvd*max((1:m)-self.vh,0)))',1,n);
            % combining bounds
            b=min(b,c);
            % infered athmospheric veil respecting w and b bounds
            v=self.p*max(min(w,b),0);

            predT = 1.0-v;

            % restoration with inverse Koschmieder's law
            factor=1.0./(1.0-v);
            r=zeros(size(img));
            if (ncol==1) 
                r=(img-v).*factor; 
                %nbr=r;
            end
            if (ncol==3) 
                r(:,:,1)= (img(:,:,1)-v).*factor; 
                r(:,:,2)= (img(:,:,2)-v).*factor; 
                r(:,:,3)= (img(:,:,3)-v).*factor; 
                % restore original light colors
                if (self.balance==0.0) 
                    r(:,:,1)=r(:,:,1).*white(1);
                    r(:,:,2)=r(:,:,2).*white(2);
                    r(:,:,3)=r(:,:,3).*white(3);
                end
                if (self.balance>0.0) 
                    r(:,:,1)=r(:,:,1).*fo(:,:,1);
                    r(:,:,2)=r(:,:,2).*fo(:,:,2);
                    r(:,:,3)=r(:,:,3).*fo(:,:,3);
                end
                %nbr=mean(r,3);
            end
            
            % final gamma correction 
            u=r.^(1.0/self.gfactor);

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

