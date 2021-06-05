%
% 17/11/2009
% Author: J.-P. Tarel
% LCPC-INRETS-IFSTTAR copyright
% completed in 06/03/2013, corrected in 08/03/2013
% Modified by Benjamin Tilbury 03/06/2021

function [predImage, predT, predA, state] = dehazeTarel(img, knowns, state)

sv = 15;
p = 0.95;
balance = 0.5;
gfactor = 1.3;

vh = 175; % Assumptions of horizon line
rcalib = 1.65*knowns.K(1);
minvd = 50.0;

[dimy, dimx, ncol]=size(img);

if (ncol==1) 
	w=img; 
	nbo=img;
end
if (ncol==3) 
	if (balance==0.0) % global white balance on clear pixels
		w=min(img,[],3); 
		ival = quantile(w(:),[.99])
		[rind,cind]=find(w>=ival);
		sel(:,1)=img(sub2ind(size(img),rind,cind,ones(size(rind))));
		sel(:,2)=img(sub2ind(size(img),rind,cind,2*ones(size(rind))));
		sel(:,3)=img(sub2ind(size(img),rind,cind,3*ones(size(rind))));
		white=mean(sel,1);
		white=white./max(white)
		img(:,:,1)=img(:,:,1)./white(1);
		img(:,:,2)=img(:,:,2)./white(2);
		img(:,:,3)=img(:,:,3)./white(3);
	end
	if (balance>0.0) % local white balance
		fo(:,:,1)=medfilt2(img(:,:,1), [sv, sv], 'symmetric');
		fo(:,:,2)=medfilt2(img(:,:,2), [sv, sv], 'symmetric');
		fo(:,:,3)=medfilt2(img(:,:,3), [sv, sv], 'symmetric');
		nbfo=mean(fo,3);
		fo(:,:,1)=(fo(:,:,1)./nbfo).^balance;
		fo(:,:,2)=(fo(:,:,2)./nbfo).^balance;
		fo(:,:,3)=(fo(:,:,3)./nbfo).^balance;
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
	nbo=mean(img,3);
end

% compute saturation bound
wm=medfilt2(w, [sv, sv], 'symmetric');
sw=abs(w-wm);
swm=medfilt2(sw, [sv, sv], 'symmetric');
b=wm-swm;
%compute planar assumption bound
c=ones(size(b));
for v=1:dimy
	ci=1-exp((log(0.05)*rcalib)/(minvd*max(v-vh,0)));
	c(v,:)=c(v,:)*ci;
end
% combining bounds
b=min(b,c);
% infered athmospheric veil respecting w and b bounds
v=p*max(min(w,b),0);

predT = 1.0-v;

% restoration with inverse Koschmieder's law
factor=1.0./(1.0-v);
r=zeros(size(img));
if (ncol==1) 
	r=(img-v).*factor; 
	nbr=r;
end
if (ncol==3) 
	r(:,:,1)= (img(:,:,1)-v).*factor; 
	r(:,:,2)= (img(:,:,2)-v).*factor; 
	r(:,:,3)= (img(:,:,3)-v).*factor; 
	% restore original light colors
	if (balance==0.0) 
		r(:,:,1)=r(:,:,1).*white(1);
		r(:,:,2)=r(:,:,2).*white(2);
		r(:,:,3)=r(:,:,3).*white(3);
	end
	if (balance>0.0) 
		r(:,:,1)=r(:,:,1).*fo(:,:,1);
		r(:,:,2)=r(:,:,2).*fo(:,:,2);
		r(:,:,3)=r(:,:,3).*fo(:,:,3);
	end
	nbr=mean(r,3);
end
assignin("base","nbr",nbr);
% final gamma correction 
u=r.^(1.0/gfactor);

% final tone mapping for a gray level between O and 1
mnbu=max(u(:));
predImage=u./(1.0+(1.0-1.0/mnbu)*u);

predA = 0;

end





