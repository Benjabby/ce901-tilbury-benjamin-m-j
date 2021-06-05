function [predImage, predT, predA, state] = dehazeAncuti(img, ~, state)
    inputA = whiteBalance(img);
    inputB = enhanceContrast(img);
    
    weightA = luminanceWeightmap(inputA).*chromaticWeightmap(inputA).*saliencyWeightmap(inputA);
    weightB = luminanceWeightmap(inputB).*chromaticWeightmap(inputB).*saliencyWeightmap(inputB);
    
    weightSum = weightA + weightB;
    weightA = weightA ./ weightSum;
    weightB = weightB ./ weightSum;

    gaussianPyramidA = genPyr(weightA,'gauss',5);
    gaussianPyramidB = genPyr(weightB,'gauss',5);
    
    fusedPyramid = cell(1,5);
    for i = 1 : 5
        tempImg = [];
        for j = 1 : size(img,3)
            laplacianPyramidA = genPyr(inputA(:,:,j),'laplace',5); %Generating Laplacian Pyramid for derrived inputs
            laplacianPyramidB = genPyr(inputB(:,:,j),'laplace',5);
            rowSize = min([size(laplacianPyramidA{i},1),size(laplacianPyramidB{i},1),size(gaussianPyramidA{i},1),size(gaussianPyramidB{i},1)]);
            columnSize = min([size(laplacianPyramidA{i},2),size(laplacianPyramidB{i},2),size(gaussianPyramidA{i},2),size(gaussianPyramidB{i},2)]);
            tempImg(:,:,j) = laplacianPyramidA{i}(1:rowSize , 1:columnSize) .* gaussianPyramidA{i}(1:rowSize, 1:columnSize) + laplacianPyramidB{i}(1:rowSize, 1:columnSize) .* gaussianPyramidB{i}(1:rowSize, 1:columnSize);
        end
        fusedPyramid{i} = tempImg;
    end

    predImage = pyrReconstruct(fusedPyramid);
    predT = 0;
    predA = 0;

end

function [ whiteBalanced ] = whiteBalance( im )
% Using Grayworld assumtion color balancing.....

R_avg = mean(mean(im(:,:,1)));           % Getting the average Of R ,G, B components
G_avg = mean(mean(im(:,:,2)));
B_avg = mean(mean(im(:,:,3)));
RGB_avg = [R_avg G_avg B_avg];

gray_value = (R_avg + G_avg + B_avg)/3;  % By Grey world, avg color of the whole image is gray
scaleValue = gray_value./RGB_avg;       % By Grey world, scale value=gray / avg of each color component

whiteBalanced(:,:,1) = scaleValue(1) * im(:,:,1);   % R,G,B components of the new white balanced new image
whiteBalanced(:,:,2) = scaleValue(2) * im(:,:,2);
whiteBalanced(:,:,3) = scaleValue(3) * im(:,:,3);

end

function [ enhancedIm ] = enhanceContrast( im )

luminance= im(:,:,1)*0.299+im(:,:,2)*0.587+im(:,:,3)*0.114;
avgLuminance =mean(luminance(:));
gamma = 2 * (0.5 + avgLuminance);
% Second input= g * (original image components - avg of the image components)
enhancedIm(:,:,1) = gamma * (im(:,:,1) - avgLuminance);
enhancedIm(:,:,2) = gamma * (im(:,:,2) - avgLuminance);
enhancedIm(:,:,3) = gamma * (im(:,:,3) - avgLuminance);
end

function [ lumWeightmap ] = luminanceWeightmap( im )
L = mean(im,3);
lumWeightmap = sqrt((1/3) * (im(:, :, 1) - L).^2 + (im(:, :, 2) - L).^2 + (im(:, :, 3) - L).^2);
end

function [ chromaticWeightmap ] = chromaticWeightmap( im )

hsvImage = rgb2hsv(im);         % Convert the image to HSV space
saturationValue = hsvImage(:,:,2);    % find saturation

% Chromatic weight map = Distnc b/w saturation value and max sat range..
saturationMax = 1 ;
sigma = .3 ;
chromaticWeightmap = exp( -1 * (((saturationValue - saturationMax).^2) / (2*(sigma^2))) );

end

function [ saliencyWeightmap ] = saliencyWeightmap( im )

if(size(im,3) > 1)
    imGray = rgb2gray(im);
else
    imGray = im;
end

kernel_1D = (1/16) * [1, 4, 6, 4, 1];
kernel_2D = kron(kernel_1D, kernel_1D');

I_mean = mean(imGray(:));

I_Whc = conv2(imGray, kernel_2D, 'same');

saliencyWeightmap = abs(I_Whc - I_mean);

end

function [ pyr ] = genPyr( img, type, level )
%GENPYR generate Gaussian or Laplacian pyramid
%   PYR = GENPYR(A,TYPE,LEVEL) A is the input image, 
%	can be gray or rgb, will be forced to double. 
%	TYPE can be 'gauss' or 'laplace'.
%	PYR is a 1*LEVEL cell array.
% Yan Ke @ THUEE, xjed09@gmail.com

pyr = cell(1,level);
pyr{1} = im2double(img);
for p = 2:level
	pyr{p} = pyr_reduce(pyr{p-1});
end
if strcmp(type,'gauss'), return; end

for p = level-1:-1:1 % adjust the image size
	osz = size(pyr{p+1})*2-1;
	pyr{p} = pyr{p}(1:osz(1),1:osz(2),:);
end

for p = 1:level-1
	pyr{p} = pyr{p}-pyr_expand(pyr{p+1});
end

end

function [ imgout ] = pyr_expand( img )
%PYR_EXPAND  Image pyramid expansion
%   B = PYR_EXPAND( A )  If A is M-by-N, then the size of B 
%	is (2*M-1)-by-(2*N-1). Support gray or rgb image.
%	B will be transformed to double class.
%	Results the same w/ MATLAB func impyramid.
% Yan Ke @ THUEE, xjed09@gmail.com

kw = 5; % default kernel width
cw = .375; % kernel centre weight, same as MATLAB func impyramid. 0.6 in the Paper
ker1d = [.25-cw/2 .25 cw .25 .25-cw/2];
kernel = kron(ker1d,ker1d')*4;

% expand [a] to [A00 A01;A10 A11] with 4 kernels
ker00 = kernel(1:2:kw,1:2:kw); % 3*3
ker01 = kernel(1:2:kw,2:2:kw); % 3*2
ker10 = kernel(2:2:kw,1:2:kw); % 2*3
ker11 = kernel(2:2:kw,2:2:kw); % 2*2

img = im2double(img);
sz = size(img(:,:,1));
osz = sz*2-1;
imgout = zeros(osz(1),osz(2),size(img,3));

for p = 1:size(img,3)
	img1 = img(:,:,p);
	img1ph = padarray(img1,[0 1],'replicate','both'); % horizontally padded
	img1pv = padarray(img1,[1 0],'replicate','both'); % horizontally padded
	
	img00 = imfilter(img1,ker00,'replicate','same');
	img01 = conv2(img1pv,ker01,'valid'); % imfilter doesn't support 'valid'
	img10 = conv2(img1ph,ker10,'valid');
	img11 = conv2(img1,ker11,'valid');
	
	imgout(1:2:osz(1),1:2:osz(2),p) = img00;
	imgout(2:2:osz(1),1:2:osz(2),p) = img10;
	imgout(1:2:osz(1),2:2:osz(2),p) = img01;
	imgout(2:2:osz(1),2:2:osz(2),p) = img11;
end
end

function [ imgout ] = pyr_reduce( img )
%PYR_REDUCE  Image pyramid reduction
%   B = PYR_REDUCE( A )  If A is M-by-N, then the size of B 
%	is ceil(M/2)-by-ceil(N/2). Support gray or rgb image.
%	B will be transformed to double class.
%	Results the same w/ MATLAB func impyramid.
% Yan Ke @ THUEE, xjed09@gmail.com

kernelWidth = 5; % default
cw = .375; % kernel centre weight, same as MATLAB func impyramid. 0.6 in the Paper
ker1d = [.25-cw/2 .25 cw .25 .25-cw/2];
kernel = kron(ker1d,ker1d');

img = im2double(img);
sz = size(img);
imgout = [];

for p = 1:size(img,3)
	img1 = img(:,:,p);
	imgFiltered = imfilter(img1,kernel,'replicate','same');
	imgout(:,:,p) = imgFiltered(1:2:sz(1),1:2:sz(2));
end
end

function [ img ] = pyrReconstruct( pyr )
%PYRRECONSTRUCT Uses a Laplacian pyramid to reconstruct a image
%   IMG = PYRRECONSTRUCT(PYR) PYR should be a 1*level cell array containing
%   the pyramid, SIZE(PYR{i}) = SIZE(PYR{i-1})*2-1
%		Yan Ke @ THUEE, xjed09@gmail.com

for p = length(pyr)-1:-1:1
	pyr{p} = pyr{p}+pyr_expand(pyr{p+1});
end
img = pyr{1};

end

