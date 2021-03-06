function [predImage, predT, predA, timeImage, timeA, state] = dehazeZhu(img, ~, state)
    r = 15;
    beta = 0.95;
    eps = 0.001;
    
    t0 = 0.1;
    t1 = 0.9;
    
    [m, n, ~] = size(img);
    
    initialTic = tic;
    [dR, dP] = calVSMap(img, r);
    refineDR = fastGuidedFilterColor(img, dP, r, eps, r/4);
    
    tR = exp(-beta*refineDR);
    %tP = exp(-beta*dP);
    
    initialTime = toc(initialTic);
    
    ATic = tic;
    predA = estA(img, dR);
    predA = reshape(predA, [1,1,3]);
    timeA = toc(ATic) + initialTime;
    
    predTic = tic;
    repAtmosphere = repmat(predA, m, n);
    predT = max(tR, t0);
    clampedTransmission = repmat(predT, [1, 1, 3]);
    predImage = ((img - repAtmosphere) ./ clampedTransmission) + repAtmosphere;
    timeImage = toc(predTic) + initialTime;
end

function [outputRegion, outputPixel] = calVSMap(I, r)
if nargin < 2
    r = 15;
end
hsvI = rgb2hsv(I);
s = hsvI(:,:,2);
v = hsvI(:,:,3);
sigma = 0.041337;
sigmaMat = normrnd(0, sigma, size(I, 1), size(I, 2));

% Here is another parameters setting
%
% output = 0.1893 + 1.0267*v  - 1.2966*s; 

output = 0.121779 + 0.959710*v  - 0.780245*s + sigmaMat;
outputPixel = output;
se = strel('square',r);
output = imerode(output,se);
% output = ordfilt2(output, 1, ones(r,r), 'symmetric');
outputRegion = output;
%imwrite(outputRegion, 'E:\vsFeature.jpg');
end

function q = fastGuidedFilterColor(I, p, r, eps, s)
%   GUIDEDFILTER_COLOR   O(1) time implementation of guided filter using a color image as the guidance.
%
%   - guidance image: I (should be a color (RGB) image)
%   - filtering input image: p (should be a gray-scale/single channel image)
%   - local window radius: r
%   - regularization parameter: eps
%   - subsampling ratio: s (try s = r/4 to s=r)

I_sub = imresize(I, 1/s, 'nearest'); % NN is often enough
p_sub = imresize(p, 1/s, 'nearest');
r_sub = r / s; % make sure this is an integer

[hei, wid] = size(p_sub);
N = windowSumFilter(ones(hei, wid), r_sub); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.

mean_I_r = windowSumFilter(I_sub(:, :, 1), r_sub) ./ N;
mean_I_g = windowSumFilter(I_sub(:, :, 2), r_sub) ./ N;
mean_I_b = windowSumFilter(I_sub(:, :, 3), r_sub) ./ N;

mean_p = windowSumFilter(p_sub, r_sub) ./ N;

mean_Ip_r = windowSumFilter(I_sub(:, :, 1).*p_sub, r_sub) ./ N;
mean_Ip_g = windowSumFilter(I_sub(:, :, 2).*p_sub, r_sub) ./ N;
mean_Ip_b = windowSumFilter(I_sub(:, :, 3).*p_sub, r_sub) ./ N;

% covariance of (I, p) in each local patch.
cov_Ip_r = mean_Ip_r - mean_I_r .* mean_p;
cov_Ip_g = mean_Ip_g - mean_I_g .* mean_p;
cov_Ip_b = mean_Ip_b - mean_I_b .* mean_p;

% variance of I in each local patch: the matrix Sigma in Eqn (14).
% Note the variance in each local patch is a 3x3 symmetric matrix:
%           rr, rg, rb
%   Sigma = rg, gg, gb
%           rb, gb, bb
var_I_rr = windowSumFilter(I_sub(:, :, 1).*I_sub(:, :, 1), r_sub) ./ N - mean_I_r .*  mean_I_r; 
var_I_rg = windowSumFilter(I_sub(:, :, 1).*I_sub(:, :, 2), r_sub) ./ N - mean_I_r .*  mean_I_g; 
var_I_rb = windowSumFilter(I_sub(:, :, 1).*I_sub(:, :, 3), r_sub) ./ N - mean_I_r .*  mean_I_b; 
var_I_gg = windowSumFilter(I_sub(:, :, 2).*I_sub(:, :, 2), r_sub) ./ N - mean_I_g .*  mean_I_g; 
var_I_gb = windowSumFilter(I_sub(:, :, 2).*I_sub(:, :, 3), r_sub) ./ N - mean_I_g .*  mean_I_b; 
var_I_bb = windowSumFilter(I_sub(:, :, 3).*I_sub(:, :, 3), r_sub) ./ N - mean_I_b .*  mean_I_b; 

a = zeros(hei, wid, 3);
for y=1:hei
    for x=1:wid        
        Sigma = [var_I_rr(y, x), var_I_rg(y, x), var_I_rb(y, x);
            var_I_rg(y, x), var_I_gg(y, x), var_I_gb(y, x);
            var_I_rb(y, x), var_I_gb(y, x), var_I_bb(y, x)];
        
        cov_Ip = [cov_Ip_r(y, x), cov_Ip_g(y, x), cov_Ip_b(y, x)];        
        
        a(y, x, :) = cov_Ip  / (Sigma + eps * eye(3)); % very inefficient. Replace this in your C++ code.
    end
end

b = mean_p - a(:, :, 1) .* mean_I_r - a(:, :, 2) .* mean_I_g - a(:, :, 3) .* mean_I_b; % Eqn. (15) in the paper;

mean_a(:, :, 1) = windowSumFilter(a(:, :, 1), r_sub)./N;
mean_a(:, :, 2) = windowSumFilter(a(:, :, 2), r_sub)./N;
mean_a(:, :, 3) = windowSumFilter(a(:, :, 3), r_sub)./N;
mean_b = windowSumFilter(b, r_sub)./N;

mean_a = imresize(mean_a, [size(I, 1), size(I, 2)], 'bilinear'); % bilinear is recommended
mean_b = imresize(mean_b, [size(I, 1), size(I, 2)], 'bilinear');
q = mean_a(:, :, 1) .* I(:, :, 1) + mean_a(:, :, 2) .* I(:, :, 2) + mean_a(:, :, 3) .* I(:, :, 3) + mean_b;
end

function sumImg = windowSumFilter(image, r)

% sum_img(x, y) = = sum(sum(image(x-r:x+r, y-r:y+r)));

[h, w] = size(image);
sumImg = zeros(size(image));

% Y axis
im_cum = cumsum(image, 1);

sumImg(1:r+1, :) = im_cum(1+r:2*r+1, :);
sumImg(r+2:h-r, :) = im_cum(2*r+2:h, :) - im_cum(1:h-2*r-1, :);
sumImg(h-r+1:h, :) = repmat(im_cum(h, :), [r, 1]) - im_cum(h-2*r:h-r-1, :);

% X axis
im_cum = cumsum(sumImg, 2);

sumImg(:, 1:r+1) = im_cum(:, 1+r:2*r+1);
sumImg(:, r+2:w-r) = im_cum(:, 2*r+2:w) - im_cum(:, 1:w-2*r-1);
sumImg(:, w-r+1:w) = repmat(im_cum(:, w), [1, r]) - im_cum(:, w-2*r:w-r-1);

end

function [ A ] = estA(img, Jdark)
    [m, n, ~] = size(img);
    nPixels = m * n;
    nSearchPixels = ceil(nPixels * 0.001);
    darkVec = reshape(Jdark, nPixels, 1);
    imageVec = reshape(img, nPixels, 3);
    [~,ind] = maxk(darkVec, nSearchPixels);
    topA = imageVec(ind,:);
    norms = sqrt(sum(topA.^2,2));
    [~, top] = max(norms);
    A = topA(top,:);
end