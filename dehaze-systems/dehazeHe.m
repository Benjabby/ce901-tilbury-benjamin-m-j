function [predImage, predT, predA, state] = dehazeHe(img, ~, state)
    omega = 0.95;
    win_size = 15;
    
    r = 15;
    res = 0.001;
    t0 = 0.1;
    [m, n, ~] = size(img);

    darkChannel = min(img,[],3);
    se = strel('square',win_size);
    darkChannel = imerode(darkChannel,se);
    
    nPixels = m * n;

    nSearchPixels = floor(nPixels * 0.01);
    darkVec = reshape(darkChannel, nPixels, 1);
    imageVec = reshape(img, nPixels, 3);
    [~,ind] = maxk(darkVec, nSearchPixels);
    predA = mean(imageVec(ind,:),1);
    predA = reshape(predA, [1, 1, 3]);
    
    repAtmosphere = repmat(predA, m, n);
    normed = img ./ repAtmosphere;
    darkChannel = min(normed,[],3);
    darkChannel = imerode(darkChannel,se);

    transEst = 1 - omega * darkChannel;

    x = guidedFilter(rgb2gray(img), transEst, r, res);

    predT = reshape(x, m, n);
    predT = max(predT, t0);
    maxTransmission = repmat(predT, [1, 1, 3]);
    
    predImage = ((img - repAtmosphere) ./ maxTransmission) + repAtmosphere;
end


function q = guidedFilter(guide, target, radius, eps)
% Guided Filter implementation from "Fast Guided Filter"
% http://arxiv.org/abs/1505.00996
% Note that this implementation is slower than Matlab's own
% 'imguidedfilter' however other dehazing methods make modifications to the
% filter than cannot easily be done with Matlab's built-in function

[h, w] = size(guide);

avgDenom = windowSumFilter(ones(h, w), radius);

mean_g = windowSumFilter(guide, radius) ./ avgDenom;
mean_t = windowSumFilter(target, radius) ./ avgDenom;

corr_gg = windowSumFilter(guide .* guide, radius) ./ avgDenom;
corr_gt = windowSumFilter(guide .* target, radius) ./ avgDenom;

var_g = corr_gg - mean_g .* mean_g;
cov_gt = corr_gt - mean_g .* mean_t;

a = cov_gt ./ (var_g + eps);
b = mean_t - a .* mean_g;

mean_a = windowSumFilter(a, radius) ./ avgDenom;
mean_b = windowSumFilter(b, radius) ./ avgDenom;

q = mean_a .* guide + mean_b;

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

