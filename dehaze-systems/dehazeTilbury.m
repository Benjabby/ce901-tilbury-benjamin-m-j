function [predImage, predT, predA, timeImage, timeA, state] = dehazeTilbury(img, ~, state)
    omega = 0.95;
    alpha = 20;
    sp = 1.5;
    r = 15;
    t0 = 0.01;
    [m, n, ~] = size(img);

    ATic = tic;
%     predA = smoothAtmLight(img);

    darkChannel = min(img,[],3);
    se = strel('square',r);
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
    if(all(size(predA)~=size(img)))
        repeatedA = repmat(predA, m, n);
    else
        repeatedA = predA;
    end
    normed = img ./ repeatedA;
    darkChannel = min(normed,[],3);
    lightChannel = max(normed,[],3);
    
    se = strel('square',r);
    ero = imerode(darkChannel,se);
    dil = imdilate(darkChannel,se);
    
    sat = (lightChannel-darkChannel)./lightChannel;
    sat = imdilate(sat,se);
    sat = 1-sat;
    darkChannel = min(darkChannel,sat.^sp);
    
    V = alpha*(dil-ero);
    expD = exp(-V.*darkChannel);
%     expD = exp(-alpha*darkChannel);
%     expD = imerode(expD,se);
    
    denom = windowSumFilter(expD,r);
    darkChannel = windowSumFilter(expD.*darkChannel,r)./denom;
%     darkChannel = darkChannel.^1.25;
    predT = 1 - omega * darkChannel;
    %predT = predT.^0.8;
    
    maxTransmission = max(predT, t0);
    
    predImage = ((img - repeatedA) ./ maxTransmission) + repeatedA;
    timeImage = toc(predTic);
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


function A = smoothAtmLight(img, alpha, low)
    if ~exist('alpha','var'), alpha = 60; end
    if ~exist('low','var'), low = 0.1; end
    
    D = min(img,[],3);
    se = strel('square',15);
    D = imerode(D,se);
    
    E = exp(D*alpha);
    A = E./sum(E,[1 2]);
    A = sum(A.*img, [1 2]);
%     assignin('base','tA',A);
%     A =  sum(A .* sum(E,3), [1 2]);
%     [m, n, ~] = size(img);
    A = reshape(A, [1 1 3]);
    %A = A./B;
    
% %     img = srgb2lrgb(img);
%     E = exp(-img*alpha);
%     A = E./sum(E,[1 2]);
%     A = sum(A.*img, [1 2]);
    
% %     E = exp(img*alpha);
% %     B = E./sum(E,3);
% %     B = sum(B.*img, 3);
%     B = max(img,[],3);
% %     B = max(B,low);
%     A = A./B;
% %     A = lrgb2srgb(A);
end
