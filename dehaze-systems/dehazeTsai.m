function [predImage, predT, predA, state] = dehazeTsai(img, ~, state)
    AThresh = 0.85;
    skyThresh = 0.6;
    skyLim = 0.3;
    ADiffThresh = 5.0/255;
    AUpdate = 0.95;
    r = 15;
    eps = 0.001;
    t0 = 0.1;
    omega = 0.95;
    gray = rgb2gray(img);
    
   
    mask = gray>=AThresh;
    % Sorry, but I hate matlab sometimes. You should be able to just do
    % img(mask,:) to get an x by 3 array. But no. You have to do it all
    % individually because you cant do partial logical indexing. 
    % If you do mask = repmat(mask,1,1,3) and index using that,
    % it removes the third dimension structure of the resultant array.
    % Seriously matlab, what are you doing. predA should be as easy as:
    % predA = mean(img(mask,:),[1 2]); with some 'keep dimensions' flag.
    % But nooooooo, it has to be overcomplicated.
    At = ones(1,1,3);
    red = img(:,:,1);
    blue = img(:,:,2);
    green = img(:,:,3);
    At(:,:,1) = mean(red(mask));
    At(:,:,2) = mean(blue(mask));
    At(:,:,3) = mean(green(mask));
    
    
    if isempty(state)
        predA = At;
    else
        Ap = state.Ap;
        
        if norm(Ap-At)>ADiffThresh
            predA = Ap;
        else
            predA = Ap.*AUpdate + (1-AUpdate).*At;
        end
    end
    
    state.Ap = predA;
    
    [m, n, ~] = size(img);

    se = strel('square',r);
    repAtmosphere = repmat(predA, m, n);
    %normed = img ./ repAtmosphere;
    darkChannel = min(img,[],3);
    darkChannel = imerode(darkChannel,se);
    transEst = 1 - omega * darkChannel;
    
    
    avgDenom = windowSumFilter(ones(m, n), r*4);
    mean_g = windowSumFilter(gray, r*4) ./ avgDenom;
    corr_gg = windowSumFilter(gray .* gray, r*4) ./ avgDenom;
    var_g = corr_gg - mean_g .* mean_g;
    mask = (mean_g-var_g)>skyThresh; % create a skymask from the grayscale guide image
    transEst(mask & transEst<skyLim) = skyLim;
%     figure;imagesc(mask);
    % There is literally no way they got the mask shown in their paper
    % if they were using a radius of 28 (4*r_min) whcih they used for single image
    % and it's even less plausible if they used 60. Unless used r_min
    % for these filters and r_min*4 for the following ones, but that would
    % completely destroy the functionality of the guided filter.
    % The more I work on this the more I really dislike this paper.

    mean_t = windowSumFilter(transEst, r*4) ./ avgDenom;
    corr_gt = windowSumFilter(gray .* transEst, r*4) ./ avgDenom;
    cov_gt = corr_gt - mean_g .* mean_t;
    a = cov_gt ./ (var_g + eps);
    b = mean_t - a .* mean_g;
    mean_a = windowSumFilter(a, r*4) ./ avgDenom;
    mean_b = windowSumFilter(b, r*4) ./ avgDenom;

    x = mean_a .* gray + mean_b;
    predT = reshape(x, m, n);
    predT = max(predT, t0);
    
    maxTransmission = repmat(predT, [1, 1, 3]);
    
    predImage = ((img - repAtmosphere) ./ maxTransmission) + repAtmosphere;
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