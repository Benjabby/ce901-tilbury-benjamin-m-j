function q = fastGuidedFilterColor(guide, target, r, eps, s)
    %   GUIDEDFILTER_COLOR   O(1) time implementation of guided filter using a color image as the guidance.
    %
    %   - guidance image: I (should be a color (RGB) image)
    %   - filtering input image: p (should be a gray-scale/single channel image)
    %   - local window radius: r
    %   - regularization parameter: eps
    %   - subsampling ratio: s (try s = r/4 to s=r)

    guide_sub = imresize(guide, 1/s, 'nearest'); % NN is often enough
    target_sub = imresize(target, 1/s, 'nearest');
    r_sub = r / s; % make sure this is an integer

    [h, w] = size(target_sub);
    N = windowSumFilter(ones(h, w), r_sub); % the size of each local patch; N=(2r+1)^2 except for boundary pixels.

    mean_I_r = windowSumFilter(guide_sub(:, :, 1), r_sub) ./ N;
    mean_I_g = windowSumFilter(guide_sub(:, :, 2), r_sub) ./ N;
    mean_I_b = windowSumFilter(guide_sub(:, :, 3), r_sub) ./ N;

    mean_p = windowSumFilter(target_sub, r_sub) ./ N;

    mean_Ip_r = windowSumFilter(guide_sub(:, :, 1).*target_sub, r_sub) ./ N;
    mean_Ip_g = windowSumFilter(guide_sub(:, :, 2).*target_sub, r_sub) ./ N;
    mean_Ip_b = windowSumFilter(guide_sub(:, :, 3).*target_sub, r_sub) ./ N;

    % covariance of (I, p) in each local patch.
    cov_Ip_r = mean_Ip_r - mean_I_r .* mean_p;
    cov_Ip_g = mean_Ip_g - mean_I_g .* mean_p;
    cov_Ip_b = mean_Ip_b - mean_I_b .* mean_p;

    % variance of I in each local patch: the matrix Sigma in Eqn (14).
    % Note the variance in each local patch is a 3x3 symmetric matrix:
    %           rr, rg, rb
    %   Sigma = rg, gg, gb
    %           rb, gb, bb
    var_I_rr = windowSumFilter(guide_sub(:, :, 1).*guide_sub(:, :, 1), r_sub) ./ N - mean_I_r .*  mean_I_r; 
    var_I_rg = windowSumFilter(guide_sub(:, :, 1).*guide_sub(:, :, 2), r_sub) ./ N - mean_I_r .*  mean_I_g; 
    var_I_rb = windowSumFilter(guide_sub(:, :, 1).*guide_sub(:, :, 3), r_sub) ./ N - mean_I_r .*  mean_I_b; 
    var_I_gg = windowSumFilter(guide_sub(:, :, 2).*guide_sub(:, :, 2), r_sub) ./ N - mean_I_g .*  mean_I_g; 
    var_I_gb = windowSumFilter(guide_sub(:, :, 2).*guide_sub(:, :, 3), r_sub) ./ N - mean_I_g .*  mean_I_b; 
    var_I_bb = windowSumFilter(guide_sub(:, :, 3).*guide_sub(:, :, 3), r_sub) ./ N - mean_I_b .*  mean_I_b; 
    
    N = h*w;
    top = cat(2,reshape(var_I_rr,1,1,[]),reshape(var_I_rg,1,1,[]),reshape(var_I_rb,1,1,[]));
    mid = cat(2,reshape(var_I_rg,1,1,[]),reshape(var_I_gg,1,1,[]),reshape(var_I_gb,1,1,[]));
    bot = cat(2,reshape(var_I_rb,1,1,[]),reshape(var_I_gb,1,1,[]),reshape(var_I_bb,1,1,[]));
    sigma = cat(1,top,mid,bot); % I'm sure there's a more efficient way to stack like this but it's not really important
    cov_Ip = cat(1,reshape(cov_Ip_r,1,[],1),reshape(cov_Ip_g,1,[],1),reshape(cov_Ip_b,1,[],1));
    reps = repmat(eps*reshape(eye(3),3,3,1),1,1,N);
    
    M = sigma + reps;
    
    cov_Ip = reshape(cov_Ip, 3*N, []); 
    I = repmat(reshape(1:3*N,3,1,N),[1 3 1]); % m x n x p
    J = repmat(reshape(1:3*N,1,3,N),[3 1 1]); % m x n x p
    A = sparse(I(:),J(:),M(:));
    a = reshape(reshape(A \ cov_Ip, [3 N])',[h w 3]);

    b = mean_p - a(:, :, 1) .* mean_I_r - a(:, :, 2) .* mean_I_g - a(:, :, 3) .* mean_I_b; % Eqn. (15) in the paper;

    mean_a(:, :, 1) = windowSumFilter(a(:, :, 1), r_sub)./N;
    mean_a(:, :, 2) = windowSumFilter(a(:, :, 2), r_sub)./N;
    mean_a(:, :, 3) = windowSumFilter(a(:, :, 3), r_sub)./N;
    mean_b = windowSumFilter(b, r_sub)./N;

    mean_a = imresize(mean_a, [size(guide, 1), size(guide, 2)], 'bilinear'); % bilinear is recommended
    mean_b = imresize(mean_b, [size(guide, 1), size(guide, 2)], 'bilinear');
    q = mean_a(:, :, 1) .* guide(:, :, 1) + mean_a(:, :, 2) .* guide(:, :, 2) + mean_a(:, :, 3) .* guide(:, :, 3) + mean_b;
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