function loss = minmin(img, A, x)

    r = (img-A)./x + A;
    darkChannel = min(r,[],3);

    expD = exp(-80.*darkChannel);

    denom = windowSumFilter(expD,30)+1e-8;
    darkChannel = windowSumFilter(expD.*darkChannel,30)./denom;
    loss = sum(1 - darkChannel,'all'); 
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

