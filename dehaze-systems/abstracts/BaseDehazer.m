classdef (Abstract) BaseDehazer < handle & matlab.mixin.Heterogeneous
    %BASEDEHAZER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract = true, Constant)
        FrameDelay(1,1) logical; % Currently only frame delay of 1 is supported for evaluator
        PredictsA(1,1) logical;
        PredictsT(1,1) logical;
    end
    
    properties (Access = public)
        SequenceState;
    end
    
    methods (Abstract)
        [predImage, predT, predA, timeImage, timeA] = dehazeFrame(self,img);
    end
    
    methods
    
        function self = newSequence(self, knowns)
            self.SequenceState = struct();
            self.SequenceState.Knowns = knowns;
        end
   
    end
    
    %% Standardized operations
    methods (Access = protected, Static = true)
        function v = clip(map, minv, maxv)
            if ~exist('minv','var'), minv = 0; end
            if ~exist('maxv','var'), maxv = 1; end
            v = min(max(map,minv),maxv);
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

            a = zeros(h, w, 3);
            for y=1:h
                for x=1:w        
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
    end
end

