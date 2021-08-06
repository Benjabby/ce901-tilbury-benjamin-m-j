global thres;
thres       = 0.02;
spacing     = 0.02;
N           = 1000 ;
K           = 40;
Amin        = [0,0.05,0.1];
Amax        = [1,1,1];
%% Dehazing parameters
gamma       = 1;
lambda      = 0.1;
leave_haze  = 1.06;
trans_min   = 0.1;


[img_ind, points] = rgb2ind(img, N);
[h,w,~] = size(img);
% Remove empty clusters
idx_in_use = unique(img_ind(:));
idx_to_remove = setdiff(0:(size(points,1)-1),idx_in_use);
points(idx_to_remove+1,:) = [];
img_ind_sequential = changem(img_ind,1:length(idx_in_use),idx_in_use);

[points_weight,~] = histcounts(img_ind_sequential(:),size(points,1));
points_weight = points_weight./(h*w);
if ~ismatrix(points), points = reshape(points,[],3); end % verify dim

%% Define arrays of candidate air-light values and angles
angle_list = reshape(linspace(0, pi, K),[],1);
% Use angle_list(1:end-1) since angle_list(end)==pi, which is the same line
% in 2D as since angle_list(1)==0
directions_all = [sin(angle_list(1:end-1)) , cos(angle_list(1:end-1)) ];

% Air-light candidates in each color channel
ArangeR = Amin(1):spacing:Amax(1);
ArangeG = Amin(2):spacing:Amax(2);
ArangeB = Amin(3):spacing:Amax(3);

%% Estimate air-light in each pair of color channels
% Estimate RG
Aall = generate_Avals(ArangeR, ArangeG);


tic;
[~, AvoteRG_orig] = vote_2D_orig(points(:,1:2), points_weight, directions_all, Aall);
toc

tic;
[~, AvoteRG] = vote_2D(points(:,1:2), points_weight, directions_all, Aall);
toc

function Aall = generate_Avals(Avals1, Avals2)
    %Generate a list of air-light candidates of 2-channels, using two lists of
    %values in a single channel each
    %Aall's length is length(Avals1)*length(Avals2)
    Avals1 = reshape(Avals1,[],1);
    Avals2 = reshape(Avals2,[],1);
    A1 = kron(Avals1, ones(length(Avals2),1));
    A2 = kron(ones(length(Avals1),1), Avals2);
    Aall = [A1, A2];
end


function [Aout, Avote2] = vote_2D(points, points_weight, directions_all, Aall)

    global thres;
    
    check = bsxfun(@and,Aall(:,1)>points(:,1)',Aall(:,2)>points(:,2)');
    
    dist = pdist2(Aall,points);
    comp = dist./sqrt(2) + 1;
    
    p1 = -points(:,1) * directions_all(:,2)' + points(:,2) * directions_all(:,1)';
    p1 = reshape(p1,[1 size(p1)]);
    
    p2 = Aall(:, 1)*directions_all(:,2)' - Aall(:, 2)*directions_all(:,1)';
    p2 = reshape(p2, [size(p2,1) 1 size(p2,2)]);
    
    val = abs(p1+p2); %bsxfun(@plus,p1,p2);
    
    accumulator_votes_idx = bsxfun(@and, check, val<2*thres.*comp);
    assignin('base','v',accumulator_votes_idx);
    
%     for i_point = 1:size(points,1)
%         % save time and ignore irelevant points from the get-go
%         idx_to_use = find( (Aall(:, 1) > points(i_point, 1)) & (Aall(:, 2) > points(i_point, 2)));
%         if isempty(idx_to_use), continue; end
% 
%         % calculate distance between all A options and the line defined by
%         % i_point and i_direction. If the distance is smaller than a thres,
%         % increase the cell in accumulator
%         dist1 = sqrt(sum([Aall(idx_to_use, 1)-points(i_point, 1), Aall(idx_to_use, 2)-points(i_point, 2)].^2,2));
%         %dist1 = dist1 - min(dist1);
%         dist1 = dist1./sqrt(2) + 1;
%         
%         for i_direction = 1:n_directions
%             
% 
%             dist =  -points(i_point, 1)*directions_all(i_direction,2) + ...
%                 points(i_point, 2)*directions_all(i_direction,1) + ...
%                 Aall(idx_to_use, 1)*directions_all(i_direction,2) - ...
%                 Aall(idx_to_use, 2)*directions_all(i_direction,1);
%             idx = abs(dist)<2*thres.*dist1;
%             if ~any(idx), continue; end
% 
%             idx_full = idx_to_use(idx);
%             accumulator_votes_idx(idx_full, i_point,i_direction) = true;
%         end
%     end
    accumulator_votes_idx2 = (sum(uint8(accumulator_votes_idx),2))>=2; 
    accumulator_votes_idx = bsxfun(@and, accumulator_votes_idx ,accumulator_votes_idx2);
%     accumulator_unique = zeros(size(Aall,1),1);
    
    points_weight_dist = (5.*exp(-dist)+1).*points_weight; 
    indx = any(accumulator_votes_idx,3);%bsxfun(@and,any(accumulator_votes_idx,3),check);
    accumulator_unique = sum(points_weight_dist.*indx,2);
    
    assignin('base','ac',accumulator_unique);
   
%     for iA = 1:size(Aall,1)
%         idx_to_use = find(Aall(iA, 1) > points(:, 1) & (Aall(iA, 2) > points(:, 2)));
%         points_dist = sqrt((Aall(iA,1) - points(idx_to_use,1)).^2+(Aall(iA,2) - points(idx_to_use,2)).^2);
%         points_weight_dist = points_weight(idx_to_use).*(5.*exp(-reshape(points_dist,1,[]))+1); 
%         accumulator_unique(iA) = sum(points_weight_dist(any(accumulator_votes_idx(iA,idx_to_use,:),3)));
%     end
    [~, Aestimate_idx] = max(accumulator_unique);
    Aout = Aall(Aestimate_idx,:);
    Avote2 = accumulator_unique; 
end

function [Aout, Avote2] = vote_2D_orig(points, points_weight, directions_all, Aall)

    global thres;
    n_directions = size(directions_all,1);
    accumulator_votes_idx = false(size(Aall,1), size(points,1), n_directions);
    for i_point = 1:size(points,1)
        for i_direction = 1:n_directions
             % save time and ignore irelevant points from the get-go
            idx_to_use = find( (Aall(:, 1) > points(i_point, 1)) & (Aall(:, 2) > points(i_point, 2)));
            if isempty(idx_to_use), continue; end

            % calculate distance between all A options and the line defined by
            % i_point and i_direction. If the distance is smaller than a thres,
            % increase the cell in accumulator
            dist1 = sqrt(sum([Aall(idx_to_use, 1)-points(i_point, 1), Aall(idx_to_use, 2)-points(i_point, 2)].^2,2));
            %dist1 = dist1 - min(dist1);
            dist1 = dist1./sqrt(2) + 1;

            dist =  -points(i_point, 1)*directions_all(i_direction,2) + ...
                points(i_point, 2)*directions_all(i_direction,1) + ...
                Aall(idx_to_use, 1)*directions_all(i_direction,2) - ...
                Aall(idx_to_use, 2)*directions_all(i_direction,1);
            idx = abs(dist)<2*thres.*dist1;
            if ~any(idx), continue; end

            idx_full = idx_to_use(idx);
            accumulator_votes_idx(idx_full, i_point,i_direction) = true;
        end
    end
    
    assignin('base','vo',accumulator_votes_idx);
    
    accumulator_votes_idx2 = (sum(uint8(accumulator_votes_idx),2))>=2; 
    accumulator_votes_idx = bsxfun(@and, accumulator_votes_idx ,accumulator_votes_idx2);
    accumulator_unique = zeros(size(Aall,1),1);
    for iA = 1:size(Aall,1)
        idx_to_use = find(Aall(iA, 1) > points(:, 1) & (Aall(iA, 2) > points(:, 2)));
        points_dist = sqrt((Aall(iA,1) - points(idx_to_use,1)).^2+(Aall(iA,2) - points(idx_to_use,2)).^2);
        points_weight_dist = points_weight(idx_to_use).*(5.*exp(-reshape(points_dist,1,[]))+1); 
        accumulator_unique(iA) = sum(points_weight_dist(any(accumulator_votes_idx(iA,idx_to_use,:),3)));
    end
    
    assignin('base','aco',accumulator_unique);
    [~, Aestimate_idx] = max(accumulator_unique);
    Aout = Aall(Aestimate_idx,:);
    Avote2 = accumulator_unique; 
end