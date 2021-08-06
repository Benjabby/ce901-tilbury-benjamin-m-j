function [points, points2] = projectPoints(depth, K, R)
    %invK = inv(K);
    [X, Y] = meshgrid(1:size(depth,2),1:size(depth,1));
    coords = cat(3,X,Y,ones(size(X)));
    
    
%     T = P(1:3,4);
%     K = P(1:3,1:3);
    
    % A
    coords = permute(coords,[3 1 2]);
    coords = reshape(coords,3,[]);
    points = (K\coords)' .* depth(:);
    
    % B
    coords(3,:) = depth(:);
    points2 = (K\coords)';
    
    % BT
%     coords = cat(3,X,Y,ones(size(X)));
%     tt = repmat(reshape(T',[1,1,3]),size(depth,1),size(depth,2));
%     coords = (coords-tt) .* depth;
%     coords = permute(coords,[3 1 2]);
%     points2 = (K\reshape(coords,3,[]))';
%     
%     points = points - T';
    
%     points = R*permute(cat(2,points,ones([size(points,1) 1])),[2 1]);
end