function transmission = createTransmission(disparity, visibility, doubleMask)
    depthMap = 1.0./(disparity.*1000);
    beta = 3.912/visibility;
    beta = repelem(beta,3);
    beta = reshape(beta,[1 1 3]);
    transmission = exp(-beta.*depthMap);
    if nargin>2
        transmission = transmission.*skyMask;
    end
end