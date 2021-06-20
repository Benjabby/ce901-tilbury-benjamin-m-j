function [predImage, predT, predA, timeImage, timeA, state] = dehazeCai(img, ~, state)
omega  = 0.7;
t0=0.1;
% ST-MRF parameter
D_r = 20;
eps = 0.0001;
f=1;
rho=0.1;
A_r=15;
s = 2;

lambda_t = reshape([0.25 0.5 0.25], [1 1 3]);

if ~exist('state','var') || isempty(fieldnames(state))
    %state.frameDelay = 1; % Going to make this implicit;
    state.frame = 0; 
    state.minFunc = @(block_struct) min(block_struct.data,[],[1 2]);
    %ATic = tic;
    minValue = blockproc(img,[A_r A_r], state.minFunc);
    minValue = reshape(minValue, [], 3);
    L = mean(minValue,2);
    predA = minValue(L==max(L),:);
    predA = mean(predA,1);
    predA = reshape(predA,1,1,3);
    
    %timeA = toc(ATic);
    
    state.A = predA;
    
    %predTic = tic;
    gray = rgb2gray(img);
    D = min(img,[],3);
    D = imresize(D, 1/s, 'nearest');
    uD = imboxfilt(D,D_r*2+1);
    
    state.dBuff = cat(3,D,D,D);
    state.udBuff = cat(3,uD,uD,uD);
    state.prevGray = gray;
    state.prevFrame = img;
    
    %timeImage = toc(predTic);
    
    % Cannot do anything first frame
    timeA = [];
    timeImage = [];
    predT = [];
    predImage = [];
    
else 
    state.frame = state.frame + 1;
    predTic = tic;
    pGray = state.prevGray;
    pFrame = state.prevFrame;
    
    if isempty(img) % Because of the frame delay, the evaluator will run this again with an empty image after all images have been processed
        D = state.dBuff(:,:,3);
        uD = state.udBuff(:,:,3);
    else
        D = min(img,[],3);
        D = imresize(D, 1/s, 'nearest');        
        uD = imboxfilt(D,D_r*2+1);
        state.prevGray = rgb2gray(img);
        state.prevFrame = img;
    end
    
    state.dBuff(:,:,1)=[];       
    state.dBuff = cat(3,state.dBuff,D);
    state.udBuff(:,:,1)=[];     
    state.udBuff = cat(3,state.udBuff,uD);
    
    refined = STMRF(pGray, state.dBuff, state.udBuff, lambda_t, D_r, s, eps);
    
    predT = 1 - omega.*(refined./state.A);
    predT = max(t0,predT);
    predImage = (pFrame-(state.A.*(1-predT)))./predT;
    
    timeImage = toc(predTic);
    
    ATic = tic;
    minValue = blockproc(img,[A_r A_r], state.minFunc);
    minValue = reshape(minValue, [], 3);
    L = mean(minValue,2);
    A = minValue(L==max(L),:);
    A = mean(A,1);
    A = reshape(A,1,1,3);
    
    predA = rho*A+(1-rho)*state.A;
    timeA = toc(ATic);
    
    state.A = predA;
end

end


function [ refined, w, b ] = STMRF(img, dBuff, udBuff, lambda_t, radius, s, eps)
    r = radius*2+1;
    V = imresize(img, 1/s, 'nearest');
    uV = imboxfilt(V,r);
    uVV = imboxfilt(V.*V, r);
    uVD_t = imboxfilt(dBuff.*V, r);

    numerator = sum((uVD_t-(uV.*udBuff)).*lambda_t,3);
    denominator = uVV - uV.^2;

    w = numerator ./ (denominator + eps);
    b = sum(udBuff.*lambda_t,3)-w .* uV;

    w = imresize(w, [size(img, 1), size(img, 2)], 'bilinear');
    b = imResample(b, [size(img, 1), size(img, 2)], 'bilinear');

    refined = w .* img + b;

end