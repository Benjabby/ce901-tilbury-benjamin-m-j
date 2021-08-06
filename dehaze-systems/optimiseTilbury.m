function result = optimiseTilbury(path, subSize, seed)

    gcp;

    if ~exist('path','var')
        path = "E:\optimise-dataset";
    end
    
    dateFolders = dir(path);
    dateFolders = dateFolders([dateFolders.isdir]);
    
    if ~exist('subSize','var')
        subSize = 20;
    end
    
    if ~exist('seed','var')
        rng(1);
    else
        rng(seed)
    end
    
    numSeqs = 0;
    
    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        if ~isfile(fullfile(datePath,'calib_cam_to_cam.txt')), continue; end

        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        numSeqs = numSeqs + length(sequenceFolders) - 2;

    end
    
    
    subset.frameNum = subSize;
    subset.numSeq   = numSeqs;
    subset.paths    = strings(1,numSeqs);
    subset.starts   = zeros(1,numSeqs);
    subset.As       = zeros(3,numSeqs);
    subset.knowns(1,numSeqs)     = struct('K',{zeros(3,4)});
    
    subset.runningMeans = zeros(3,subSize,numSeqs);
    subset.IJSSIMs      = zeros(1,subSize,numSeqs);
    
% % %     subset.Js       = cell(1,numSeqs);
% % %     subset.Is       = cell(1,numSeqs);
% % %     subset.Ts       = cell(1,numSeqs);
    
    n = 1;
    
    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        if ~isfile(fullfile(datePath,'calib_cam_to_cam.txt')), continue; end

        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        calib = loadCalibrationCamToCam(fullfile(datePath,'calib_cam_to_cam.txt'));
        K = calib.P_rect{3};
        knowns.K = K;

        for j = 3:length(sequenceFolders)
            sequenceName = sequenceFolders(j).name;
            sequencePath = fullfile(datePath, sequenceName);
            
            fid = fopen(fullfile(sequencePath,'haze','props.json'));
            if fid<0
                fprintf("props.json not found for sequence %s. Skipping\n",sequenceName)
                continue;
            end
            rawJSON = fread(fid);
            str = char(rawJSON');
            fclose(fid);
            hazeProps = jsondecode(str);

            trueA = hazeProps.A;
            A = reshape(trueA, [1 1 3]);
            
            JPath = fullfile(sequencePath,'image_02','data');
            IPath = fullfile(sequencePath,'haze','image');
            mPath = fullfile(sequencePath,'sky_mask','default_refined','img'); 
            
            JSize = length(dir(JPath))-2;
            
            % Don't include first and last frame.
            start = randi(JSize-subSize-1);
            
            for frameID = 0:subSize-1
                framePath = num2str(frameID+start,'%010.f')+".png";
                J = imread(fullfile(JPath,framePath));
                J = im2double(J);
                
                m = imread(fullfile(mPath,framePath));
                m = im2double(m);
                m = repmat(m,1,1,3);
                J = m.*A + (1-m).*J;
                
                I = imread(fullfile(IPath,framePath));
                I = im2double(I);
                
                subset.runningMeans(:,frameID+1,n) = squeeze(mean(rgb2ycbcr(J),[1 2]));
                subset.IJSSIMs(:,frameID+1,n) = evaluateSSIM(I,J,true);
            end
            
            subset.paths(n) = sequencePath;
            subset.starts(n) = start;
            subset.knowns(n) = knowns;
            subset.As(:,n) = trueA;
            
            n = n + 1;

        end
    end
    
    alpha = optimizableVariable('alpha',[1,100],'Transform','log');
    lambda = optimizableVariable('lambda',[0.8,1.4]);
    omega = optimizableVariable('omega',[0.7,1]);
    r = optimizableVariable('r',[3,45],'Type','integer');
    rA = optimizableVariable('rA',[3,120],'Type','integer');
    rM = optimizableVariable('rM',[3,120],'Type','integer');
    t0 = optimizableVariable('t0',[1e-5,0.2],'Transform','log');

%     minvd = optimizableVariable('minvd',[1,100]);
%     method = optimizableVariable('method',{'default','sat','global'},'Type','categorical');

    vars = [alpha,lambda,omega,r,rA,rM,t0];
    
    %% A previously good result:
    % omega = 0.956081121476979
    % alpha = 9.97097095180161
    % lambda = 1.37912519131371
    % r = 29
    % rA = 42
    % t0 = 0.192701658635203
    
    %% Another good result:

    init = cell2table({20,1.25,0.95,15,30,15,0.01},'VariableNames',{'alpha','lambda','omega','r','rA','rM','t0'});
    
    result = cell(1,3);
    
    % IMPORTANT NOTE.
%      Due to how bayesopt's default algorithm (expected-improvement-per-second-plus) works, the results are not repeatable.
%      This is because the time each evaluation takes influences the selection of a new point, which is highly variable even on the same machine
%      If you run this you will end up with different (though likely very similar) values to that found in the paper.

    %% Optimise using the SSIMprovement loss,
%     lossType = "ssim_t_a";
    lossType = "ssim_t";
    
    %% Optimise Default Mask Type
    lfg = @(a, b)(@(b)lossFunc('default',lossType,b,a));
    fun = lfg(subset);
    result{1} = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'UseParallel',true);
    
    %% Optimise Saturation Mask Type
%     lfg = @(a, b)(@(b)lossFunc('sat',lossType,b,a));
%     fun = lfg(subset);
%     result{2} = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'UseParallel',true);
%     
%     %% Optimise Global Alpha Type
%     vars = [alpha,omega,r,rA,t0];
%     init = cell2table({20, 0.95, 15, 30, 0.01},'VariableNames',{'alpha','omega','r','rA','t0'});
%     
%     lfg = @(a, b)(@(b)lossFunc('global',lossType,b,a));
%     fun = lfg(subset);
%     result{3} = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'UseParallel',true);
end

function loss = lossFunc(method, lossType, vars, subset)

%     method = vars.method;
    a = vars.alpha;
    
    if method=="global"
        l = 0;
        rM = 0;
    else
        l = vars.lambda;
        rM = vars.rM;
    end
    
    o = vars.omega;
    r = vars.r;
    rA = vars.rA;
    t = vars.t0;
    
    lt = lossType;
    
    frameNum    = subset.frameNum;
    seqPaths    = subset.paths;
    seqKnowns   = subset.knowns;
    seqStarts   = subset.starts;
    seqAs       = subset.As;
    
    runningMeans = subset.runningMeans;
    IJSSIMs = subset.IJSSIMs;
    
    seqERRS   = zeros(1,subset.numSeq);
    
    parfor i = 1:length(seqPaths)
        
        dehazer = TilburyDehazer('default',a,l,o,r,rA,rM,t);
        
        sequencePath = seqPaths(i);
        JPath = fullfile(sequencePath,'image_02','data');
        mPath = fullfile(sequencePath,'sky_mask','default_refined','img'); 
        IPath = fullfile(sequencePath,'haze','image');
        TPath = fullfile(sequencePath,'haze','transmission');
        
        A = reshape(seqAs(:,i),[1 1 3]);

        errs              = zeros(1,frameNum);
        if lt == "multiple"
            predRunningMeans   = zeros(3,frameNum);
        else
            predRunningMeans   = []; % To appease the parallel loop.
        end
        dehazer.newSequence(seqKnowns(i));
        
        for frameID = 0:frameNum-1
            
            framePath = num2str(frameID+seqStarts(i),'%010.f')+".png";
            
            indexA = frameID + 1;

            J = imread(fullfile(JPath,framePath));
            J = im2double(J);
            
            m = imread(fullfile(mPath,framePath));
            m = im2double(m);
            m = repmat(m,1,1,3);
            J = m.*A + (1-m).*J;

            I = imread(fullfile(IPath,framePath));
            I = im2double(I);            
            
            [predImage, predT, predA, ~, ~] = dehazer.dehazeFrame(I);
            
            if lt == "multiple" % Does not work well. FADE is really bad as an objective function
                trueT = imread(fullfile(TPath,framePath));
                trueT = im2double(trueT);
                predT = max(predT,t); 
                tError     = sqrt(mean((predT - trueT).^2,'all'));
                
                s = evaluateSSIM(predImage,J,true);
                
                fd = FADE(predImage);
                
                predRunningMeans(:,indexA) = squeeze(mean(rgb2ycbcr(predImage),[1 2]));
                
                errs(indexA) = tError*s*fd;
            elseif lt == "ssim_t_a"
                trueT = imread(fullfile(TPath,framePath));
                trueT = im2double(trueT);
                predT = max(predT,t); 
                tError     = sqrt(mean((predT - trueT).^2,'all'));
                s = evaluateSSIM(predImage,J,true);
                aError     = sqrt(mean((predA - A).^2,'all'));
                errs(indexA) = s * tError * aError;
            elseif lt == "ssim_t"
                trueT = imread(fullfile(TPath,framePath));
                trueT = im2double(trueT);
                predT = max(predT,t); 
                tError     = sqrt(mean((predT - trueT).^2,'all'));
                s = evaluateSSIM(predImage,J,true);
                errs(indexA) = s * tError;
            elseif lt == "fade" % Does not work well. FADE is really bad as an objective function
                errs(indexA) = FADE(predImage);
            elseif lt == "mae"
                errs(indexA) = mean(abs(predImage-J),'all');
            else % SSIM as default
                errs(indexA) = evaluateSSIM(predImage,J,true);
            end
            
        end
        
        
        if lt == "multiple"
            seqERRS(i) = mean(errs) * evaluateMIC(runningMeans(:,:,i),predRunningMeans,true);
        elseif lt == "ssimprovement"
            seqERRS(i) = mean(errs./IJSSIMs(:,:,i));
        else 
            seqERRS(i) = mean(errs);
        end
    end
    
    loss = mean(seqERRS);
end

function q = evaluateSSIM(predImage, trueImage, asLoss)

    ssimr = ssim(predImage(:,:,1),trueImage(:,:,1));
    ssimg = ssim(predImage(:,:,2),trueImage(:,:,2));
    ssimb = ssim(predImage(:,:,3),trueImage(:,:,3));
    q = (ssimr+ssimg+ssimb)/3;
    
    if nargin>2 && asLoss
        q = 1.0-q;
    end
end

function q = evaluateMIC(predMeans, trueMeans, asLoss)
    corrY = corrcoef(predMeans(1,:), trueMeans(1,:));
    corrY = corrY(1,2);
    corrCb = corrcoef(predMeans(2,:), trueMeans(2,:));
    corrCb = corrCb(1,2);
    corrCr = corrcoef(predMeans(3,:), trueMeans(3,:));
    corrCr = corrCr(1,2);
    
    q = (corrY*2+corrCb+corrCr)/4;
    
    if nargin>2 && asLoss
        q = 1.0-q;
    end
end