function table = temp(path, subSize, seed)

    gcp;

    if ~exist('path','var')
        path = "E:\dataset";
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
            JSize = length(dir(JPath))-2;
            
            mPath = fullfile(sequencePath,'sky_mask','default_refined','img'); 
            
            % Don't include first and last frame.
            start = randi(JSize-subSize-1);
            
% % %             IPath = fullfile(sequencePath,'haze','image');
% % %             TPath = fullfile(sequencePath,'haze','transmission');
% % %             
% % %             framePath = num2str(0,'%010.f')+".png";
% % %             temp = imread(fullfile(JPath,framePath));
% % %             
% % %             imsize = size(temp);
% % %             subset.Js{n} = zeros([subSize imsize]);
% % %             subset.Is{n} = zeros([subSize imsize]);
% % %             subset.Ts{n} = zeros([subSize imsize]);
% % %             
            for frameID = 0:subSize-1
% % %             
                framePath = num2str(frameID+start,'%010.f')+".png";
% % %                 
                J = imread(fullfile(JPath,framePath));
                J = im2double(J);
                
                m = imread(fullfile(mPath,framePath));
                m = im2double(m);
                m = repmat(m,1,1,3);
                J = m.*A + (1-m).*J;
                
                subset.runningMeans(:,frameID+1,n) = squeeze(mean(rgb2ycbcr(J),[1 2]));
% % % 
% % %                 I = imread(fullfile(IPath,framePath));
% % %                 I = im2double(I);
% % %             
% % %                 trueT = imread(fullfile(TPath,framePath));
% % %                 trueT = im2double(trueT);
% % %                 
% % %                 subset.Js{n}(frameID+1, :, :, :) = J;
% % %                 subset.Is{n}(frameID+1, :, :, :) = I;
% % %                 subset.Ts{n}(frameID+1, :, :, :) = trueT;
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
    rA = optimizableVariable('rA',[3,60],'Type','integer');
    rM = optimizableVariable('rM',[3,60],'Type','integer');
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
    init = cell2table({33.712,1.2009,0.90867,26,53,51,0.1494},'VariableNames',{'alpha','lambda','omega','r','rA','rM','t0'});

    % IMPORTANT NOTE.
%      Due to how bayesopt's default algorithm (expected-improvement-per-second-plus) works, the results are not repeatable.
%      This is because the time each evaluation takes influences the selection of a new point, which is highly variable even on the same machine
%      If you run this you will end up with different (though likely very similar) values to that found in the paper.

%     %% Optimise with SSIM only
%     lfg = @(a, b)(@(b)lossFunc('default','ssim',b,a));
%     fun = lfg(subset);
%     ssimOptim = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'UseParallel',true);
    
    %% Optimise with everything
    lfg = @(a, b)(@(b)lossFunc('default','multall',b,a));
    fun = lfg(subset);
    allOptim = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'UseParallel',true);
    
    %% Optimise with FADE
    lfg = @(a, b)(@(b)lossFunc('default','fade',b,a));
    fun = lfg(subset);
    fadeOptim = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'UseParallel',true);
    
    %% Optimise with SSIM and TError
%     lfg = @(a, b)(@(b)lossFunc('default','ssim_t',b,a));
%     fun = lfg(subset);
%     optim4 = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'UseParallel',true);
%     
%     %% Optimise with MAE
%     lfg = @(a, b)(@(b)lossFunc('default','mae',b,a));
%     fun = lfg(subset);
%     optim5 = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'UseParallel',true);
    
%     d1 = TilburyDehazer('default', 45.792,1.2559 , 0.9151, 21, 60, 60, 0.0055583);
%     d1.rename("RMSE");
%     
%     d2 = TilburyDehazer('default',9.97097095180161,1.37912519131371,0.956081121476979,29,42,45,0.192701658635203);
%     d2.rename("SSIM_TError_NoClip");
    
    ssimDehazer = TilburyDehazer('default',33.712,1.2009,0.90867,26,53,51,0.1494);
    ssimDehazer.rename("SSIM_2");
    
    fadeParams = fadeOptim.bestPoint;
    fadeDehazer = TilburyDehazer('default',fadeParams.alpha,fadeParams.lambda,fadeParams.omega,fadeParams.r,fadeParams.rA,fadeParams.rM,fadeParams.t0);
    fadeDehazer.rename("FADE");
    
    allParams = allOptim.bestPoint;
    allDehazer = TilburyDehazer('default',allParams.alpha,allParams.lambda,allParams.omega,allParams.r,allParams.rA,allParams.rM,allParams.t0);
    allDehazer.rename("All");
    
%     params4 = optim4.bestPoint;
%     d4 = TilburyDehazer('default',params4.alpha,params4.lambda,params4.omega,params4.r,params4.rA,params4.rM,params4.t0);
%     d4.rename("SSIM_TError_Clip");
%     
%     params5 = optim5.bestPoint;
%     d5 = TilburyDehazer('default',params5.alpha,params5.lambda,params5.omega,params5.r,params5.rA,params5.rM,params5.t0);
%     d5.rename("MAE");
    
%     he = HeDehazer;
    
    params.Systems = [ssimDehazer fadeDehazer allDehazer];
    params.Sample = 20;
    
    result = evaluateDehazers(path,params);
    
    table = groupResults(result);
    
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
        if lt == "multall"
            predRunningMeans   = zeros(3,frameNum);
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
            
            [predImage, predT, ~, ~, ~] = dehazer.dehazeFrame(I);
            
            if lt == "multall"
                trueT = imread(fullfile(TPath,framePath));
                trueT = im2double(trueT);
                predT = max(predT,t); 
                tError     = sqrt(mean((predT - trueT).^2,'all'));
                
                ssimr = ssim(predImage(:,:,1),J(:,:,1));
                ssimg = ssim(predImage(:,:,2),J(:,:,2));
                ssimb = ssim(predImage(:,:,3),J(:,:,3));
                s = 1.0 - (ssimr+ssimg+ssimb)/3;
                
                fd = FADE(predImage);
                
                predRunningMeans(:,indexA) = squeeze(mean(rgb2ycbcr(predImage),[1 2]));
                
                errs(indexA) = tError*s*fd;
            elseif lt == "ssim_t"
                trueT = imread(fullfile(TPath,framePath));
                trueT = im2double(trueT);
                predT = max(predT,t); 
                tError     = sqrt(mean((predT - trueT).^2,'all'));
                ssimr = ssim(predImage(:,:,1),J(:,:,1));
                ssimg = ssim(predImage(:,:,2),J(:,:,2));
                ssimb = ssim(predImage(:,:,3),J(:,:,3));
                s = (ssimr+ssimg+ssimb)/3;
                errs(indexA) = 1.0-s + tError;
%                 errs(indexA) = errs(indexA)./0.248143367922327;
            elseif lt == "ssim"
                ssimr = ssim(predImage(:,:,1),J(:,:,1));
                ssimg = ssim(predImage(:,:,2),J(:,:,2));
                ssimb = ssim(predImage(:,:,3),J(:,:,3));
                s = (ssimr+ssimg+ssimb)/3;
                errs(indexA) = 1.0-s;
%                 errs(indexA) = errs(indexA)./0.106019897218699;
            elseif lt == "fade"
                errs(indexA) = FADE(predImage);
            else
                errs(indexA) = mean(abs(predImage-J),'all');
%                 errs(indexA) = errs(indexA)./0.063010440792390;
            end
            
%             errs(indexA) = sqrt(mean((predImage - J).^2,'all'));
        end
        
        seqERRS(i) = mean(errs)
        if lt == "multall"
            seqERRS(i) = seqERRS(i) * (1.0-evaluateMIC(runningMeans(:,:,i),predRunningMeans));
        end
    end
    
    loss = mean(seqERRS);
end

%% Temporal Consistency
function q = evaluateMIC(predMeans, trueMeans)
    corrY = corrcoef(predMeans(1,:), trueMeans(1,:));
    corrY = corrY(1,2);
    corrCb = corrcoef(predMeans(2,:), trueMeans(2,:));
    corrCb = corrCb(1,2);
    corrCr = corrcoef(predMeans(3,:), trueMeans(3,:));
    corrCr = corrCr(1,2);
    
    q = (corrY*2+corrCb+corrCr)/4;
end
