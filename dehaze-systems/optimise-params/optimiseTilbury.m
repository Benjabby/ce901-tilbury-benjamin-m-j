function optimResults = optimiseTilbury(lossType, path, subSize, seed)
    
    gcp;

    if ~exist('lossType','var')
        lossType = "msssim";
    end
    
    if ~exist('path','var')
        path = "E:\optimise-dataset";
    end
    
%     validatePath = "E:\dataset";
    
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
    
    % These need to have the full structure, even if they are never used,
    % in order to appease the parfor loop in the objective function.
    % The par foor loop is annoying. This is the easiest way to deal with it
    subset.runningMeans = zeros(3,subSize,numSeqs);
    subset.IJSSIMs      = zeros(1,subSize,numSeqs);
    subset.IJMSEs       = zeros(1,subSize,numSeqs);
    
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
            
            if any(strcmp(["badCombo" "msssimprovement" "mseimprovement" "msssimmic"], lossType))
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
                    
                    if lossType == "badCombo" || lossType == "msssimmic"
                        subset.runningMeans(:,frameID+1,n) = squeeze(mean(rgb2ycbcr(J),[1 2]));
                    elseif lossType == "msssimprovement"
                        subset.IJSSIMs(:,frameID+1,n) = evaluateMSSSIM(I,J,true);
                    elseif lossType == "mseimprovement"
                       subset.IJMSEs(:,frameID+1,n)  = sqrt(mean((I-J).^2,'all'));
                    end
                end
            end
            
            subset.paths(n) = sequencePath;
            subset.starts(n) = start;
            subset.knowns(n) = knowns;
            subset.As(:,n) = trueA;
            
            n = n + 1;

        end
    end
    
    alpha = optimizableVariable('alpha',[1,100],'Transform','log');
    gamma = optimizableVariable('gamma',[0.7,1.3]);
    omega = optimizableVariable('omega',[0.7,1]);
    t0 = optimizableVariable('t0',[1e-5,0.2],'Transform','log');
    r = optimizableVariable('r',[3,45],'Type','integer');
    rA = optimizableVariable('rA',[3,60],'Type','integer');
    rM = optimizableVariable('rM',[3,60],'Type','integer');
    alphaM = optimizableVariable('alphaM',[1,100],'Transform','log');
    
    % IMPORTANT NOTE.
%      Due to how the bayesopt's algorithm (expected-improvement-per-second-plus) works, the results are not repeatable.
%      This is because the time each evaluation takes influences the selection of a new point, which is highly variable even on the same machine
%      If you run this you will end up with different (though likely very similar) values to that found in the paper.
%      Also, although bayesian optimisation is a global optimisation method, what I have found is that when running
%      the same type and loss multiple times, the optimiser converges on often very different parameters
%      (though with very very similar objective function values). 
%      I suspect this is due to the the aformentioned timing impact combined with the complexity of the
%      objective function space meaning that certain regions don't get considered at all even



    %% Optimise Local Type
    vars = [alpha,gamma,omega,t0,r,rA,rM];
    init = cell2table({20,1.25,0.95,0.1,15,30,15},'VariableNames',{'alpha','gamma','omega','t0','r','rA','rM'});
    
    lfg = @(a, b)(@(b)lossFunc('local',lossType,b,a));
    fun = lfg(subset);
    optimResults.local = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'Verbose',0,'PlotFcn',[]);

    %% Optimise Opening Type
    vars = [alpha,gamma,omega,t0,r,rA,rM,alphaM];
    init = cell2table({20,1.25,0.95,0.1,15,30,15,20},'VariableNames',{'alpha','gamma','omega','t0','r','rA','rM','alphaM'});
    
    lfg = @(a, b)(@(b)lossFunc('opening',lossType,b,a));
    fun = lfg(subset);
    optimResults.opening = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'Verbose',0,'PlotFcn',[]);

    %% Optimise Global Type
    vars = [alpha,gamma,omega,t0,r,rA];
    init = cell2table({20, 1.25, 0.95, 0.1, 15, 30},'VariableNames',{'alpha','gamma','omega','t0','r','rA'});

    lfg = @(a, b)(@(b)lossFunc('global',lossType,b,a));
    fun = lfg(subset);
    optimResults.global = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100);%,'Verbose',0,'PlotFcn',[]);

%     save("dehaze-systems\optimise-params\optims2.mat","optimResults");    
end



%     num = 7*3;
%     dehazers = HeDehazer; % HeDehazer for comparison
% 
%     i = 1;
%     for lossType = ["msssim","msssimT","msssimTA","msssimprovement","mae","rmse"]
%   
%         vars = [alpha,gamma,omega,r,rA,rM,t0];
%         init = cell2table({20,1.25,0.95,15,30,15,0.1},'VariableNames',{'alpha','gamma','omega','r','rA','rM','t0'});
%         
%         %% Optimise Default Mask Type
%         lfg = @(a, b)(@(b)lossFunc('default',lossType,b,a));
%         fun = lfg(subset);
%         result.(lossType).default = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100,'Verbose',0,'PlotFcn',[]);
%         
%         params = result.(lossType).default.bestPoint;
%         dh = TilburyDehazer('default',params.alpha,params.gamma,params.omega,params.r,params.rA,params.rM,params.t0);
%         dh.rename(append(lossType,"_default"));
%         dehazers = [dehazers dh];
%         fprintf("\n\n%d of %d done\n\n",i,num);
%         i = i+1;
%         
%         %% Optimise Saturation Mask Type
%         lfg = @(a, b)(@(b)lossFunc('sat',lossType,b,a));
%         fun = lfg(subset);
%         result.(lossType).sat = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100,'Verbose',0,'PlotFcn',[]);
%         
%         params = result.(lossType).sat.bestPoint;
%         dh = TilburyDehazer('sat',params.alpha,params.gamma,params.omega,params.r,params.rA,params.rM,params.t0);
%         dh.rename(append(lossType,"_sat"));
%         dehazers = [dehazers dh];
%         fprintf("\n\n%d of %d done\n\n",i,num);
%         i = i+1;
%         
%         %% Optimise Global Type
%         vars = [alpha,gamma,omega,r,rA,t0];
%         init = cell2table({20, 1.25, 0.95, 15, 30, 0.1},'VariableNames',{'alpha','gamma','omega','r','rA','t0'});
% 
%         lfg = @(a, b)(@(b)lossFunc('global',lossType,b,a));
%         fun = lfg(subset);
%         result.(lossType).global = bayesopt(fun,vars,'IsObjectiveDeterministic',true,'InitialX',init,'MaxObjectiveEvaluations',100,'Verbose',0,'PlotFcn',[]);
%         
%         params = result.(lossType).global.bestPoint;
%         dh = TilburyDehazer('global',params.alpha,params.gamma,params.omega,params.r,params.rA,1,params.t0);
%         dh.rename(append(lossType,"_global"));
%         dehazers = [dehazers dh];
%         fprintf("\n\n%d of %d done\n\n",i,num);
%         i = i+1;
%         
%         save("dehaze-systems\optimise-params\optims.mat","result");
%     end
%     
%     evalParams.Systems = dehazers;
%     evalParams.Sample = 20;
%     
%     dehazeResults = evaluateDehazers(validatePath,evalParams);
%     
%     table = groupResults(dehazeResults);
%     
%     save("dehaze-systems\optimise-params\tables.mat","table");
%     

function loss = lossFunc(method, lossType, vars, subset)

%     method = vars.method;
    a = vars.alpha;
    g = vars.gamma;
    o = vars.omega;
    t = vars.t0;
    r = vars.r;
    rA = vars.rA;
    if method=="global"
        rM = 1;
    else
        rM = vars.rM;
    end
    
    if method=="opening"
        alphaM = vars.alphaM;
    else 
        alphaM = 0;
    end
    
    lt = lossType;
    
    frameNum    = subset.frameNum;
    seqPaths    = subset.paths;
    seqKnowns   = subset.knowns;
    seqStarts   = subset.starts;
    seqAs       = subset.As;
    
    runningMeans = subset.runningMeans;
    IJSSIMs = subset.IJSSIMs;    
    IJMSEs = subset.IJMSEs;    
  
    seqERRS   = zeros(1,subset.numSeq);
    
%     for i = 1:length(seqPaths)
    parfor i = 1:length(seqPaths)
        
        dehazer = TilburyDehazer(method,a,g,o,t,r,rA,rM,alphaM);
        
        sequencePath = seqPaths(i);
        JPath = fullfile(sequencePath,'image_02','data');
        mPath = fullfile(sequencePath,'sky_mask','default_refined','img'); 
        IPath = fullfile(sequencePath,'haze','image');
        TPath = fullfile(sequencePath,'haze','transmission');
        
        A = reshape(seqAs(:,i),[1 1 3]);

        errs              = zeros(1,frameNum);
        
        if lt == "badCombo" || lt == "msssimmic"
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
            
            if lt == "msssimCol"
                errs(indexA) = evaluateMSSSIM(predImage,J,true)+colBlowout(predImage,J);
            elseif lt == "badCombo" % NOT RECOMMENDED. FADE is really bad as an objective function and the scales for the metrics all out of whack.
                trueT = imread(fullfile(TPath,framePath));
                trueT = im2double(trueT);
                predT = max(predT,t); 
                tError     = sqrt(mean((predT - trueT).^2,'all'));
                
                s = evaluateMSSSIM(predImage,J,true);
                
                fd = FADE(predImage);
                
                predRunningMeans(:,indexA) = squeeze(mean(rgb2ycbcr(predImage),[1 2]));
                
                errs(indexA) = tError*s*fd;
            elseif lt == "msssimmic"
                predRunningMeans(:,indexA) = squeeze(mean(rgb2ycbcr(predImage),[1 2]));
                errs(indexA) = evaluateMSSSIM(predImage,J,true);
            elseif lt == "msssimTA"
                trueT = imread(fullfile(TPath,framePath));
                trueT = im2double(trueT);
                predT = max(predT,t); 
                tError     = sqrt(mean((predT - trueT).^2,'all'));
                s = evaluateMSSSIM(predImage,J,true);
                aError     = sqrt(mean((predA - A).^2,'all'));
                errs(indexA) = s * tError * aError;
            elseif lt == "msssimT"
                trueT = imread(fullfile(TPath,framePath));
                trueT = im2double(trueT);
                predT = max(predT,t); 
                tError     = sqrt(mean((predT - trueT).^2,'all'));
                s = evaluateMSSSIM(predImage,J,true);
                errs(indexA) = s * tError;
            elseif lt == "fade" % NOT RECOMMENDED. FADE is really bad as an objective function
                errs(indexA) = FADE(predImage);
            elseif lt == "rmse" %|| lt == "rmseimprovement"
                errs(indexA) = sqrt(mean((predImage-J).^2,'all'));
            elseif lt == "mae"
                errs(indexA) = mean(abs(predImage-J),'all');
            else % MS-SSIM as default
                errs(indexA) = evaluateMSSSIM(predImage,J,true);
            end
            
        end
        
        
        if lt == "badCombo"
            seqERRS(i) = mean(errs) * evaluateMIC(runningMeans(:,:,i),predRunningMeans,true);
        elseif lt == "msssimmic"
            seqERRS(i) = mean(errs) + evaluateMIC(runningMeans(:,:,i),predRunningMeans,true);
        elseif lt == "msssimprovement"
            seqERRS(i) = mean(errs./IJSSIMs(:,:,i));
        elseif lt == "mseimprovement"
            seqERRS(i) = mean(errs./IJMSEs(:,:,i));
        else 
            seqERRS(i) = mean(errs);
        end
    end
    
    loss = mean(seqERRS);
end

function err = colBlowout(predImage,trueImage)
 
    hsvPred = rgb2hsv(predImage);
    hsvTrue = rgb2hsv(trueImage);
    
    v = (max(hsvPred(:,:,2)-hsvTrue(:,:,2),0).^2).*hsvPred(:,:,3);
    err = sqrt(mean(v,'all'));
end

function q = evaluateMSSSIM(predImage, trueImage, asLoss)
    q = mean(multissim(predImage,trueImage));
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