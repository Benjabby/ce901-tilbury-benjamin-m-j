function results = evaluateDehaze(path, systems, metrics, subset)
    warning('on','all');
        
    if ~exist('path','var') || isempty(path), path = fullfile(fileparts(mfilename('fullpath')),"..","haze-video-dataset","dataset"); end
    
    if ~exist('subset','var') || isempty(subset)
       fid = fopen(fullfile(fileparts(mfilename('fullpath')),"..","haze-video-dataset","haze-subset.txt"));
       S = textscan(fid,"%s");
       subset = string(S{1});
       fclose(fid);
       clearvars fid S;
    end
    

    if ~exist('systems','var') || isempty(systems)
        systems = defaultSystems();
    elseif length(fieldnames(systems))<2 || ~isfield(systems,'prepped') || ~systems.prepped
        isHandles = true;
        names = fieldnames(systems);
        for s = 1:length(names)
            name = names{s};
           isHandles = isHandles & isa(systems.(name),'function_handle');
        end
        
        if isHandles
            warning("'systems' has not been prepared properly, attempting to prepare now")
            systems = prepareSystems(systems);
        else
            error("Incorrectly formatted 'systems' structure")
        end
    end
    
    results.names = string(fieldnames(systems))';
    results.names(strcmp(results.names,"prepped")) = [];
    
    allMetrics = ["mppsImage" "mppsTotal" "mppsA" "AError" "psnr" "ssim" "fade" "colDiff" "mic" "tcm" "btcm" "TError", "brisque", "piqe"];
    % TODO Properly include only required metrics when calculating
    if ~exist('metrics','var') || isempty(metrics)
        results.metrics = allMetrics;
    else
        removals = {};
        for i = 1:length(metrics)
            metric = metrics(i);
            if ~any(strcmp(allMetrics,metric))
                warning("Unknown metric '%s' found in 'metrics'. Removing...",metric);
                removals{end+1} = i;
            end
        end
        for i = 1:length(removals)
            metrics(removals{i}) = [];
        end
        results.metrics = metrics;
    end
    
    tempMetrics = results.metrics;
    metrics = struct();
    for metric = allMetrics
        metrics.(metric) = any(strcmp(tempMetrics, metric));
    end
     
    dateFolders = dir(path);
    dateFolders = dateFolders([dateFolders.isdir]);
    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        calib = loadCalibrationCamToCam(fullfile(datePath,'calib_cam_to_cam.txt'));
        if isempty(calib)
            fprintf("Calibration file not found for date %s. Skipping\n",dateName)
            continue;
        end
        knowns.K = calib.P_rect{3};
        
        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        for j = 3:length(sequenceFolders)
            sequenceName = sequenceFolders(j).name;
            sequencePath = fullfile(datePath, sequenceName);
            
            if subset~="all"
                if ~any(strcmp(subset,sequenceName)), continue;  end
            end
            
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
            
            JPath = fullfile(sequencePath,'image_02','data');
            JSize = length(dir(JPath))-2;
            
            IPath = fullfile(sequencePath,'haze','image');
            
            if metrics.TError, TPath = fullfile(sequencePath,'haze','transmission'); end
            
            
            seqStructStr = ['S' sequenceName(1:end-5)];
            
            results.(seqStructStr).visibility = hazeProps.visibility; % Will be useful for comparing how different dehazing methods do on different fog densities
            results.(seqStructStr).frames = JSize;
            
            % Doing it frameX->system1,system2,... rather than systemX->frame1,frame2,...
            % to save having to load the images multiple times.
            % the only downside is needing to store state and 
            % results for each system which isn't a big deal
            % it also means that for delayed methods the previous J's and T's
            % need to be kept for as long as the delay, and that saving video is more of a hassle
            % To put it simply, doing it this way is slightly faster at the expense of requiring more memory.
            for frameID = 0:JSize-1
                framePath = num2str(frameID,'%010.f')+".png";

                J = imread(fullfile(JPath,framePath));
                J = im2double(J);

                I = imread(fullfile(IPath,framePath));
                I = im2double(I);
                                
                if frameID==0
                    results.(seqStructStr).megapixels = (size(J,1)*size(J,2))/(10^6);
                    
                    if metrics.btcm, runningSSIM_I = zeros(1,JSize-1); end
                    if metrics.tcm,  runningSSIM_J = zeros(1,JSize-1); end
                    
                    if metrics.mic
                        runningMeans.H = zeros(1,JSize);
                        runningMeans.S = zeros(1,JSize);
                        runningMeans.V = zeros(1,JSize);
                    end
%                     runningSSIM_J = zeros(1,JSize-1);
                else
                    if metrics.btcm, runningSSIM_I(frameID) = evaluateSSIM(prevI,I); end
                    if metrics.tcm,  runningSSIM_J(frameID) = evaluateSSIM(prevJ,J);  end
                end
                
                if metrics.mic
                    [hue, sat, val] = rgb2hsv(J);
                    runningMeans.H(frameID+1) = mean(hue,'all');
                    runningMeans.S(frameID+1) = mean(sat,'all');
                    runningMeans.V(frameID+1) = mean(val,'all');
                end
                
                if metrics.TError
                    trueT = imread(fullfile(TPath,framePath));
                    trueT = im2double(trueT);
                end
                
                textwaitbar(frameID,JSize,sprintf("Evaluating on %s",seqStructStr))
                
                for name = string(fieldnames(systems))'
                    if name == "prepped", continue; end
                    
                    if frameID==0
                        
                        if metrics.mppsTotal, sequenceTotals.(name).timeTotal    = 0; end
                        
                        if systems.(name).predictsA
                            if metrics.mppsA, sequenceTotals.(name).timeA    = 0; end
                            if metrics.AError, sequenceTotals.(name).AError   = 0; end
                        end
                            
                        if metrics.mppsImage, sequenceTotals.(name).timeImage    = 0; end
                        
                        if metrics.psnr, sequenceTotals.(name).psnr         = 0; end
                        if metrics.ssim, sequenceTotals.(name).ssim         = 0; end
                        if metrics.fade, sequenceTotals.(name).fade         = 0; end
                        
                        if metrics.piqe, sequenceTotals.(name).piqe         = 0; end
                        if metrics.brisque, sequenceTotals.(name).brisque      = 0; end
                        if metrics.colDiff, sequenceTotals.(name).colDiff      = 0; end
                        
                        if metrics.btcm || metrics.tcm, sequenceTotals.(name).runningSSIM  = zeros(1,JSize-1); end
                        
                        if metrics.mic
                            sequenceTotals.(name).runningMeans.H = zeros(1,JSize);
                            sequenceTotals.(name).runningMeans.S = zeros(1,JSize);
                            sequenceTotals.(name).runningMeans.V = zeros(1,JSize);
                        end
                    
                        if systems.(name).predictsT && metrics.TError
                           sequenceTotals.(name).TError    = 0;
                        end
                        
                        systems.(name).state = struct();
                    end
                    
                    totalTic = tic;
%                     try
                        [predImage, predT, predA, timeImage, timeA, systems.(name).state] = systems.(name).function(I,knowns,systems.(name).state);
%                     catch err
%                         sprintf("Error for %s:\n%s", name, err.message);
%                         continue
%                     end
                    timeTotal = toc(totalTic);
                    
                    if frameID < systems.(name).frameDelay
                        if metrics.mppsTotal, sequenceTotals.(name).timeTotal    = sequenceTotals.(name).timeTotal      + timeTotal; end
                        continue
                    end
                    
                    predImage = min(max(predImage,0),1);
                    predT = min(max(predT,0),1);
                    predA = min(max(predA,0),1);
                    
                    if systems.(name).frameDelay == 0
                        targetJ = J;
                        if metrics.TError, targetT = trueT; end
                    else
                        targetJ = prevJ;
                        if metrics.TError, targetT = prevT; end
                    end
                    
                    if metrics.mppsTotal, sequenceTotals.(name).timeTotal    = sequenceTotals.(name).timeTotal      + timeTotal; end
                    if metrics.mppsImage, sequenceTotals.(name).timeImage    = sequenceTotals.(name).timeImage      + timeImage; end
                    
                    if systems.(name).predictsA
                        if metrics.mppsA, sequenceTotals.(name).timeA    = sequenceTotals.(name).timeA          + timeA; end
                        if metrics.AError, sequenceTotals.(name).AError   = sequenceTotals.(name).AError         + evaluateAError(predA,trueA); end
                    end
                    
                    if metrics.psnr, sequenceTotals.(name).psnr         = sequenceTotals.(name).psnr           + evaluatePSNR(predImage,targetJ); end
                    if metrics.ssim, sequenceTotals.(name).ssim         = sequenceTotals.(name).ssim           + evaluateSSIM(predImage,targetJ); end
                    if metrics.fade, sequenceTotals.(name).fade         = sequenceTotals.(name).fade           + evaluateFADE(predImage); end
                    
                    if metrics.piqe, sequenceTotals.(name).piqe         = sequenceTotals.(name).piqe           + evaluatePIQE(predImage); end
                    if metrics.brisque, sequenceTotals.(name).brisque      = sequenceTotals.(name).brisque        + evaluateBRISQUE(predImage); end
                    if metrics.colDiff, sequenceTotals.(name).colDiff      = sequenceTotals.(name).colDiff        + evaluateColDiff(predImage,targetJ); end
                    
                    if systems.(name).predictsT && metrics.TError
                        sequenceTotals.(name).TError   = sequenceTotals.(name).TError         + evaluateTError(predT,targetT);
                    end
                    
                    if metrics.mic
                        [hue,sat,val] = rgb2hsv(predImage);
                        sequenceTotals.(name).runningMeans.H(frameID+1-systems.(name).frameDelay) = mean(hue,'all');
                        sequenceTotals.(name).runningMeans.S(frameID+1-systems.(name).frameDelay) = mean(sat,'all');
                        sequenceTotals.(name).runningMeans.V(frameID+1-systems.(name).frameDelay) = mean(val,'all');
                    end
                    
                    if frameID > systems.(name).frameDelay && (metrics.tcm || metrics.btcm)
                        sequenceTotals.(name).runningSSIM(frameID-systems.(name).frameDelay) = evaluateSSIM(sequenceTotals.(name).prevPred, predImage);
                    end
                    
                    sequenceTotals.(name).prevPred = predImage;
                end
                
                prevJ = J;
                prevI = I;
                if metrics.TError, prevT = trueT; end
            end
            textwaitbar(JSize,JSize,sprintf("Evaluating on %s",seqStructStr))
                
%             btcmDenom = evaluateBTCM(runningSSIM_J, runningSSIM_I)
            
            %assignin('base','totals',frameTotals);
            for name = string(fieldnames(systems))'
                if name == "prepped", continue; end
                
                %% Deal with frame delay by passing in empty images
                if systems.(name).frameDelay > 0
                    % Currently assuming delay of one
                    
                    [predImage, predT, predA, timeImage, timeA, systems.(name).state] = systems.(name).function([],knowns,systems.(name).state);

                    targetJ = J;
                    if metrics.TError, targetT = trueT; end
                    
                    if metrics.mppsTotal, sequenceTotals.(name).timeTotal    = sequenceTotals.(name).timeTotal      + timeTotal; end
                    if metrics.mppsImage, sequenceTotals.(name).timeImage    = sequenceTotals.(name).timeImage      + timeImage; end
                    
                    if systems.(name).predictsA
                        if metrics.mppsA, sequenceTotals.(name).timeA    = sequenceTotals.(name).timeA          + timeA; end
                        if metrics.AError, sequenceTotals.(name).AError   = sequenceTotals.(name).AError         + evaluateAError(predA,trueA); end
                    end
                    
                    if metrics.psnr, sequenceTotals.(name).psnr         = sequenceTotals.(name).psnr           + evaluatePSNR(predImage,targetJ); end
                    if metrics.ssim, sequenceTotals.(name).ssim         = sequenceTotals.(name).ssim           + evaluateSSIM(predImage,targetJ); end
                    if metrics.fade, sequenceTotals.(name).fade         = sequenceTotals.(name).fade           + evaluateFADE(predImage); end
                    
                    if metrics.piqe, sequenceTotals.(name).piqe         = sequenceTotals.(name).piqe           + evaluatePIQE(predImage); end
                    if metrics.brisque, sequenceTotals.(name).brisque      = sequenceTotals.(name).brisque        + evaluateBRISQUE(predImage); end
                    if metrics.colDiff, sequenceTotals.(name).colDiff      = sequenceTotals.(name).colDiff        + evaluateColDiff(predImage,targetJ); end
                    
                    if systems.(name).predictsT && metrics.TError
                        sequenceTotals.(name).TError   = sequenceTotals.(name).TError         + evaluateTError(predT,targetT);
                    end
                    
                    if metrics.mic
                        [hue,sat,val] = rgb2hsv(predImage);
                        sequenceTotals.(name).runningMeans.H(frameID+1-systems.(name).frameDelay) = mean(hue,'all');
                        sequenceTotals.(name).runningMeans.S(frameID+1-systems.(name).frameDelay) = mean(sat,'all');
                        sequenceTotals.(name).runningMeans.V(frameID+1-systems.(name).frameDelay) = mean(val,'all');
                    end
                    
                    if frameID > systems.(name).frameDelay && (metrics.tcm || metrics.btcm)
                        sequenceTotals.(name).runningSSIM(frameID-systems.(name).frameDelay) = evaluateSSIM(sequenceTotals.(name).prevPred, predImage);
                    end
                end
%                 if systems.(name).frameDelay > 0
%                     for fd = 0:systems.(name).frameDelay
%                         [predImage, predT, predA, timeImage, timeA, systems.(name).state] = systems.(name).function([],knowns,systems.(name).state);     
%                     end
%                 end

%%
                
                mp = results.(seqStructStr).megapixels;
                
                if metrics.mppsImage, results.(seqStructStr).(name).mppsImage     = (1/(sequenceTotals.(name).timeImage/JSize))*mp;end
                if metrics.mppsTotal, results.(seqStructStr).(name).mppsTotal     = (1/(sequenceTotals.(name).timeTotal/JSize))*mp;end
                
                if systems.(name).predictsA
                    if metrics.mppsA, results.(seqStructStr).(name).mppsA     = (1/(sequenceTotals.(name).timeA/JSize))*mp; end
                    if metrics.AError, results.(seqStructStr).(name).AError    = sequenceTotals.(name).AError/JSize; end
                else
                    if metrics.mppsA, results.(seqStructStr).(name).mppsA     = NaN; end
                    if metrics.AError, results.(seqStructStr).(name).AError    = NaN;end
                end
                
                if metrics.psnr, results.(seqStructStr).(name).psnr          =  sequenceTotals.(name).psnr / JSize; end
                if metrics.ssim, results.(seqStructStr).(name).ssim          =  sequenceTotals.(name).ssim / JSize;end
                if metrics.fade, results.(seqStructStr).(name).fade          =  sequenceTotals.(name).fade / JSize;end

                if metrics.mic, results.(seqStructStr).(name).mic            = evaluateMIC(sequenceTotals.(name).runningMeans,runningMeans);end
                if metrics.tcm, results.(seqStructStr).(name).tcm           = evaluateTCM(sequenceTotals.(name).runningSSIM, runningSSIM_J);end
                if metrics.btcm, results.(seqStructStr).(name).btcm          = evaluateTCM(sequenceTotals.(name).runningSSIM, runningSSIM_I);end
                
                if metrics.piqe, results.(seqStructStr).(name).piqe          = sequenceTotals.(name).piqe / JSize;end
                if metrics.brisque, results.(seqStructStr).(name).brisque       = sequenceTotals.(name).brisque / JSize;end
                if metrics.colDiff, results.(seqStructStr).(name).colDiff       =  sequenceTotals.(name).colDiff / JSize;end
                
                if systems.(name).predictsT
                    if metrics.TError, results.(seqStructStr).(name).TError    = sequenceTotals.(name).TError/JSize; end
                else
                    if metrics.TError, results.(seqStructStr).(name).TError    = NaN;end
                end
            end
        end
    end
%     assignin('base','ssimI',runningSSIM_I);

end

%% Evaluation Metrics
% A lot of these are superfluous wrappers

%% Real-time Performance
% - Atmospheric light FPS
% - Dehaze FPS
% - Total FPS

%% Visibility Improvement
function q = evaluatePSNR(pred, trueImage)
    if any(size(pred,[1 2])~=size(trueImage,[1 2]))
       wind = centerCropWindow2d(size(trueImage),size(pred));
       trueImage = imcrop(trueImage, wind);
    end
    q = psnr(pred, trueImage);
end

function q = evaluateSSIM(pred, trueImage)
    if any(size(pred,[1 2])~=size(trueImage,[1 2]))
       wind = centerCropWindow2d(size(trueImage),size(pred));
       trueImage = imcrop(trueImage, wind);
    end
    q = ssim(pred, trueImage);
end

function q = evaluateFADE(pred)
    q = real(FADE(pred));
end

%% Image Quality
function q = evaluatePIQE(pred)
    q = piqe(pred);
end

function q = evaluateBRISQUE(pred)
    q = brisque(pred);
end

function q = evaluateColDiff(pred, trueImage)
    if any(size(pred,[1 2])~=size(trueImage,[1 2]))
       wind = centerCropWindow2d(size(trueImage),size(pred));
       trueImage = imcrop(trueImage, wind);
    end
    q = mean(imcolordiff(pred,trueImage,"Standard","CIEDE2000"),'all');
end

%% Temporal Consistency
function q = evaluateMIC(predMeans, trueMeans)
    corrH = corrcoef(predMeans.H, trueMeans.H);
    corrH = corrH(1,2);
    corrS = corrcoef(predMeans.S, trueMeans.S);
    corrS = corrS(1,2);
    corrV = corrcoef(predMeans.V, trueMeans.V);
    corrV = corrV(1,2);
    
    q = (corrH+corrS+corrV)/3;
end

function q = evaluateTCM(ssimA, ssimB)
    ssimCorr = corrcoef(ssimA,ssimB);
    q = ssimCorr(1,2);
end


%% Transmission Error
function q = evaluateTError(predT, trueT)
    if any(size(predT,[1 2])~=size(trueT,[1 2]))
       wind = centerCropWindow2d(size(trueT),size(pred));
       trueT = imcrop(trueT, wind);
    end
    q = sqrt(mean((predT - trueT).^2,'all'));
end

%% Atmospheric Light Error
function q = evaluateAError(predA, trueA)
    predA = squeeze(predA);
    trueA = squeeze(trueA);
%     if size(predA,1) ~= 1
%         trueA = repmat(trueA, m, n);
    q = sqrt(mean((predA - trueA).^2,'all'));
end