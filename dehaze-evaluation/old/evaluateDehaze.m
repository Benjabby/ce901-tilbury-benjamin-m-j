function results = evaluateDehaze(path, params)
    warning('on','all');
        
    %if ~exist('path','var') || isempty(path), path = fullfile(fileparts(mfilename('fullpath')),"..","haze-video-dataset","dataset"); end
    
    [results, systems, metrics, subset] = validateInputStruct(params);
     
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
                        runningMeans.Y = zeros(1,JSize);
                        runningMeans.Cb = zeros(1,JSize);
                        runningMeans.Cr = zeros(1,JSize);
                    end
%                     runningSSIM_J = zeros(1,JSize-1);
                else
                    if metrics.btcm, runningSSIM_I(frameID) = evaluateSSIM(prevI,I); end
                    if metrics.tcm,  runningSSIM_J(frameID) = evaluateSSIM(prevJ,J);  end
                end
                
                if metrics.mic
                    ycbcr = rgb2ycbcr(J);
                    runningMeans.Y(frameID+1) = mean(ycbcr(:,:,1),'all');
                    runningMeans.Cb(frameID+1) = mean(ycbcr(:,:,2),'all');
                    runningMeans.Cr(frameID+1) = mean(ycbcr(:,:,3),'all');
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
                            sequenceTotals.(name).runningMeans.Y = zeros(1,JSize);
                            sequenceTotals.(name).runningMeans.Cb = zeros(1,JSize);
                            sequenceTotals.(name).runningMeans.Cr = zeros(1,JSize);
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
                        ycbcr = rgb2ycbcr(predImage);
                        sequenceTotals.(name).runningMeans.Y(frameID+1-systems.(name).frameDelay) = mean(ycbcr(:,:,1),'all');
                        sequenceTotals.(name).runningMeans.Cb(frameID+1-systems.(name).frameDelay) = mean(ycbcr(:,:,2),'all');
                        sequenceTotals.(name).runningMeans.Cr(frameID+1-systems.(name).frameDelay) = mean(ycbcr(:,:,3),'all');
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
                        ycbcr = rgb2ycbcr(predImage);
                        sequenceTotals.(name).runningMeans.Y(end) = mean(ycbcr(:,:,1),'all');
                        sequenceTotals.(name).runningMeans.Cb(end) = mean(ycbcr(:,:,2),'all');
                        sequenceTotals.(name).runningMeans.Cr(end) = mean(ycbcr(:,:,3),'all');
                    end
                    
                    if frameID > systems.(name).frameDelay && (metrics.tcm || metrics.btcm)
                        sequenceTotals.(name).runningSSIM(end) = evaluateSSIM(sequenceTotals.(name).prevPred, predImage);
                    end
                end
%                 if systems.(name).frameDelay > 0
%                     for fd = 0:systems.(name).frameDelay
%                         [predImage, predT, predA, timeImage, timeA, systems.(name).state] = systems.(name).function([],knowns,systems.(name).state);     
%                     end
%                 end

%%
                
                mp = results.(seqStructStr).megapixels;
                
                if metrics.mppsImage
                    results.(seqStructStr).(name).mppsImage     = (1/(sequenceTotals.(name).timeImage/JSize))*mp;
                end
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
%%

function [results, systems, metrics, subset] = validateInputStruct(params)
    if isempty(params) || ~isstruct(params)
        params = struct();
    end
    
    if ~isfield(params,'Subset') || isempty(params.Subset)
       fid = fopen(fullfile(fileparts(mfilename('fullpath')),"..","haze-video-dataset","haze-subset.txt"));
       S = textscan(fid,"%s");
       subset = string(S{1});
       fclose(fid);
    else
        subset = params.Subset;
    end
    
    if ~isfield(params,'Systems') || isempty(params.Systems)
        systems = defaultSystems();
    elseif length(fieldnames(params.Systems))<2 || ~isfield(params.Systems,'prepped') || ~params.Systems.prepped
        isHandles = true;
        names = fieldnames(params.Systems);
        for s = 1:length(names)
            name = names{s};
           isHandles = isHandles & isa(params.Systems.(name),'function_handle');
        end
        
        if isHandles
            warning("'systems' has not been prepared properly, attempting to prepare now")
            systems = prepareSystems(params.Systems);
        else
            error("Incorrectly formatted 'Systems' structure")
        end
    else
        systems = params.Systems;
    end
    
    allMetrics = ["mppsImage" "mppsTotal" "mppsA" "AError" "psnr" "ssim" "fade" "colDiff" "mic" "tcm" "btcm" "TError", "brisque", "piqe"];
    useMetrics = ["mppsImage" "mppsTotal" "mppsA" "AError" "psnr" "ssim" "fade" "colDiff" "mic" "TError", "brisque", "piqe"];
    
    if ~isfield(params,'Metrics') || isempty(params.Metrics)
        textrics = useMetrics;
    else
        removals = {};
        for i = 1:length(params.Metrics)
            metric = params.Metrics(i);
            if ~any(strcmp(allMetrics,metric))
                warning("Unknown metric '%s' found in 'metrics'. Removing...",metric);
                removals{end+1} = i;
            end
        end
        for i = 1:length(removals)
            params.Metrics(removals{i}) = [];
        end
        textrics = params.Metrics;
    end
    
    metrics = struct();
    for metric = allMetrics
        metrics.(metric) = any(strcmp(textrics, metric));
    end
    
    if ~isfield(params,'Results') || isempty(params.Results)
        results = struct();
        results.names = string(fieldnames(systems))';
        results.metrics = textrics;
        results.colCode = struct();
        results.colCode.index = 1;
    else
        results.names = union(params.Results.names,string(fieldnames(systems))');
        results.metrics = union(params.Results.metrics,textrics);
    end
    
    results.names(strcmp(results.names,"prepped")) = [];
    
    fullCols = [
        31 119 180;
        255 127 14;
        44 160 44;
        214 39 40;
        148 103 189;
        140 86 75;
        227 119 194;
        127 127 127;
        188 189 34;
        23 190 207] / 255;
    
    for name = results.names
        if ~isfield(results.colCode,name)
            results.colCode.(name) = fullCols(results.colCode.index,:);
            results.colCode.index = results.colCode.index + 1;
        end
    end
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
    corrY = corrcoef(predMeans.Y, trueMeans.Y);
    corrY = corrY(1,2);
    corrCb = corrcoef(predMeans.Cb, trueMeans.Cb);
    corrCb = corrCb(1,2);
    corrCr = corrcoef(predMeans.Cr, trueMeans.Cr);
    corrCr = corrCr(1,2);
    
    q = (corrY*2+corrCb+corrCr)/4;
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