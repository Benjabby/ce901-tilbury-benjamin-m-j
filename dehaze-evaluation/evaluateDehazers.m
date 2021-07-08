function results = evaluateDehazers(path, params)
    warning('on','all');
 
    if ~exist('params','var'), params = []; end
    [results, systems, metrics, subset, autosave] = validateInputStruct(params);
     
    dateFolders = dir(path);
    dateFolders = dateFolders([dateFolders.isdir]);
    
    remainingFrames = 0;
    remainingSequences = 0;
    
    %% Preliminary look to find total number of frames
    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        if ~isfile(fullfile(datePath,'calib_cam_to_cam.txt')), continue; end
        
        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        for j = 3:length(sequenceFolders)
            sequenceName = sequenceFolders(j).name;
            sequencePath = fullfile(datePath, sequenceName);
            
            if ~any(strcmp(subset,sequenceName)), continue;  end
            
            remainingSequences = remainingSequences + 1;
            
            JPath = fullfile(sequencePath,'image_02','data');
            JSize = length(dir(JPath))-2;
            
            remainingFrames = remainingFrames+JSize;
        end
    end
    
    frameTimeEst = 0;
    seqOverheadEst = 0;
    timeEst = 0;
    
    frameTimeTotal = 0;
    seqOverheadTotal = 0;
    
    completedFrames = 0;
    completedSequences = 0;
    
    fprintf("Evaluating on %d sequences totaling %d frames\n", remainingSequences, remainingFrames);
    
    %% Loop Drives
    for i = 3:length(dateFolders)
        seqTic = tic;
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
        
        seqOverheadTotal = seqOverheadTotal + toc(seqTic);
        
        %% Loop Sequence
        for j = 3:length(sequenceFolders)
            
            seqTic = tic;
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
            
            results.(seqStructStr).props.visibility = hazeProps.visibility; % Will be useful for comparing how different dehazing methods do on different fog densities
            results.(seqStructStr).props.frames = JSize;
            
            seqOverheadTotal = seqOverheadTotal + toc(seqTic);
            
            % Doing it frameX->system1,system2,... rather than systemX->frame1,frame2,...
            % to save having to load the images multiple times.
            % the only downside is needing to store state and 
            % results for each system which isn't a big deal
            % it also means that for delayed methods the previous J's and T's
            % need to be kept for as long as the delay, and that saving video is more of a hassle
            % To put it simply, doing it this way is slightly faster at the expense of requiring more memory.
            %% Loop Frames
            for frameID = 0:JSize-1
                frameTic = tic;
                framePath = num2str(frameID,'%010.f')+".png";

                J = imread(fullfile(JPath,framePath));
                J = im2double(J);

                I = imread(fullfile(IPath,framePath));
                I = im2double(I);
                                
                if frameID==0
                    results.(seqStructStr).props.megapixels = (size(J,1)*size(J,2))/(10^6);
                    results.(seqStructStr).props.height = size(J,1);
                    
                    if metrics.mic
                        runningMeans.Y = zeros(1,JSize);
                        runningMeans.Cb = zeros(1,JSize);
                        runningMeans.Cr = zeros(1,JSize);
                    end
%                     runningSSIM_J = zeros(1,JSize-1);
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
                
                %% Loop Dehazers
                for dehazer = systems
                    name = erase(class(dehazer),"Dehazer");
                    
                    if frameID==0
                        
                        if metrics.mppsTotal,           results.(seqStructStr).(name).timeTotal.all         = zeros(1,JSize+dehazer.FrameDelay); end
                        
                        if dehazer.PredictsA
                            if metrics.mppsA,           results.(seqStructStr).(name).timeA.all             = zeros(1,JSize); end
                            if metrics.AError,          results.(seqStructStr).(name).AError.all            = zeros(1,JSize); end
                        end
                            
                        if metrics.mppsImage,           results.(seqStructStr).(name).timeImage.all         = zeros(1,JSize); end
                        
                        if metrics.psnr,                results.(seqStructStr).(name).psnr.all              = zeros(1,JSize); end
                        if metrics.ssim,                results.(seqStructStr).(name).ssim.all              = zeros(1,JSize); end
                        if metrics.fade,                results.(seqStructStr).(name).fade.all              = zeros(1,JSize); end
                        
                        if metrics.piqe,                results.(seqStructStr).(name).piqe.all              = zeros(1,JSize); end
                        if metrics.brisque,             results.(seqStructStr).(name).brisque.all           = zeros(1,JSize); end
                        if metrics.colDiff,             results.(seqStructStr).(name).colDiff.all           = zeros(1,JSize); end
                        
                        
                        if metrics.mic
                                                        sequenceTotals.(name).runningMeans.Y = zeros(1,JSize);
                                                        sequenceTotals.(name).runningMeans.Cb = zeros(1,JSize);
                                                        sequenceTotals.(name).runningMeans.Cr = zeros(1,JSize);
                        end
                    
                        if dehazer.PredictsT && metrics.TError
                                                        results.(seqStructStr).(name).TError.all            = zeros(1,JSize);
                        end
                        
                        dehazer.newSequence(knowns);
                    end
                    
                    totalTic = tic;
%                     try
                        [predImage, predT, predA, timeImage, timeA] = dehazer.dehazeFrame(I);
%                     catch err
%                         sprintf("Error for %s:\n%s", name, err.message);
%                         continue
%                     end
                    timeTotal = toc(totalTic);
                    
                    indexA = frameID+1;
                    indexB = frameID+1-dehazer.FrameDelay;
                    
                    if frameID < dehazer.FrameDelay
                        if metrics.mppsTotal,           results.(seqStructStr).(name).timeTotal.all(indexA)    = timeTotal; end
                        continue
                    end
                    
                    
%                     predImage = min(max(predImage,0),1);
%                     predT = min(max(predT,0),1);
%                     predA = min(max(predA,0),1);
                    
                    if dehazer.FrameDelay == 0
                        targetJ = J;
                        if metrics.TError, targetT = trueT; end
                    else
                        targetJ = prevJ;
                        if metrics.TError, targetT = prevT; end
                    end
                    
                    if metrics.mppsTotal,               results.(seqStructStr).(name).timeTotal.all(indexA)     = timeTotal; end
                    if metrics.mppsImage,               results.(seqStructStr).(name).timeImage.all(indexB)     = timeImage; end
                    
                    if dehazer.PredictsA
                        if metrics.mppsA,               results.(seqStructStr).(name).timeA.all(indexB)         = timeA; end
                        if metrics.AError,              results.(seqStructStr).(name).AError.all(indexB)        = evaluateAError(predA,trueA); end
                    end
                    
                    if metrics.psnr,                    results.(seqStructStr).(name).psnr.all(indexB)          = evaluatePSNR(predImage,targetJ); end
                    if metrics.ssim,                    results.(seqStructStr).(name).ssim.all(indexB)          = evaluateSSIM(predImage,targetJ); end
                    if metrics.fade,                    results.(seqStructStr).(name).fade.all(indexB)          = evaluateFADE(predImage); end
                    
                    if metrics.piqe,                    results.(seqStructStr).(name).piqe.all(indexB)          = evaluatePIQE(predImage); end
                    if metrics.brisque,                 results.(seqStructStr).(name).brisque.all(indexB)       = evaluateBRISQUE(predImage); end
                    if metrics.colDiff,                 results.(seqStructStr).(name).colDiff.all(indexB)       = evaluateColDiff(predImage,targetJ); end
                    
                    if dehazer.PredictsT && metrics.TError
                                                        results.(seqStructStr).(name).TError.all(indexB)        = evaluateTError(predT,targetT);
                    end
                    
                    if metrics.mic
                        ycbcr = rgb2ycbcr(predImage);
                                                        sequenceTotals.(name).runningMeans.Y(indexB)            = mean(ycbcr(:,:,1),'all');
                                                        sequenceTotals.(name).runningMeans.Cb(indexB)           = mean(ycbcr(:,:,2),'all');
                                                        sequenceTotals.(name).runningMeans.Cr(indexB)           = mean(ycbcr(:,:,3),'all');
                    end
                    
                end
                %% End Loop Dehazers
                prevJ = J;
                prevI = I;
                if metrics.TError, prevT = trueT; end
                
                frameTimeTotal = frameTimeTotal + toc(frameTic);
                
                remainingFrames = remainingFrames - 1;
                completedFrames = completedFrames + 1;
                
            end %% End Loop Frames
            
            seqTic = tic;
            textwaitbar(JSize,JSize,sprintf("Evaluating on %s",seqStructStr))
            
%             btcmDenom = evaluateBTCM(runningSSIM_J, runningSSIM_I)
            
            %assignin('base','totals',frameTotals);
            
            %% Post-Sequence Systems Loop 
            for dehazer = systems
                name = erase(class(dehazer),"Dehazer");
                
                % Deal with frame delay by passing in empty images
                if dehazer.FrameDelay > 0
                    % Currently assuming max delay of 1
                    
                    [predImage, predT, predA, timeImage, timeA] = dehazer.dehazeFrame([]);

                    targetJ = J;
                    if metrics.TError, targetT = trueT; end
                    
                    if metrics.mppsTotal,               results.(seqStructStr).(name).timeTotal.all(end)       = timeTotal; end
                    if metrics.mppsImage,               results.(seqStructStr).(name).timeImage.all(end)       = timeImage; end
                    
                    if dehazer.PredictsA
                        if metrics.mppsA,               results.(seqStructStr).(name).timeA.all(end)         = timeA; end
                        if metrics.AError,              results.(seqStructStr).(name).AError.all(end)        = evaluateAError(predA,trueA); end
                    end
                    
                    if metrics.psnr,                    results.(seqStructStr).(name).psnr.all(end)          = evaluatePSNR(predImage,targetJ); end
                    if metrics.ssim,                    results.(seqStructStr).(name).ssim.all(end)          = evaluateSSIM(predImage,targetJ); end
                    if metrics.fade,                    results.(seqStructStr).(name).fade.all(end)          = evaluateFADE(predImage); end
                    
                    if metrics.piqe,                    results.(seqStructStr).(name).piqe.all(end)          = evaluatePIQE(predImage); end
                    if metrics.brisque,                 results.(seqStructStr).(name).brisque.all(end)       = evaluateBRISQUE(predImage); end
                    if metrics.colDiff,                 results.(seqStructStr).(name).colDiff.all(end)       = evaluateColDiff(predImage,targetJ); end
                    
                    if dehazer.PredictsT && metrics.TError
                                                        results.(seqStructStr).(name).TError.all(end)        = evaluateTError(predT,targetT);
                    end
                    
                    if metrics.mic
                        ycbcr = rgb2ycbcr(predImage);
                                                        sequenceTotals.(name).runningMeans.Y(end)            = mean(ycbcr(:,:,1),'all');
                                                        sequenceTotals.(name).runningMeans.Cb(end)           = mean(ycbcr(:,:,2),'all');
                                                        sequenceTotals.(name).runningMeans.Cr(end)           = mean(ycbcr(:,:,3),'all');
                    end
                end
%                 if systems.(name).frameDelay > 0
%                     for fd = 0:systems.(name).frameDelay
%                         [predImage, predT, predA, timeImage, timeA, systems.(name).state] = systems.(name).function([],knowns,systems.(name).state);     
%                     end
%                 end

%%
                mp = results.(seqStructStr).props.megapixels;
                
                if metrics.mppsTotal
                            results.(seqStructStr).(name).mppsTotal.all             = mp./results.(seqStructStr).(name).timeTotal.all;
                            results.(seqStructStr).(name).mppsTotal.mean            = mean(results.(seqStructStr).(name).mppsTotal.all);
                            results.(seqStructStr).(name).mppsTotal.stdv            = std(results.(seqStructStr).(name).mppsTotal.all);
                            
                            results.(seqStructStr).(name).timeTotal.mean            = mean(results.(seqStructStr).(name).timeTotal.all);
                            results.(seqStructStr).(name).timeTotal.stdv            = std(results.(seqStructStr).(name).timeTotal.all);
                end
                
                if metrics.mppsImage
                            results.(seqStructStr).(name).mppsImage.all             = mp./results.(seqStructStr).(name).timeImage.all;
                            results.(seqStructStr).(name).mppsImage.mean            = mean(results.(seqStructStr).(name).mppsImage.all);
                            results.(seqStructStr).(name).mppsImage.stdv            = std(results.(seqStructStr).(name).mppsImage.all);
                            
                            results.(seqStructStr).(name).timeImage.mean            = mean(results.(seqStructStr).(name).timeImage.all);
                            results.(seqStructStr).(name).timeImage.stdv            = std(results.(seqStructStr).(name).timeImage.all);
                end
                
                if dehazer.PredictsA
                            results.(seqStructStr).(name).mppsA.all                 = mp./results.(seqStructStr).(name).timeA.all;
                            results.(seqStructStr).(name).mppsA.mean                = mean(results.(seqStructStr).(name).mppsA.all);
                            results.(seqStructStr).(name).mppsA.stdv                = std(results.(seqStructStr).(name).mppsA.all);
                            
                            results.(seqStructStr).(name).timeA.mean                = mean(results.(seqStructStr).(name).timeA.all);
                            results.(seqStructStr).(name).timeA.stdv                = std(results.(seqStructStr).(name).timeA.all);
                    
                    if metrics.AError
                            results.(seqStructStr).(name).AError.mean                = mean(results.(seqStructStr).(name).AError.all);
                            results.(seqStructStr).(name).AError.stdv                = std(results.(seqStructStr).(name).AError.all);
                    end
                else
                    if metrics.mppsA
                        results.(seqStructStr).(name).mppsA.all = NaN;
                        results.(seqStructStr).(name).mppsA.mean = NaN;
                        results.(seqStructStr).(name).mppsA.stdv = NaN;
                    end
                    if metrics.AError
                        results.(seqStructStr).(name).AError.all = NaN;
                        results.(seqStructStr).(name).AError.mean = NaN;
                        results.(seqStructStr).(name).AError.stdv = NaN;
                    end
                end
                
                if metrics.psnr
                        results.(seqStructStr).(name).psnr.mean                     =  mean(results.(seqStructStr).(name).psnr.all);
                        results.(seqStructStr).(name).psnr.stdv                     =  std(results.(seqStructStr).(name).psnr.all);
                end
                if metrics.ssim
                        results.(seqStructStr).(name).ssim.mean                     =  mean(results.(seqStructStr).(name).ssim.all);
                        results.(seqStructStr).(name).ssim.stdv                     =  std(results.(seqStructStr).(name).ssim.all);
                end
                if metrics.fade
                        results.(seqStructStr).(name).fade.mean                     =  mean(results.(seqStructStr).(name).fade.all);
                        results.(seqStructStr).(name).fade.stdv                     =  std(results.(seqStructStr).(name).fade.all);
                end
                
                
                if metrics.piqe
                        results.(seqStructStr).(name).piqe.mean                     =  mean(results.(seqStructStr).(name).piqe.all);
                        results.(seqStructStr).(name).piqe.stdv                     =  std(results.(seqStructStr).(name).piqe.all);
                end
                if metrics.brisque
                        results.(seqStructStr).(name).brisque.mean                     =  mean(results.(seqStructStr).(name).brisque.all);
                        results.(seqStructStr).(name).brisque.stdv                     =  std(results.(seqStructStr).(name).brisque.all);
                end
                if metrics.colDiff
                        results.(seqStructStr).(name).colDiff.mean                     =  mean(results.(seqStructStr).(name).colDiff.all);
                        results.(seqStructStr).(name).colDiff.stdv                     =  std(results.(seqStructStr).(name).colDiff.all);
                end
                
                if metrics.mic, results.(seqStructStr).(name).mic                   = evaluateMIC(sequenceTotals.(name).runningMeans,runningMeans);end
                
                if dehazer.PredictsT
                    if metrics.TError
                            results.(seqStructStr).(name).TError.mean                = mean(results.(seqStructStr).(name).TError.all);
                            results.(seqStructStr).(name).TError.stdv                = std(results.(seqStructStr).(name).TError.all);
                    end
                else
                    if metrics.TError
                            results.(seqStructStr).(name).TError.mean                = NaN;
                            results.(seqStructStr).(name).TError.stdv                = NaN;
                    end
                end
            end %% End Post-Sequence Systems Loop
            
            if ~isempty(autosave)
                save(autosave,'results');
                fprintf("Current results saved to %s\n",autosave);
            end
            
            seqOverheadTotal = seqOverheadTotal + toc(seqTic);
            remainingSequences = remainingSequences - 1;
            completedSequences = completedSequences + 1;
            
            frameTimeEst = (frameTimeTotal / completedFrames) * remainingFrames;
            seqOverheadEst = (seqOverheadTotal / completedSequences) * remainingSequences;
            
            timeEst = seconds(frameTimeEst+seqOverheadEst);
            timeEst.Format = "hh:mm:ss";
            
            fprintf("%d of %d frames completed\n",completedFrames,completedFrames+remainingFrames);
            fprintf("Estimated time remaining: %s\n",timeEst);
        end %% End Sequence Loop
    end %% End Driver Loop

end
%%

function [results, systems, metrics, subset, autosave] = validateInputStruct(params)
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
        systems = [BermanDehazer, CaiDehazer, ChenDehazer, HeDehazer, LiDehazer, TarelDehazer, TilburyDehazer, TsaiDehazer, ZhuDehazer];
    elseif isa(params.Systems, "BaseDehazer")
        systems = params.Systems;
    elseif isa(params.Systems, "string")
        systems = [];
        for s = params.Systems
            if ismember("BaseDehazer", superclasses("s"))
                systems(end+1) = feval(s)
            end
        end
    else
        systems = [BermanDehazer, CaiDehazer, ChenDehazer, HeDehazer, LiDehazer, TarelDehazer, TilburyDehazer, TsaiDehazer, ZhuDehazer];
    end
    
    allMetrics = ["mppsImage" "mppsTotal" "mppsA" "AError" "psnr" "ssim" "fade" "colDiff" "mic" "TError", "brisque", "piqe"];
%     useMetrics = ["mppsImage" "mppsTotal" "mppsA" "AError" "psnr" "ssim" "fade" "colDiff" "mic" "TError", "brisque", "piqe"];
    
    if ~isfield(params,'Metrics') || isempty(params.Metrics)
        textrics = allMetrics;
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
%         results.metrics = textrics;
    else
%         results.names = union(params.Results.names,string(fieldnames(systems))');
%         results.metrics = union(params.Results.metrics,textrics);
    end
    
    if ~isfield(params,'Autosave') || isempty(params.Autosave) || ~isstring(params.Autosave)
        autosave = [];
    else
        if ~endsWith(params.Autosave,".mat"), params.Autosave = append(params.Autosave,".mat"); end
        autosave = fullfile(fileparts(mfilename('fullpath')),params.Autosave);
    end
%     fullCols = [
%         31 119 180;
%         255 127 14;
%         44 160 44;
%         214 39 40;
%         148 103 189;
%         140 86 75;
%         227 119 194;
%         127 127 127;
%         188 189 34;
%         23 190 207] / 255;

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