function results = evaluateDehazers(path, params)
    warning('on','all');
 
    if ~exist('params','var'), params = []; end
    [results, systems, metrics, sequences, autosave, sample, seed] = validateInputStruct(params);   
    
    rng(seed);
    
    %% Preliminary look to find total number of frames and deal with sampling
       
    dateFolders = dir(path);
    dateFolders = dateFolders([dateFolders.isdir]);
    
    remainingFrames = 0;
    remainingSequences = 0;
    
    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        if ~isfile(fullfile(datePath,'calib_cam_to_cam.txt')), continue; end
        
        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        for j = 3:length(sequenceFolders)
            
            sequenceName = sequenceFolders(j).name;
            sequencePath = fullfile(datePath, sequenceName);
            seqStr = ['S' sequenceName(1:end-5)];
            
            if ~any(strcmp(sequences,sequenceName)), continue;  end
            
            JPath = fullfile(sequencePath,'image_02','data');
            JSize = length(dir(JPath))-2;
            
            if sample>0
                numFrames =  min(sample,JSize);
                frameStarts.(seqStr) = randi(JSize - numFrames + 1) - 1;
            else
                numFrames = JSize;
                frameStarts.(seqStr) = 0;
            end
            
            frameLengths.(seqStr) = numFrames;
            
            remainingSequences = remainingSequences + 1;
            
            remainingFrames = remainingFrames + numFrames;
        end
    end
    
%     frameTimeEst = 0;
%     seqOverheadEst = 0;
%     timeEst = 0;
    
    frameTimeTotal = 0;
    seqOverheadTotal = 0;
    
    completedFrames = 0;
    completedSequences = 0;
    
    if sample>0
        fprintf("Evaluating %d dehazing systems on %d sequences using a maximum of %d frames per sequence. Total frames: %d\r\n", length(systems), remainingSequences, sample, remainingFrames);
    else
        fprintf("Evaluating %d dehazing systems on %d sequences. Total frames: %d\r\n", length(systems), remainingSequences, remainingFrames);
    end
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
        K = calib.P_rect{3};
        
        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        seqOverheadTotal = seqOverheadTotal + toc(seqTic);
        
        %% Loop Sequence
        for j = 3:length(sequenceFolders)
            
            seqTic = tic;
            sequenceName = sequenceFolders(j).name;
            sequencePath = fullfile(datePath, sequenceName);
            seqStr = ['S' sequenceName(1:end-5)];
            
            if sequences~="all"
                if ~any(strcmp(sequences,sequenceName)), continue;  end
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
            
            trueA = reshape(hazeProps.A, [1 1 3]);
            
            JPath = fullfile(sequencePath,'image_02','data');
            
            start = frameStarts.(seqStr);
            JSize = frameLengths.(seqStr);
            
            mPath = fullfile(sequencePath,'sky_mask','default_refined','img'); 
            
            IPath = fullfile(sequencePath,'haze','image');
            
            if metrics.TError, TPath = fullfile(sequencePath,'haze','transmission'); end
            
            results.(seqStr).props.visibility = hazeProps.visibility; % Will be useful for comparing how different dehazing methods do on different fog densities
            results.(seqStr).props.evaluatedFrames = [start, JSize-1];
            
            if hazeProps.visibility<0.4
                results.(seqStr).props.group = 'dense';
            elseif hazeProps.visibility<0.7
                results.(seqStr).props.group = 'moderate';
            else
                results.(seqStr).props.group = 'light';
            end
            
            seqOverheadTotal = seqOverheadTotal + toc(seqTic);
            
            % Doing it frameX->system1,system2,... rather than systemX->frame1,frame2,...
            % to save having to load the images multiple times.
            % the only downside is needing to store state and 
            % results for each system which isn't a big deal
            % it also means that for delayed methods the previous J's and T's
            % need to be kept for as long as the delay, and that saving video is more of a hassle
            % To put it simply, doing it this way is slightly faster at the expense of requiring more memory.
            %% Loop Frames
            prevJ = [];
            prevPrevJ = [];
            for frameID = 0:JSize-1
                frameTic = tic;
                frameFileID = frameID + start;
                framePath = num2str(frameFileID,'%010.f')+".png";

                J = imread(fullfile(JPath,framePath));
                J = im2double(J);
                
                % Apply sky mask
                m = imread(fullfile(mPath,framePath));
                m = im2double(m);
                m = repmat(m,1,1,3);
                J = m.*trueA + (1-m).*J;
                %%%%%

                I = imread(fullfile(IPath,framePath));
                I = im2double(I);
                                
                if frameID==0
                    results.(seqStr).props.megapixels = (size(J,1)*size(J,2))/(10^6);
                    results.(seqStr).props.height = size(J,1);
                    results.(seqStr).props.width = size(J,2);
                    
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
                
                textwaitbar(frameID,JSize,sprintf("Evaluating on %s",seqStr))
                
                %% Loop Dehazers
                for dehazer = systems
                    name = dehazer.Name;
                    
                    if frameID==0
                        
                        if metrics.timeTotal,           results.(seqStr).(name).timeTotal.all         = zeros(1,JSize+dehazer.FrameDelay); end
                        
                        if dehazer.PredictsA
                            if metrics.timeA,           results.(seqStr).(name).timeA.all             = zeros(1,JSize); end
                            if metrics.AError,          results.(seqStr).(name).AError.all            = zeros(1,JSize); end
                        end
                            
                        if metrics.timeImage,           results.(seqStr).(name).timeImage.all         = zeros(1,JSize); end
                        
                        if metrics.psnr,                results.(seqStr).(name).psnr.all              = zeros(1,JSize); end
                        if metrics.ssim,                results.(seqStr).(name).ssim.all              = zeros(1,JSize); end
                        if metrics.msssim,              results.(seqStr).(name).msssim.all            = zeros(1,JSize); end
                        if metrics.fade,                results.(seqStr).(name).fade.all              = zeros(1,JSize); end
                        
                        if metrics.piqe,                results.(seqStr).(name).piqe.all              = zeros(1,JSize); end
                        if metrics.niqe,                results.(seqStr).(name).niqe.all              = zeros(1,JSize); end
                        if metrics.brisque,             results.(seqStr).(name).brisque.all           = zeros(1,JSize); end
                        if metrics.colDiff,             results.(seqStr).(name).colDiff.all           = zeros(1,JSize); end
                        
                        
                        if metrics.mic
                                                        sequenceTotals.(name).runningMeans.Y                = zeros(1,JSize);
                                                        sequenceTotals.(name).runningMeans.Cb               = zeros(1,JSize);
                                                        sequenceTotals.(name).runningMeans.Cr               = zeros(1,JSize);
                        end
                        if metrics.speedqa
                                                        results.(seqStr).(name).speedqa.all           = zeros(1,JSize-1);

                        end
                    
                        if dehazer.PredictsT && metrics.TError
                                                        results.(seqStr).(name).TError.all            = zeros(1,JSize);
                        end
                        
                        knowns.K = K;
                        knowns.imgSize = [size(J,1) size(J,2)];
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
                    indexC = frameID-dehazer.FrameDelay;
                    
                    if frameID < dehazer.FrameDelay
                        if metrics.timeTotal,           results.(seqStr).(name).timeTotal.all(indexA)    = timeTotal; end
                        continue
                    end
                    
                    
                    if dehazer.FrameDelay == 0
                        targetJ = J;
                        targetPrevJ = prevJ;
                        if metrics.TError, targetT = trueT; end
                    else
                        targetJ = prevJ;
                        targetPrevJ = prevPrevJ;
                        if metrics.TError, targetT = prevT; end
                    end
                    
                    if metrics.speedqa 
                        if frameID>dehazer.FrameDelay
                                                        results.(seqStr).(name).speedqa.all(indexC)      = evaluateSpEEDQA(predImage, prevPreds.(name), targetJ, targetPrevJ);
                        end
                        prevPreds.(name) = predImage;
                    end
                    
                    if metrics.timeTotal,               results.(seqStr).(name).timeTotal.all(indexA)     = timeTotal; end
                    if metrics.timeImage,               results.(seqStr).(name).timeImage.all(indexB)     = timeImage; end
                    
                    if dehazer.PredictsA
                        if metrics.timeA,               results.(seqStr).(name).timeA.all(indexB)         = timeA; end
                        if metrics.AError,              results.(seqStr).(name).AError.all(indexB)        = evaluateAError(predA,trueA); end
                    end
                    
                    if metrics.psnr,                    results.(seqStr).(name).psnr.all(indexB)          = evaluatePSNR(predImage,targetJ); end
                    if metrics.ssim,                    results.(seqStr).(name).ssim.all(indexB)          = evaluateSSIM(predImage,targetJ); end
                    if metrics.msssim,                  results.(seqStr).(name).msssim.all(indexB)        = evaluateMSSSIM(predImage,targetJ); end
                    if metrics.fade,                    results.(seqStr).(name).fade.all(indexB)          = evaluateFADE(predImage); end
                    
                    if metrics.piqe,                    results.(seqStr).(name).piqe.all(indexB)          = evaluatePIQE(predImage); end
                    if metrics.niqe,                    results.(seqStr).(name).niqe.all(indexB)          = evaluateNIQE(predImage); end
                    if metrics.brisque,                 results.(seqStr).(name).brisque.all(indexB)       = evaluateBRISQUE(predImage); end
                    if metrics.colDiff,                 results.(seqStr).(name).colDiff.all(indexB)       = evaluateColDiff(predImage,targetJ); end
                    
                    if dehazer.PredictsT && metrics.TError
                                                        results.(seqStr).(name).TError.all(indexB)        = evaluateTError(predT,targetT);
                    end
                    
                    if metrics.mic
                        ycbcr = rgb2ycbcr(predImage);
                        sequenceTotals.(name).runningMeans.Y(indexB)            = mean(ycbcr(:,:,1),'all');
                        sequenceTotals.(name).runningMeans.Cb(indexB)           = mean(ycbcr(:,:,2),'all');
                        sequenceTotals.(name).runningMeans.Cr(indexB)           = mean(ycbcr(:,:,3),'all');
                    end
                    
                end
                %% End Loop Dehazers
                prevPrevJ = prevJ; % Dumb hacky way of doing this.
                prevJ = J;
                prevI = I;
%                 prevPred = predImage;
                
                if metrics.TError, prevT = trueT; end
                
                frameTimeTotal = frameTimeTotal + toc(frameTic);
                
                remainingFrames = remainingFrames - 1;
                completedFrames = completedFrames + 1;
                
            end %% End Loop Frames
            
            seqTic = tic;
            
            
%             btcmDenom = evaluateBTCM(runningSSIM_J, runningSSIM_I)
            
            %assignin('base','totals',frameTotals);
            
            %% Post-Sequence Systems Loop 
            for dehazer = systems
                name = dehazer.Name;
                
                % Deal with frame delay by passing in empty images
                if dehazer.FrameDelay > 0
                    % Currently assuming max delay of 1
                    
                    [predImage, predT, predA, timeImage, timeA] = dehazer.dehazeFrame([]);

                    targetJ = J;
                    targetPrevJ = prevJ;
                    
                    if metrics.TError, targetT = trueT; end
                    
                    if metrics.speedqa,                 results.(seqStr).(name).speedqa.all(indexC)      = evaluateSpEEDQA(predImage, dehazer.SequenceState.prevFrame, targetJ, targetPrevJ); end
                    
                    if metrics.timeTotal,               results.(seqStr).(name).timeTotal.all(end)       = timeTotal; end
                    if metrics.timeImage,               results.(seqStr).(name).timeImage.all(end)       = timeImage; end
                    
                    if dehazer.PredictsA
                        if metrics.timeA,               results.(seqStr).(name).timeA.all(end)         = timeA; end
                        if metrics.AError,              results.(seqStr).(name).AError.all(end)        = evaluateAError(predA,trueA); end
                    end
                    
                    if metrics.psnr,                    results.(seqStr).(name).psnr.all(end)          = evaluatePSNR(predImage,targetJ); end
                    if metrics.ssim,                    results.(seqStr).(name).ssim.all(end)          = evaluateSSIM(predImage,targetJ); end
                    if metrics.msssim,                  results.(seqStr).(name).msssim.all(end)        = evaluateMSSSIM(predImage,targetJ); end
                    if metrics.fade,                    results.(seqStr).(name).fade.all(end)          = evaluateFADE(predImage); end
                    
                    if metrics.piqe,                    results.(seqStr).(name).piqe.all(end)          = evaluatePIQE(predImage); end
                    if metrics.niqe,                    results.(seqStr).(name).niqe.all(end)          = evaluateNIQE(predImage); end
                    if metrics.brisque,                 results.(seqStr).(name).brisque.all(end)       = evaluateBRISQUE(predImage); end
                    if metrics.colDiff,                 results.(seqStr).(name).colDiff.all(end)       = evaluateColDiff(predImage,targetJ); end
                    
                    if dehazer.PredictsT && metrics.TError
                                                        results.(seqStr).(name).TError.all(end)        = evaluateTError(predT,targetT);
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

                if metrics.speedqa
                    results.(seqStr).(name).speedqa.mean                     =  mean(results.(seqStr).(name).speedqa.all,'omitnan');
                    results.(seqStr).(name).speedqa.stdv                     =  std(results.(seqStr).(name).psnr.all,'omitnan');
                end
                
                % megapixels, seconds-per-megapixel, and
                % megapixels-per-second stuff will be dealt with elsewhere
%                 mp = results.(seqStr).props.megapixels;
                
                if metrics.timeTotal
                        if dehazer.FrameDelay
                            % Distribute the extra 'set-up' frame time accross the the rest of the times. Won't change stdev
                            results.(seqStr).(name).timeTotal.all = results.(seqStr).(name).timeTotal.all(2:end) + results.(seqStr).(name).timeTotal.all(1)./JSize;
                        end

                            results.(seqStr).(name).timeTotal.mean            = mean(results.(seqStr).(name).timeTotal.all);
                            results.(seqStr).(name).timeTotal.stdv            = std(results.(seqStr).(name).timeTotal.all);
                            
                end
                
                if metrics.timeImage 
                            results.(seqStr).(name).timeImage.mean            = mean(results.(seqStr).(name).timeImage.all);
                            results.(seqStr).(name).timeImage.stdv            = std(results.(seqStr).(name).timeImage.all);
                end
                
                if dehazer.PredictsA
                    if metrics.timeA
                            results.(seqStr).(name).timeA.mean                = mean(results.(seqStr).(name).timeA.all);
                            results.(seqStr).(name).timeA.stdv                = std(results.(seqStr).(name).timeA.all);
                    end
                    
                    if metrics.AError
                            results.(seqStr).(name).AError.mean                = mean(results.(seqStr).(name).AError.all);
                            results.(seqStr).(name).AError.stdv                = std(results.(seqStr).(name).AError.all);
                    end
                else
                    if metrics.timeA
                        results.(seqStr).(name).timeA.mean = NaN;
                    end
                    if metrics.AError
                        results.(seqStr).(name).AError.mean = NaN;
                        results.(seqStr).(name).AError.stdv = NaN;
                    end
                end
                
                if metrics.psnr
                        results.(seqStr).(name).psnr.mean                     =  mean(results.(seqStr).(name).psnr.all);
                        results.(seqStr).(name).psnr.stdv                     =  std(results.(seqStr).(name).psnr.all);
                end
                if metrics.ssim
                        results.(seqStr).(name).ssim.mean                     =  mean(results.(seqStr).(name).ssim.all);
                        results.(seqStr).(name).ssim.stdv                     =  std(results.(seqStr).(name).ssim.all);
                end
                if metrics.msssim
                        results.(seqStr).(name).msssim.mean                     =  mean(results.(seqStr).(name).msssim.all);
                        results.(seqStr).(name).msssim.stdv                     =  std(results.(seqStr).(name).msssim.all);
                end
                if metrics.fade
                        results.(seqStr).(name).fade.mean                     =  mean(results.(seqStr).(name).fade.all);
                        results.(seqStr).(name).fade.stdv                     =  std(results.(seqStr).(name).fade.all);
                end
                
                
                if metrics.piqe
                        results.(seqStr).(name).piqe.mean                     =  mean(results.(seqStr).(name).piqe.all);
                        results.(seqStr).(name).piqe.stdv                     =  std(results.(seqStr).(name).piqe.all);
                end
                if metrics.niqe
                        results.(seqStr).(name).niqe.mean                     =  mean(results.(seqStr).(name).niqe.all);
                        results.(seqStr).(name).niqe.stdv                     =  std(results.(seqStr).(name).niqe.all);
                end
                if metrics.brisque
                        results.(seqStr).(name).brisque.mean                     =  mean(results.(seqStr).(name).brisque.all);
                        results.(seqStr).(name).brisque.stdv                     =  std(results.(seqStr).(name).brisque.all);
                end
                if metrics.colDiff
                        results.(seqStr).(name).colDiff.mean                     =  mean(results.(seqStr).(name).colDiff.all);
                        results.(seqStr).(name).colDiff.stdv                     =  std(results.(seqStr).(name).colDiff.all);
                end
                
                if metrics.mic, results.(seqStr).(name).mic                   = evaluateMIC(sequenceTotals.(name).runningMeans,runningMeans);end
                
                if dehazer.PredictsT
                    if metrics.TError
                            results.(seqStr).(name).TError.mean                = mean(results.(seqStr).(name).TError.all);
                            results.(seqStr).(name).TError.stdv                = std(results.(seqStr).(name).TError.all);
                    end
                else
                    if metrics.TError
                            results.(seqStr).(name).TError.mean                = NaN;
                    end
                end
                
                textwaitbar(JSize,JSize,sprintf("Evaluating on %s",seqStr))
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
            timeEst.Format = "dd:hh:mm:ss";
            
            timeEstString = string(timeEst);
            if strlength(timeEstString)==8 % So Matlabs duration apparently ignores the days format completely when duration is less than a day. Doesn't do the same for hours or minutes though. If I specify an output format it should follow that output format not whatever matlab fancies. Rant over. 
                timeEstString = append("00:",timeEstString);
            end
            
            stringDay = extractBetween(timeEstString,1,2);
            stringHour = extractBetween(timeEstString,4,5);
            stringMin = extractBetween(timeEstString,7,8);
            stringSec = extractBetween(timeEstString,10,11);
            
            fprintf("%d of %d frames completed\n",completedFrames,completedFrames+remainingFrames);
            fprintf("Estimated Time Remaining: %s Days %s Hours %s Minutes and %s Seconds\r\n",stringDay,stringHour,stringMin,stringSec);
        end %% End Sequence Loop
    end %% End Driver Loop

end
%%

function [results, systems, metrics, sequences, autosave, sample, seed] = validateInputStruct(params)
    if isempty(params) || ~isstruct(params)
        params = struct();
    end

    if ~isfield(params,'Sequences') || isempty(params.Sequences)
        % Use only the dataset's subset of KITTI data 
        fid = fopen(fullfile(fileparts(mfilename('fullpath')),"..","haze-subset.txt"));
        S = textscan(fid,"%s");
        sequences = string(S{1});
        fclose(fid);

    else
        sequences = params.Sequences;
    end
    
    
    if ~isfield(params,'Seed') || isempty(params.Seed)
        seed = 0;
    else
        seed = params.Seed;
    end
    
    
    if ~isfield(params,'Sample') || isempty(params.Sample)
        sample = 0;
    else
        sample = params.Sample;
    end
    
    
    if ~isfield(params,'Systems') || isempty(params.Systems)
        % Unoptimised TilburyDehazer parameters
%         unoptTilbury = TilburyDehazer.initial;
        % Optimised TilburyDehazer parameters (default values in class)
        % Currently excluded, need to re-optimise...
%         optTilbury = TilburyDehazer;
%         optTilbury.rename("TilburyOptimised"); % Rename to avoid name clash
        
        systems = [BermanDehazer, CaiDehazer, ChenDehazer, HeDehazer, LiDehazer, TarelDehazer, TsaiDehazer, ZhuDehazer, TilburyDehazer.bestLocal, TilburyDehazer.bestGlobal, TilburyDehazer.bestOpening];
    elseif isa(params.Systems, "BaseDehazer")
        systems = params.Systems;
    elseif isa(params.Systems, "string")
        systems = [];
        for s = params.Systems
            if ismember("BaseDehazer", superclasses("s"))
                systems(end+1) = feval(s);
            end
        end
    else
        systems = [BermanDehazer, CaiDehazer, ChenDehazer, HeDehazer, LiDehazer, TarelDehazer, TilburyDehazer, TsaiDehazer, ZhuDehazer];
    end
    
    allMetrics = ["timeImage" "timeTotal" "timeA" "AError" "psnr" "msssim" "ssim" "fade" "colDiff" "mic" "TError", "brisque", "piqe", "niqe", "speedqa"];
%     useMetrics = ["timeImage" "timeTotal" "timeA" "AError" "psnr" "ssim" "fade" "colDiff" "mic" "TError", "brisque", "piqe"];
    
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
        results = params.Results;
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

function q = evaluateMSSSIM(pred, trueImage)
    if any(size(pred,[1 2])~=size(trueImage,[1 2]))
       wind = centerCropWindow2d(size(trueImage),size(pred));
       trueImage = imcrop(trueImage, wind);
    end
    
    q = mean(multissim(pred,trueImage));
end

function q = evaluateSSIM(pred, trueImage)
    if any(size(pred,[1 2])~=size(trueImage,[1 2]))
       wind = centerCropWindow2d(size(trueImage),size(pred));
       trueImage = imcrop(trueImage, wind);
    end
    
    ssimr = ssim(pred(:,:,1),trueImage(:,:,1));
    ssimg = ssim(pred(:,:,2),trueImage(:,:,2));
    ssimb = ssim(pred(:,:,3),trueImage(:,:,3));
            
    q = (ssimr+ssimg+ssimb)/3;
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

function q = evaluateNIQE(pred)
    q = niqe(pred);
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

function q = evaluateSpEEDQA(pred, prevPred, trueImage, prevTrue)
    % Default parameters
    blk_speed = 5;
    sigma_nsq = 0.1;
    down_size = 4;
    window = fspecial('gaussian', 7, 7/6);
    window = window/sum(sum(window));
    
    pred = rgb2ycbcr(pred);
    pred = pred(:,:,1);
    
    prevPred = rgb2ycbcr(prevPred);
    prevPred = prevPred(:,:,1);
    
    trueImage = rgb2ycbcr(trueImage);
    trueImage = trueImage(:,:,1);
    
    prevTrue = rgb2ycbcr(prevTrue);
    prevTrue = prevTrue(:,:,1);

    [speed_s, ~, speed_t, ~] = Single_Scale_Video_SPEED(prevTrue, trueImage, prevPred, pred, down_size, window, blk_speed, sigma_nsq);  
    
    q = speed_s*speed_t;
end

% function q = evaluateTCM(ssimA, ssimB)
%     ssimCorr = corrcoef(ssimA,ssimB);
%     q = ssimCorr(1,2);
% end


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