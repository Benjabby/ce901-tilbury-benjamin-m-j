function results = evaluateDehaze(path, subset, systems)
    warning('on','all');
    if ~exist('path','var'), path = fullfile(mfilename('fullpath'),"..","..","haze-video-dataset","dataset"); end

    if ~exist('subset','var')
       % To do. Load default subset file 
       % For now set to 'all'
       subset = "all";
    end

    if ~exist('systems','var')
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
    results.metrics = ["mppsImage" "mppsTotal" "mppsA" "AError" "psnr" "ssim" "fade" "colDiff" "ic" "tcm" "btcm" "TError"];
    
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
        if subset~="all"
            % do subset
        end
        
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
            
            JPath = fullfile(sequencePath,'image_02','data');
            JSize = length(dir(JPath))-2;
            
            IPath = fullfile(sequencePath,'haze','image');
%             ISize = length(dir(IPath))-2;
            
            TPath = fullfile(sequencePath,'haze','transmission');
%             TSize = length(dir(IPath))-2;
            
            
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
                    runningSSIM_I = zeros(1,JSize-1);
                    runningMeans.H = zeros(1,JSize);
                    runningMeans.S = zeros(1,JSize);
                    runningMeans.V = zeros(1,JSize);
%                     runningSSIM_J = zeros(1,JSize-1);
                else
                    runningSSIM_I(frameID) = evaluateSSIM(prevI,I);
%                     runningSSIM_J(frameID) = evaluateSSIM(prevJ,J);
%                     frameDiff = prevJ-J;
                    corrJ = corrcoef(J,prevJ);
                    corrJ = corrJ(1,2);
                end
                
                [hue, sat, val] = rgb2hsv(J);
                runningMeans.H(frameID+1) = mean(hue,'all');
                runningMeans.S(frameID+1) = mean(sat,'all');
                runningMeans.V(frameID+1) = mean(val,'all');

                trueT = imread(fullfile(TPath,framePath));
                trueT = im2double(trueT);
                
                textwaitbar(frameID,JSize,sprintf("Evaluating on %s",seqStructStr))
                
                for name = string(fieldnames(systems))'
                    if name == "prepped", continue; end
                    
                    if frameID==0
                        
                        sequenceTotals.(name).timeTotal    = 0;
                        
                        if systems.(name).predictsA
                            sequenceTotals.(name).timeA    = 0;
                            sequenceTotals.(name).AError   = 0;
                        end
                            
                        sequenceTotals.(name).timeImage    = 0;
                        
                        sequenceTotals.(name).psnr         = 0;
                        sequenceTotals.(name).ssim         = 0;
                        sequenceTotals.(name).fade         = 0;
                        
                        sequenceTotals.(name).piqe         = 0;
                        sequenceTotals.(name).brisque      = 0;
                        sequenceTotals.(name).colDiff      = 0;
                        
                        sequenceTotals.(name).tcm          = 0;
                        sequenceTotals.(name).runningSSIM  = zeros(1,JSize-1);
                        
                        sequenceTotals.(name).runningMeans.H = zeros(1,JSize);
                        sequenceTotals.(name).runningMeans.S = zeros(1,JSize);
                        sequenceTotals.(name).runningMeans.V = zeros(1,JSize);
                    
                        if systems.(name).predictsT
                           sequenceTotals.(name).TError    = 0;
                        end
                        
%                         if systems.(name).frameDelay > 1
%                             sequenceTotals.JBuffer = zeros([size(J) systems.(name).frameDelay])
%                             sequenceTotals.TBuffer = zeros([size(trueT) systems.(name).frameDelay])
%                         end
                        
                        systems.(name).state = struct();
                    end
                    
                    totalTic = tic;
                    [predImage, predT, predA, timeImage, timeA, systems.(name).state] = systems.(name).function(I,knowns,systems.(name).state);
                    timeTotal = toc(totalTic);
                    
                    if frameID < systems.(name).frameDelay
                        sequenceTotals.(name).timeTotal    = sequenceTotals.(name).timeTotal      + timeTotal;
                        continue
                    end
                    
                    predImage = min(max(predImage,0),1);
                    predT = min(max(predT,0),1);
                    predA = min(max(predA,0),1);
                    
                    if systems.(name).frameDelay == 0
                        targetJ = J;
                        targetT = trueT;
                        if frameID ~= 0
                            targetCorrJ = corrJ;
                        end
                    else
                        targetJ = prevJ;
                        targetT = prevT;
                        if frameID > systems.(name).frameDelay
                            targetCorrJ = prevCorrJ;
                        end
                    end
                    
                    sequenceTotals.(name).timeTotal    = sequenceTotals.(name).timeTotal      + timeTotal;
                    sequenceTotals.(name).timeImage    = sequenceTotals.(name).timeImage      + timeImage;
                    
                    if systems.(name).predictsA
                        sequenceTotals.(name).timeA    = sequenceTotals.(name).timeA          + timeA;
                        sequenceTotals.(name).AError   = sequenceTotals.(name).AError         + evaluateAError(predA,trueA);
                    end
                    
                    sequenceTotals.(name).psnr         = sequenceTotals.(name).psnr           + evaluatePSNR(predImage,targetJ);
                    sequenceTotals.(name).ssim         = sequenceTotals.(name).ssim           + evaluateSSIM(predImage,targetJ);
                    sequenceTotals.(name).fade         = sequenceTotals.(name).fade           + evaluateFADE(predImage);
                    
                    sequenceTotals.(name).piqe         = sequenceTotals.(name).piqe           + evaluatePIQE(predImage); 
                    sequenceTotals.(name).brisque      = sequenceTotals.(name).brisque        + evaluateBRISQUE(predImage);
                    sequenceTotals.(name).colDiff      = sequenceTotals.(name).colDiff        + evaluateColDiff(predImage,targetJ);
                    
                    if systems.(name).predictsT
                        sequenceTotals.(name).TError   = sequenceTotals.(name).TError         + evaluateTError(predT,targetT);
                    end
                    
                    [hue,sat,val] = rgb2hsv(predImage);
                    
                    sequenceTotals.(name).runningMeans.H(frameID+1-systems.(name).frameDelay) = mean(hue,'all');
                    sequenceTotals.(name).runningMeans.S(frameID+1-systems.(name).frameDelay) = mean(sat,'all');
                    sequenceTotals.(name).runningMeans.V(frameID+1-systems.(name).frameDelay) = mean(val,'all');
                                        
                    if frameID > systems.(name).frameDelay
                        sequenceTotals.(name).tcm      = sequenceTotals.(name).tcm            + evaluateTCM(predImage, sequenceTotals.(name).prevPred, targetCorrJ);
                        sequenceTotals.(name).runningSSIM(frameID-systems.(name).frameDelay) = evaluateSSIM(sequenceTotals.(name).prevPred, predImage);
                    end
                    
                    sequenceTotals.(name).prevPred = predImage;
                end
                
                prevJ = J;
                prevI = I;
                prevT = trueT;
                if frameID>0
                    prevCorrJ = corrJ;
                end
            end
            textwaitbar(JSize,JSize,sprintf("Evaluating on %s",seqStructStr))
                
%             btcmDenom = evaluateBTCM(runningSSIM_J, runningSSIM_I)
            
            %assignin('base','totals',frameTotals);
            for name = string(fieldnames(systems))'
                if name == "prepped", continue; end
                
                %% Deal with frame delay by passing in empty images
                if systems.(name).frameDelay > 0
                    % Currently assuming delay of one
                    targetJ = J;
                    targetT = trueT;
                    targetCorrJ = corrJ;
                        
                    sequenceTotals.(name).timeTotal    = sequenceTotals.(name).timeTotal      + timeTotal;
                    sequenceTotals.(name).timeImage    = sequenceTotals.(name).timeImage      + timeImage;
                    
                    if systems.(name).predictsA
                        sequenceTotals.(name).timeA    = sequenceTotals.(name).timeA          + timeA;
                        sequenceTotals.(name).AError   = sequenceTotals.(name).AError         + evaluateAError(predA,trueA);
                    end
                    
                    sequenceTotals.(name).psnr         = sequenceTotals.(name).psnr           + evaluatePSNR(predImage,targetJ);
                    sequenceTotals.(name).ssim         = sequenceTotals.(name).ssim           + evaluateSSIM(predImage,targetJ);
                    sequenceTotals.(name).fade         = sequenceTotals.(name).fade           + evaluateFADE(predImage);
                    
                    sequenceTotals.(name).piqe         = sequenceTotals.(name).piqe           + evaluatePIQE(predImage); 
                    sequenceTotals.(name).brisque      = sequenceTotals.(name).brisque        + evaluateBRISQUE(predImage);
                    sequenceTotals.(name).colDiff      = sequenceTotals.(name).colDiff        + evaluateColDiff(predImage,targetJ);
                    
                    if systems.(name).predictsT
                        sequenceTotals.(name).TError   = sequenceTotals.(name).TError         + evaluateTError(predT,targetT);
                    end
                    
                    [hue,sat,val] = rgb2hsv(predImage);
                    
                    sequenceTotals.(name).runningMeans.H(end) = mean(hue,'all');
                    sequenceTotals.(name).runningMeans.S(end) = mean(sat,'all');
                    sequenceTotals.(name).runningMeans.V(end) = mean(val,'all');
                                        
                    sequenceTotals.(name).tcm      = sequenceTotals.(name).tcm            + evaluateTCM(predImage, sequenceTotals.(name).prevPred, targetCorrJ);
                    sequenceTotals.(name).runningSSIM(end) = evaluateSSIM(sequenceTotals.(name).prevPred, predImage);
                end
%                 if systems.(name).frameDelay > 0
%                     for fd = 0:systems.(name).frameDelay
                        [predImage, predT, predA, timeImage, timeA, systems.(name).state] = systems.(name).function([],knowns,systems.(name).state);
                        
%                     end
%                 end

%%
                
                mp = results.(seqStructStr).megapixels;
                results.(seqStructStr).(name).mppsImage     = (1/(sequenceTotals.(name).timeImage/JSize))*mp;
                results.(seqStructStr).(name).mppsTotal     = (1/(sequenceTotals.(name).timeTotal/JSize))*mp;
                
                if systems.(name).predictsA
                    results.(seqStructStr).(name).mppsA     = (1/(sequenceTotals.(name).timeA/JSize))*mp;
                    results.(seqStructStr).(name).AError    = sequenceTotals.(name).AError/JSize;
                else
                    results.(seqStructStr).(name).mppsA     = NaN;
                    results.(seqStructStr).(name).AError    = NaN;
                end
                
                results.(seqStructStr).(name).psnr          =  sequenceTotals.(name).psnr / JSize;
                results.(seqStructStr).(name).ssim          =  sequenceTotals.(name).ssim / JSize;
                results.(seqStructStr).(name).fade          =  sequenceTotals.(name).fade / JSize;

                results.(seqStructStr).(name).colDiff       =  sequenceTotals.(name).colDiff / JSize;

                results.(seqStructStr).(name).ic            = evaluateIC(sequenceTotals.(name).runningMeans,runningMeans);
                results.(seqStructStr).(name).tcm           = sequenceTotals.(name).tcm / (JSize-1);
                results.(seqStructStr).(name).btcm          = evaluateBTCM(sequenceTotals.(name).runningSSIM, runningSSIM_I);
                
                if systems.(name).predictsT
                    results.(seqStructStr).(name).TError    = sequenceTotals.(name).TError/JSize;
                else
                    results.(seqStructStr).(name).TError    = NaN;
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
function q = evaluateIC(predMeans, trueMeans)
    corrH = corrcoef(predMeans.H, trueMeans.H);
    corrH = corrH(1,2);
    corrS = corrcoef(predMeans.S, trueMeans.S);
    corrS = corrS(1,2);
    corrV = corrcoef(predMeans.V, trueMeans.V);
    corrV = corrV(1,2);
    
    q = (corrH+corrS+corrV)/3;
end

function q = evaluateTCM(pred, prevPred, corrTrue)
    corrPred = corrcoef(pred,prevPred);
    corrPred = corrPred(1,2);
    q = (corrPred-corrTrue).^2;
end

function q = evaluateBTCM(ssimPred,ssimHaze)
    ssimCorr = corrcoef(ssimPred,ssimHaze);
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