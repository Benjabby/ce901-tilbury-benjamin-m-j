function [scores, scaling] = normalisedScore(results, scaling)
    
%     asError = (nargin>1 && asError);
    newScaling = ~(nargin>1 && ~isempty(scaling));
    
    varNames = ["mppsTotal","psnr","msssim","fade","niqe","colDiff","mic","speedqa","TError","AError"];

    seqs = string(fieldnames(results))';
    
    names = string(fieldnames(results.(seqs(1))))';
    names(strcmp(names,'props')) = [];
    
    numDehazers = length(names);
    numPoints = 0;
    numSeqs = length(seqs);
    
    for seq = seqs
        n = results.(seq).props.evaluatedFrames(2) + 1;
        numPoints = numPoints + n*numDehazers;
    end

    spmps = zeros(1,numPoints); % Will be converted to MPPS later
    psnrs = zeros(1,numPoints);
    msssims = zeros(1,numPoints);
    fades = zeros(1,numPoints);
    niqes = zeros(1,numPoints);
    colDiffs = zeros(1,numPoints);
    mics = zeros(1,numPoints);
    speedqas = zeros(1,numPoints);
    TErrors = zeros(1,numPoints);
    AErrors = zeros(1,numPoints);

    start = 1;
    for seq = string(fieldnames(results))'
        n = results.(seq).props.evaluatedFrames(2) + 1;

        for name = names
            stop = start+n-1;
            
            mic = repelem(results.(seq).(name).mic,n);
            speedqa = [results.(seq).(name).speedqa.all, mean(results.(seq).(name).speedqa.mean)];

            spmps(start:stop) = results.(seq).(name).timeTotal.all./results.(seq).props.megapixels;
            psnrs(start:stop) = results.(seq).(name).psnr.all;
            msssims(start:stop) = results.(seq).(name).msssim.all;
            fades(start:stop) = results.(seq).(name).fade.all;
            niqes(start:stop) = results.(seq).(name).niqe.all;
            colDiffs(start:stop) = results.(seq).(name).colDiff.all;
            mics(start:stop) = mic;
            speedqas(start:stop) = speedqa;

            if isnan(results.(seq).(name).TError.mean)
                TErrors(start:stop) = repelem(NaN,n);
            else
                TErrors(start:stop) = results.(seq).(name).TError.all;
            end

            if isnan(results.(seq).(name).AError.mean)
                AErrors(start:stop) = repelem(NaN,n);
            else
                AErrors(start:stop) = results.(seq).(name).AError.all;
            end

            start = start + n;
        end

    end


    % MPPS done seperately.
    full = [psnrs; msssims; fades; niqes; colDiffs; mics; speedqas; TErrors; AErrors]';
    flipScoreError = [1, 1, -1, -1, -1, 1, -1, -1, -1];
    full = full.*flipScoreError;
    
%     if nargin>1 && asError, full = -full; end
    
    %% An initial scaling is required to ensure that the combined forms are equally weighted
    % A second scaling is also required to ensure the combined values are
    % within the proper scales.
    if newScaling
        scaling.initial = cat(1,min(full,[],1),max(full,[],1));
    end
    
%     if asError
%         big = scaling.initial(1,:);
%         small = scaling.initial(2,:);
%     else
        big = scaling.initial(2,:);
        small = scaling.initial(1,:);
%     end
    
    full = (full-small)./(big-small);
    
    VIs = mean(full(:,1:3),2); % Mean of normalised psnr, msssim and fade
    NQs = mean(full(:,4:5),2); % Mean of normalised niqe and colour difference
    VCs = mean(full(:,6:7),2); % Mean of normalised mic and speedqa
    MCs = mean(full(:,8:9),2,'omitnan'); % Mean of normalised TError and AError

    % All wi.
    full = [spmps; VIs'; NQs'; VCs'; MCs']';

    %% Start process of reducing down to singular values for each dehazer
    mini = zeros(numSeqs*numDehazers,5);
    start = 1;
    m = 1;
    for seq = string(fieldnames(results))'
        n = results.(seq).props.evaluatedFrames(2) + 1;
        for name = names
            stop = start+n-1;
            mini(m,:) = mean(full(start:stop,:),1);
            start = start+n;
            m = m + 1;
        end
    end

    compact = zeros(numDehazers,5);

    for i = 1:numDehazers
        compact(i,:) = mean(mini(i:numDehazers:end,:),1);
    end

    %% Turn the SPMP into MPPS and renormalise everything
    compact(:,1) = 1./compact(:,1);
    if newScaling
        scaling.final = cat(1,min(compact,[],1),max(compact,[],1));
    end
    
%     if asError
%         big = scaling.final(1,:);
%         small = scaling.final(2,:);
%     else
        big = scaling.final(2,:);
        small = scaling.final(1,:);
%     end
    
    scores = (compact-small)./(big-small);

    
end
