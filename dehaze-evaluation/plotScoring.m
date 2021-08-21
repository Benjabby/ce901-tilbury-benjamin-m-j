function plotScoring(results, comparings, excludeComparings, firstNormOnly)
    % Comparings are the names of methods excluded from the scale calculation.
    % In my case this will be my methods as they are to be compared against
    % the scales of the already evaluated methods.
    
    seqs = string(fieldnames(results))';
    names = string(fieldnames(results.(seqs(1))))';
    names(strcmp(names,'props')) = [];
    
    firstOnly = (nargin>3 && ~isempty(firstNormOnly) && firstNormOnly);
    
    if nargin>1 && ~isempty(comparings)
        partA = struct();
        partB = struct();
        
        partANames = setdiff(names,comparings);
        for seq = seqs
            partA.(seq) = results.(seq);
            partB.(seq) = results.(seq);
            
            for pa = partANames
                partB.(seq) = rmfield(partB.(seq),pa);
            end
            
            for pb = comparings
                partA.(seq) = rmfield(partA.(seq),pb);
            end
        end
        
        % Get the order right
        if nargin>2 && excludeComparings
            names = partANames;
            igrey = 0;
        else
            names = [partANames, comparings];
            igrey = length(partANames);
        end
    else
        partA = results;
        partB = [];
        igrey = 0;
    end
    
    [scores, scaling] = normalisedScore(partA,[],firstOnly);
    
    
    fullCols = [
            31 119 180;
            214 39 40;
            44 160 44;
            255 127 14;
            148 103 189;
            140 86 75;
            227 119 194;
%             127 127 127;
            23 190 207
            188 189 34;
            0 0 0
            0 0 0
            0 0 0 ] / 255;
    
    cols = zeros(length(names),5,3);
    for i = 1:length(names)
        k = i-igrey;
        if k<1
            cols(i,:,:) = repmat(fullCols(i,:).*0.5+[0.5,0.5,0.5],5,1);
        else
            cols(i,:,:) = repmat(fullCols(k,:),5,1);
        end
    end
    
    if ~isempty(partB) && (nargin<=2 || ~excludeComparings)
        scoresB = normalisedScore(partB,scaling,firstOnly);
        scores = cat(1,scores,scoresB);
    end

    bars=bar(scores','FaceColor','flat');

    
    for i=1:length(bars)
       bars(i).CData = squeeze(cols(i,:,:));
    end
    set(gca,'XTickLabel',["Realtime Capability","Visibility Improvement","Naturalness Quality","Video Consistency","Model Conformity"]);
    legend(names);
end
