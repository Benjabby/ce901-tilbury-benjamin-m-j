function plotScoring(results, comparings)
    % Comparings are the names of methods excluded from the scale calculation.
    % In my case this will be my methods as they are to be compared against
    % the scales of the already evaluated methods.
    
    seqs = string(fieldnames(results))';
    names = string(fieldnames(results.(seqs(1))))';
    names(strcmp(names,'props')) = [];
    
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
        names = [partANames, comparings];
    else
        partA = results;
        partB = [];
    end
    
    [scores, scaling] = normalisedScore(partA);
    
    if ~isempty(partB)
        scoresB = normalisedScore(partB,scaling);
        scores = cat(1,scores,scoresB);
    end

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

    bars=bar(scores','FaceColor','flat');

    for i=1:length(bars)
       bars(i).CData = repmat(fullCols(i,:),5,1);
    end
    set(gca,'XTickLabel',["Realtime Performance","Visibility Improvement","Naturalness Quality","Video Consistency","Model Conformity"]);
    legend(names);
end
