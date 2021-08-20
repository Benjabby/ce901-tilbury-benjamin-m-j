function [means, meanStrings, stdStrings] = groupResults(results)

    groupTotals = struct();
    groupTotals.dense = struct();
    groupTotals.moderate = struct();
    groupTotals.light = struct();
    groupTotals.all = struct();
    
    seqs = string(fieldnames(results))';
    
    names = string(fieldnames(results.(seqs(1))))';
    names(strcmp(names,'props')) = [];
    
    metrics = string(fieldnames(results.(seqs(1)).(names(1))))';
    
%     numDehazers = length(names);
% %     numPoints = 0;
%     numSeqs = length(seqs);
%     
    
%     num.dense = 0;
%     num.moderate = 0;
%     num.light = 0;
    
%     for seq = seqs
%         n = results.(seq).props.evaluatedFrames(2) + 1;
% %         numPoints = numPoints + n*numDehazers;
%         group = results.(seq).props.group;
%         num.(group) = num.(group) + n;
%     end
%     
%     for name = names
%         for g = ["dense", "moderate", "light"]
%             groupTotals
%         end
%     end
    
    for seq = seqs
        group   = results.(seq).props.group;
        mp      = results.(seq).props.megapixels;
        for name = names
            if ~isfield(groupTotals.all, name), groupTotals.all.(name) = struct(); end
            if ~isfield(groupTotals.(group), name), groupTotals.(group).(name) = struct(); end
            for metric = metrics
                if ~isfield(results.(seq).(name),metric)
                    val = NaN;
                elseif isstruct(results.(seq).(name).(metric))
                    if ~isnan(results.(seq).(name).(metric).mean)
                        val = results.(seq).(name).(metric).all;
                    else
                        val = NaN;
                    end
                else
                    val = results.(seq).(name).(metric);
                end
                
                % Seconds to seconds-per-megapixel
                if contains(metric,"time"), val = val./mp; end
                
                if isfield(groupTotals.(group).(name), metric)
                    groupTotals.(group).(name).(metric) = [groupTotals.(group).(name).(metric) val];
                else
                    groupTotals.(group).(name).(metric) = val;
                end
                
                if isfield(groupTotals.all.(name), metric)
                    groupTotals.all.(name).(metric) = [groupTotals.all.(name).(metric) val];
                else
                    groupTotals.all.(name).(metric) = val;
                end
            end
        end
    end

    metrics = ["mppsTotal","psnr","msssim","fade","colDiff","niqe","mic","speedqa","TError","AError"];

%     strtables = struct();
    meanStrings = struct();
    stdStrings = struct();
    means = struct();

    for group = string(fieldnames(groupTotals))'
        means.(group) = table('Size', [length(names) length(metrics)], 'VariableTypes',repelem("double",length(metrics)),'VariableNames',metrics,'RowNames',names);
        meanStrings.(group) = table('Size', [length(names) length(metrics)], 'VariableTypes',repelem("string",length(metrics)),'VariableNames',metrics,'RowNames',names);
        stdStrings.(group) = table('Size', [length(names) length(metrics)], 'VariableTypes',repelem("string",length(metrics)),'VariableNames',metrics,'RowNames',names);
        for name = names
            for metric = metrics
                
                if metric=="mppsTotal"
                    meanstr = string(round(1./mean(groupTotals.(group).(name).timeTotal,'omitnan'),4,'significant'));
                    val = 1./mean(groupTotals.(group).(name).timeTotal,'omitnan');

                    % Okay so... 
                    % Initially I was going to take the mean of the
                    % MPPS, as in mean((1./timeTotals)). [in fact this was originaly calculated as part of the evaluator so there was an additional mppsTotal metric that was just megapixels/timeTotal]
                    % However that doesn't make sense.
                    % For example if two things took 0.1 seconds and 1.2
                    % seconds, the mean time is 0.65
                    % the first thing would have a framerate of 10 and
                    % the second of 0.8333. 
                    % the mean of this is 5.4167 FPS.
                    % The FPS of the mean is 1.5385 which makes more
                    % sense because this is done frame-by-frame so the
                    % FPS should be calculated on an averaged time not
                    % averaged on calculated FPS. 
                    % ANYWAY. Because of that, there's no easy way to
                    % get the standard deviation of the MPPS.
                    stdstr = "-";
                elseif ~isfield(groupTotals.(group).(name),metric) || all(isnan(groupTotals.(group).(name).(metric)))
%                     str = "N/A";
                    val = NaN;
                    meanstr = "N/A";% NaN;
                    stdstr = "N/A";% NaN;
                else
                    meanstr = string(round(mean(groupTotals.(group).(name).(metric),'omitnan'),4,'significant'));
                    val = mean(groupTotals.(group).(name).(metric),'omitnan');
                    stdstr = string(round(std(groupTotals.(group).(name).(metric),'omitnan'),4,'significant'));
%                     stdev = std(groupTotals.(group).(sys).(metric),'omitnan');
%                     str = append(meanstr," ",char(177)," ",stdstr);

                end

                meanStrings.(group)(name,metric) = {meanstr};
                stdStrings.(group)(name,metric) = {stdstr};
                means.(group)(name,metric) = {val};
%                 strtables.(group)(sys,metric) = {str};
            end
        end
    end
end