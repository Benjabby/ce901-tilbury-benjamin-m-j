function tables = groupResults(results)
    groupTotals = struct();

    for seq = string(fieldnames(results))'
        group = results.(seq).props.group;
        if ~isfield(groupTotals, group), groupTotals.(group) = struct(); end
        for sys = string(fieldnames(results.(seq)))'
            if sys == "props", continue; end

            if ~isfield(groupTotals.(group), sys), groupTotals.(group).(sys) = struct(); end

            for metric = string(fieldnames(results.(seq).(sys)))'

                if isstruct(results.(seq).(sys).(metric))
                    if isfield(results.(seq).(sys).(metric), 'all')
                        val = results.(seq).(sys).(metric).all;
                    else
                        val = NaN;
                    end
                else
                    val = results.(seq).(sys).(metric);
                end

                if isfield(groupTotals.(group).(sys), metric)
                    groupTotals.(group).(sys).(metric) = [groupTotals.(group).(sys).(metric) val];
                else
                    groupTotals.(group).(sys).(metric) = val;
                end
            end
        end
    end

    names = string(fieldnames(groupTotals.dense))';
    metrics = string(fieldnames(groupTotals.dense.(names(1))))';

    strtables = struct();
    tables = struct();

    for group = string(fieldnames(groupTotals))'
        strtables.(group) = table('Size', [length(names) length(metrics)], 'VariableTypes',repelem("string",length(metrics)),'VariableNames',metrics,'RowNames',names);
        tables.(group) = table('Size', [length(names) length(metrics)], 'VariableTypes',repelem("double",length(metrics)),'VariableNames',metrics,'RowNames',names);
        for sys = names
            for metric = metrics

                if ~isfield(groupTotals.(group).(sys),metric) || all(isnan(groupTotals.(group).(sys).(metric)))
                    str = "N/A";
                    val = NaN;
                else
                    meanstr = string(round(mean(groupTotals.(group).(sys).(metric),'omitnan'),4,'significant'));
                    val = mean(groupTotals.(group).(sys).(metric),'omitnan');
                    stdstr = string(round(std(groupTotals.(group).(sys).(metric),'omitnan'),4,'significant'));

                    str = append(meanstr," ",char(177)," ",stdstr);
                end

                tables.(group)(sys,metric) = {val};
                strtables.(group)(sys,metric) = {str};
            end
        end
    end
end