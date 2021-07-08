function totals = calcTotals(results)
    totals = struct();
    i = 1;
    for name = results.names
        totals(i).name = name;
        totals(i).col = results.colCode.(name);
        for metric = results.metrics
            totals(i).(metric) = 0;
        end
        i = i+1;
    end
    
    totalSample = 0;
    
    for seq = string(fieldnames(results))'
        if seq == "names" || seq == "metrics" || seq == "colCode", continue; end
        sample = results.(seq).frames;
        totalSample = totalSample + sample;
        i = 1;
        for name = results.names
            for metric = results.metrics
                totals(i).(metric) = totals(i).(metric) + results.(seq).(name).(metric) * sample;
            end
            i = i+1;
        end
    end
    
    i = 1;
    for name = results.names
        for metric = results.metrics
            totals(i).(metric) = totals(i).(metric)/totalSample;
        end
        i = i+1;
    end
    
end