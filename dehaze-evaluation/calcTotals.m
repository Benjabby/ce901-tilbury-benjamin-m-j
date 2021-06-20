function totals = calcTotals(results)
    totals = struct();
    
    for name = results.names
        for metric = results.metrics
            totals.(name).(metric) = 0;
        end
    end
    
    totalSample = 0;
    
    for seq = string(fieldnames(results))'
        sample = results.(seq).frames;
        totalSample = totalSample + sample;
        for name = results.names
            for metric = results.metrics
                totals.(name).(metric) = totals.(name).(metric) + results.(seq).(name).(metric) * sample;
            end
        end
    end
    
    for name = results.names
        for metric = results.metrics
            totals.(name).(metric) = totals.(name).(metric)/totalSample;
        end
    end
    
end