function plotVisibility(results,metric)
    N = length(fieldnames(results))-3;
    x = zeros(1,N);
    vals = zeros(length(results.names),N);
    cols = zeros(length(results.names),3);
    
    i = 1;
    for seq = string(fieldnames(results))'
        if seq == "names" || seq == "metrics" || seq == "colCode", continue; end
        x(i) = results.(seq).visibility * 1000; % Visibility in meters
        for s = 1:length(results.names)
            name = results.names(s);
            vals(s,i) = results.(seq).(name).(metric);
            cols(s,:) = results.colCode.(name);
        end
        i = i + 1;
    end
    
    [x, sortidx] = sort(x);
    assignin('base','vis',x);
    vals = vals(:,sortidx);
    
    h = plot(x,vals);
    set(h, {'color'}, num2cell(cols, 2));
    
end

