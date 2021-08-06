function labelScatter(table, metricX, metricY,type)
    if ~exist('type','var')
        type='new';
    end
    if type=="new"
        figure;
    elseif type=="append"
        hold on;
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
    
    names = string(table.Row)';
    
    cols = fullCols(1:length(names),:);
    
    scatter([table.(metricX)],[table.(metricY)],46,cols,'filled');
    labelpoints([table.(metricX)],[table.(metricY)],names,'N',0.15,1);
    
    if type=="append"
        hold off
    end
end

