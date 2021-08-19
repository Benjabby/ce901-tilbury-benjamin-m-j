function labelScatter(table, metricX, metricY,type)
    if ~exist('type','var')
        type='new';
    end
    if type=="new"
        figure;
    elseif type=="append"
        hold on;
    end
    
    existingCols = [
            31 119 180;
            214 39 40;
            44 160 44;
            255 127 14;
            148 103 189;
            140 86 75;
            227 119 194;
            23 190 207
             ] / 255;
     
    existingCols = existingCols.*0.5 + 0.5;
    myCols = [ 31 119 180;
            214 39 40;
            44 160 44; ] / 255;
    
    fullCols = cat(1,existingCols,myCols);
    
    names = string(table.Row)';
    
    cols = fullCols(1:length(names),:);
    
    scatter([table.(metricX)],[table.(metricY)],64,cols,'filled');
%     labelpoints([table.(metricX)],[table.(metricY)],names,'NE',0.15,1);
    
    if type=="append"
        hold off
    end
end

