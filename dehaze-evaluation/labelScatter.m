function labelScatter(table, metricX, metricY, ignores, type)
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
    
    names = string(table.Row);
    
    cols = fullCols(1:length(names),:);
    
%     if ~exist('ignores','var') && ~isempty(ignores)
    
%     gscatter([table.(metricX)],[table.(metricY)],names,cols,'.',52);
    scatter([table.(metricX)],[table.(metricY)],120,cols,'filled');
    labelpoints([table.(metricX)],[table.(metricY)],names','N',0.15,1);
%     legend(names);
    
    if type=="append"
        hold off
    end
end

