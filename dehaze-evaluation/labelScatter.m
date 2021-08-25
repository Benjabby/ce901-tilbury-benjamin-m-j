function labelScatter(table, metricX, metricY, ignore, type)
    if ~exist('type','var')
        type='new';
    end
    if type=="new"
        figure;
    elseif type=="append"
        hold on;
    end
    
    if nargin>3 && ~isempty(ignore)
        names = setdiff(string(table.Row),ignore');
        table = table(names,:);
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
    
    if ~(nargin>3 && ~isempty(ignore))
        existingCols = existingCols.*0.5 + 0.5;
    end
    
    myCols = [ 31 119 180;
            214 39 40;
            44 160 44; ] / 255;
    
    fullCols = cat(1,existingCols,myCols);
    
    names = table.Row;
%     names = append(string(table.Row)," et al.");
    
    cols = fullCols(1:length(names),:);
    
    
%     gscatter([table.(metricX)],[table.(metricY)],names,cols,'.',52);
    scatter([table.(metricX)],[table.(metricY)],120,cols,'filled');
    labelpoints([table.(metricX)],[table.(metricY)],names','N',0.15,1);
%     legend(names);
    
    if type=="append"
        hold off
    end
end

