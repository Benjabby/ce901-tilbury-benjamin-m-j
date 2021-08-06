function results = resetCols(results)
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
    
    results.colCode = struct();
    results.colCode.index = 1;
    for name = results.names
        results.colCode.(name) = fullCols(results.colCode.index,:);
        results.colCode.index = results.colCode.index + 1;
    end
end

