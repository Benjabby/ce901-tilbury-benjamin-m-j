function results = mergeResults(original, new, replace)
    if ~exist('replace','var'), replace = false; end
    
%     if replace
%         newMetrics  = new.metrics;
%         newNames    = new.names;
%         newSeqs     = string(fieldnames(new))';
%         newSeqs(strcmp(newSeqs, 'names')) = [];
%         newSeqs(strcmp(newSeqs, 'metrics')) = [];
%     else
%         newMetrics  = setdiff(new.metrics, original.metrics);
%         if isempty(newMetrics)
%             newNames    = setdiff(new.names, original.names);
%         else
%             newNames = union(original,
%         end
%         
%         if isempty(newNames) && isempty(newMetrics)
%             newSeqs     = setdiff(string(fieldnames(new))',string(fieldnames(original))');
%         else
%             newSeqs     = union(string(fieldnames(new))',string(fieldnames(original))');
%             newSeqs(strcmp(newSeqs, 'names')) = [];
%             newSeqs(strcmp(newSeqs, 'metrics')) = [];
%         end
%     end
    
    for seq = string(fieldnames(new))'
        if seq=="names" || seq=="metrics" || seq == "colCode", continue; end
        original.(seq).frames = new.(seq).frames;
        original.(seq).visibility = new.(seq).visibility;
        original.(seq).megapixels = new.(seq).megapixels;
        for name = new.names
            for metric = new.metrics
               if replace || ~(isfield(original,seq) && isfield(original.(seq),name) && isfield(original.(seq).(name),metric))
                   original.(seq).(name).(metric) = new.(seq).(name).(metric);
               end
            end
        end
    end
    original.names      = union(original.names, new.names);
    original.metrics    = union(original.metrics, new.metrics);
    
    results = original;
end