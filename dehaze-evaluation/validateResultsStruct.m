function valid = validateResultsStruct(results)
%     for seq = string(fieldnames(results))'
%         if seq == "names" || seq == "metrics", continue; end
%         
%         sysNames = string(fieldnames(results.(seq)))';
%         sysNames = erase(sysNames, 'visibility');
%         sysNames = erase(sysNames, 'megapixels');
%         sysNames = erase(sysNames, 'frames');
%         
%         missingNames = setdiff
%     end
end

