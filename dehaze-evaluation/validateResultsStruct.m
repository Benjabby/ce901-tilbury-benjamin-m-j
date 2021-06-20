function valid = validateResultsStruct(results)
    for seq = string(fieldnames(results))'
        if seq == "names" || seq == "metrics", continue; end
        
    end
end

