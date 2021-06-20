function systems = prepareSystems(systems)
    % Input arguments:
    %   systems - A structure containing the systems to be prepared; each
    %             field must be the name of a system with its value being the
    %             function. The function must conform to the correct format. See the
    %             readme file in 'dehaze-systems' for the correct format.
    % Outputs:
    %   systems - A structure containing the systems to be evaluated,
    %             modified into a proper format include properties of the system 
    %%
    
    isHandles = true;
    names = fieldnames(systems);
    for s = 1:length(names)
        name = names{s};
       isHandles = isHandles & isa(systems.(name),'function_handle');
    end
    
    if ~isHandles
        error("Incorrectly formatted 'systems' structure. 'systems' as passed to 'prepareSystems' should contain only function handles")
    end
        
    dummyImg = ones([512,512,3]);
    dummyKnowns.K = ones(3,4);
    
    removals = {};
    
    for name = string(fieldnames(systems))'
        state = struct();
        
        sysFunc = systems.(name);
        
        try
            [predImage, predT, predA, timeImage, timeA, state] = sysFunc(dummyImg,dummyKnowns,state);
        catch err
            warning("System %s%s produced the following error. Removing from systems.\n\n%s",name,strtrim(mlreportgen.utils.toString(sysFunc)),err.msgtext)
            removals{end+1} = name;
            continue
        end
        
        frameDelay = 0;
        
        isError = false;
        
        while isempty(predImage)
            % There is a frame delay so keep feeding it until it gives output
            frameDelay = frameDelay + 1;
            if frameDelay > 1 %30
                warning("Currently only frame delays of 1 are supported in evaluation; function %s%s has a larger frame delay and will not be evaluated", name, strtrim(mlreportgen.utils.toString(sysFunc)))
%                 error("Function %s%s has output too many consequetive null values to be considered an exceptable delay. Please check the code", name, strtrim(mlreportgen.utils.toString(sysFunc)))
                isError = true;
                break
            end
            
            try
                [predImage, predT, predA, timeImage, timeA, state] = sysFunc(dummyImg,dummyKnowns,state);
            catch err
                warning("System %s%s produced the following error. Removing from systems.\n\n%s",name,strtrim(mlreportgen.utils.toString(sysFunc)),err.msgtext)
                isError = true;
                break
            end
        end
        
        if isError
            removals{end+1} = name;
            continue;
        end
        
        systems.(name) = struct();
        systems.(name).function = sysFunc;
        systems.(name).frameDelay = frameDelay;
        systems.(name).predictsA = ~isempty(predA);
        systems.(name).predictsT = ~isempty(predT);
        
        assert(~xor(isempty(predA),isempty(timeA)),"Function %s%s outputs a predA but not timeA", name, mlreportgen.utils.toString(sysFunc)) 
        assert(~isempty(timeImage),"Function %s%s fails to provide timeImage", name, sysFunc)
    end
    
    for r = string(removals)
        systems.(r) = [];
    end
    
    systems.prepped = true;
end

