function saveVideos(sequencePath, outPath, calib, systems)
    if ~exist('systems','var')
        systems = defaultSystems();
    elseif length(fieldnames(systems))<2 || ~isfield(systems,'prepped') || ~systems.prepped
        isHandles = true;
        names = fieldnames(systems);
        for s = 1:length(names)
            name = names{s};
           isHandles = isHandles & isa(systems.(name),'function_handle');
        end
        
        if isHandles
            warning("'systems' has not been prepared properly, attempting to prepare now")
            systems = prepareSystems(systems);
        else
            error("Incorrectly formatted 'systems' structure")
        end
    end
    
    for name = string(fieldnames(systems))'
        if name == "prepped", continue; end
        writers.(name) = VideoWriter(fullfile(outPath,name+".mp4"),"MPEG-4");
        writers.(name).FrameRate = 10;
        open(writers.(name));
    end
    
    knowns.K = calib.P_rect{3};
    IPath = fullfile(sequencePath,'haze','image');
    ISize = length(dir(IPath))-2;
    

    for frameID = 0:ISize-1
        framePath = num2str(frameID,'%010.f')+".png";
        I = imread(fullfile(IPath,framePath));
        I = im2double(I);
        %Progress_bar(ISize);
        textwaitbar(frameID,ISize,"Creating Videos")
        for name = string(fieldnames(systems))'
            if name == "prepped", continue; end
                    
            if frameID==0
                systems.(name).state = struct();
            end
            [pred, ~, ~, ~, ~, systems.(name).state] = systems.(name).function(I, knowns, systems.(name).state);
            if pred
                pred = min(max(pred,0),1);
                writeVideo(writers.(name),pred);
            end
        end
        
    end
    
    for name = string(fieldnames(systems))'
        if name == "prepped", continue; end
                    
        if systems.(name).frameDelay>0
            for fd = 0:systems.(name).frameDelay
                [pred, ~, ~, ~, ~, systems.(name).state] = systems.(name).function(I, knowns, systems.(name).state);
                pred = min(max(pred,0),1);
                writeVideo(writers.(name),pred);
            end
        end
        close(writers.(name));
    end
    
    textwaitbar(ISize,ISize,"Creating Videos")
end

