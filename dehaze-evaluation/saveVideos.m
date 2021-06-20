function saveVideos(sequencePath, outPath, calib, systems)
    if ~exist('systems','var')
        systems = defaultSystems();
    end
    
    for s = 1:length(systems)
        name = systems(s).name;
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
        textwaitba1r(frameID,ISize,"Creating Videos")
        for s = 1:length(systems)
                name = systems(s).name;
                if frameID==0
                    systems(s).state = struct();
                end
                [pred, ~, ~, ~, ~, systems(s).state] = systems(s).function(I, knowns, systems(s).state);
                pred = min(max(pred,0),1);
                writeVideo(writers.(name),pred);
        end
        
    end
    textwaitbar(ISize,ISize,"Creating Videos")
    
    for s = 1:length(systems)
        name = systems(s).name;
        close(writers.(name));
    end
end

