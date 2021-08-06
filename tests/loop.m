function [lines] = loop(path) 
    dateFolders = dir(path);
    dateFolders = dateFolders([dateFolders.isdir]);
    
    
    remainingSequences = 0;
    remainingFrames = 0;
    
    %% Preliminary look to find total number of frames
    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        
        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        
        for j = 3:length(sequenceFolders)
            sequenceName = sequenceFolders(j).name;
            sequencePath = fullfile(datePath, sequenceName);
            
            IPath = fullfile(sequencePath,'haze','image');
            ISize = length(dir(IPath))-2;
            
            remainingFrames = remainingFrames+ISize;
        end
    end
    
    lines = struct();

    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        
        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        for j = 3:length(sequenceFolders)
            sequenceName = sequenceFolders(j).name;
            sequencePath = fullfile(datePath, sequenceName);
            
            
            IPath = fullfile(sequencePath,'haze','image');
            ISize = length(dir(IPath))-2;
            
            eh = ['S' sequenceName];
            lines.(eh) = zeros(1,ISize);
            
            fid = fopen(fullfile(sequencePath,'haze','props.json'));
            if fid<0
                fprintf("props.json not found for sequence %s. Skipping\n",sequenceName)
                continue;
            end
            rawJSON = fread(fid);
            str = char(rawJSON');
            fclose(fid);
            hazeProps = jsondecode(str);
            
            vis = hazeProps.visibility;
            vis = 3.912/vis; % turn into beta
            
            prevVH = 175;
            
            for frameID = 0:ISize-1
                framePath = num2str(frameID,'%010.f')+".png";
                I = imread(fullfile(IPath,framePath));
                I = im2double(I);
                
                darkChannel = min(I,[],3);
                
                bigSE = strel('square',31);
                h = (imdilate(darkChannel,bigSE)-imerode(darkChannel,bigSE))==0;
                VH = find(max(h,[],2), 1, 'last');
                if isempty(VH)
                    VH = prevVH;
                end
                lines.(eh)(frameID+1)=VH;
                
                prevVH = VH;
            end
        end
    end