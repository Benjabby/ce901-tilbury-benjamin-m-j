function saveVideos(systems, path, outPath)

    saveIJ = (nargin==0 || isempty(systems));
        
    if ~exist('path','var')
        path = "haze-video-dataset\dataset";
    end
    
    if ~exist('outPath','var')
        outPath = "tests\videos";
    end

    dateFolders = dir(path);
    dateFolders = dateFolders([dateFolders.isdir]);

    remainingSequences = 0;
    
    fid = fopen(fullfile(fileparts(mfilename('fullpath')),"..","haze-video-dataset","haze-subset.txt"));
    S = textscan(fid,"%s");
    sequences = string(S{1});
    fclose(fid);

    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        if ~isfile(fullfile(datePath,'calib_cam_to_cam.txt')), continue; end

        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        for j = 3:length(sequenceFolders)
            sequenceName = sequenceFolders(j).name;
            if ~any(strcmp(sequences,sequenceName)), continue;  end
            remainingSequences = remainingSequences + 1;
        end

    end

    names    = strings(1,remainingSequences);
    paths    = strings(1,remainingSequences);
    As      =  zeros(3,remainingSequences);
    vises      =  zeros(1,remainingSequences);
    knowns(1,remainingSequences)     = struct('K',{zeros(3,4)});
    numFrames = zeros(1,remainingSequences);
    n = 1;

    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        if ~isfile(fullfile(datePath,'calib_cam_to_cam.txt')), continue; end

        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);

        calib = loadCalibrationCamToCam(fullfile(datePath,'calib_cam_to_cam.txt'));
        K = calib.P_rect{3};
        kn.K = K;

        for j = 3:length(sequenceFolders)
            sequenceName = sequenceFolders(j).name;
            sequencePath = fullfile(datePath, sequenceName);
            
            if ~any(strcmp(sequences,sequenceName)), continue;  end

            fid = fopen(fullfile(sequencePath,'haze','props.json'));
            if fid<0
                fprintf("props.json not found for sequence %s. Skipping\n",sequenceName)
                continue;
            end
            rawJSON = fread(fid);
            str = char(rawJSON');
            fclose(fid);
            hazeProps = jsondecode(str);

            JPath = fullfile(sequencePath,'image_02','data');

            JSize = length(dir(JPath))-2;

            numFrames(n) = JSize;
            names(n) = ['S' sequenceName(1:end-5)];
            paths(n) = sequencePath;
            knowns(n) = kn;
            As(:,n) = hazeProps.A;
            vises(n) = hazeProps.visibility;

            n = n + 1;

        end
    end

    seqTab = table(paths',As',vises',knowns',numFrames','RowNames',names,'VariableNames',["path","A","vis","knowns","frames"]);
    seqTab = sortrows(seqTab,'vis','descend');

    remainingFrames = sum(seqTab.frames);

    if saveIJ
        IWriter = VideoWriter(fullfile(outPath,"I.mp4"),"MPEG-4");
        IWriter.FrameRate = 20; % To be able to watch them all faster. 
        open(IWriter);
        JWriter = VideoWriter(fullfile(outPath,"J.mp4"),"MPEG-4");
        JWriter.FrameRate = 20; % To be able to watch them all faster. 
        open(JWriter);
        
        for irow = 1:height(seqTab)
            sequencePath = seqTab(irow,:).path;
            IPath = fullfile(sequencePath,'haze','image');
            JPath = fullfile(sequencePath,'image_02','data');
            ISize = length(dir(IPath))-2;
            
            for frameID = 0:ISize-1
                framePath = num2str(frameID,'%010.f')+".png";
                I = imread(fullfile(IPath,framePath));
                I = im2double(I);
                I = imresize(I, [375, 1242],'bilinear');
                writeVideo(IWriter,I);

                J = imread(fullfile(JPath,framePath));
                J = im2double(J);
                J = imresize(J, [375, 1242],'bilinear');
                writeVideo(JWriter,J);

                textwaitbar(frameID,ISize,"Creating Videos")

            end
        end
        
        close(IWriter);
        close(JWriter);
    else
        for system = systems
            name = system.Name;
            writer = VideoWriter(fullfile(outPath,name+".mp4"),"MPEG-4");
            writer.FrameRate = 20; % To be able to watch them all faster. 
            open(writer);
            i = 1;

            for irow = 1:height(seqTab)
                sequencePath = seqTab(irow,:).path;
                IPath = fullfile(sequencePath,'haze','image');

                system.newSequence(seqTab(irow,:).knowns);
                for frameID = 0:seqTab(irow,:).frames-1
                    framePath = num2str(frameID,'%010.f')+".png";
                    I = imread(fullfile(IPath,framePath));
                    I = im2double(I);
                    pred = system.dehazeFrame(I);
                    if frameID-system.FrameDelay>=0
                        pred = imresize(pred, [375, 1242],'bilinear');
                        writeVideo(writer,pred);
                    end

                    textwaitbar(i,remainingFrames,"Creating Videos")
                    i = i+1;
                end
            end
            if system.FrameDelay
                pred = system.dehazeFrame([]);
                pred = imresize(pred, [375, 1242],'bilinear');
                writeVideo(writer,pred);
            end

            close(writer);
        end
    end

    textwaitbar(remainingFrames,remainingFrames,"Creating Videos")
end
