function [timings, sizes] = empiricalComplexity(dehazers, path, seed)
    
    if ~exist('seed','var')
        seed = 0;
    end
    
    rng(seed);
    
    %% Prepare testing data
    dateFolders = dir(path);
    dateFolders = dateFolders([dateFolders.isdir]);
    numSeqs = 0;
    
%     fid = fopen(fullfile(fileparts(mfilename('fullpath')),"..","haze-video-dataset","haze-subset.txt"));
%     S = textscan(fid,"%s");
%     sequences = string(S{1});
%     fclose(fid);
    
    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        if ~isfile(fullfile(datePath,'calib_cam_to_cam.txt')), continue; end

        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        numSeqs = numSeqs + length(sequenceFolders) - 2;

    end
    
    
    % I image will be done at these scales (in terms of total pixels h*w) 8x, 4x, 2x, 1x, 0.5x, 0.25x, 0.125x, 0.0625x
    % plus each of these will be done at normal aspect and a cropped square aspect.
    % So 16*numSeq readings in total. Last dimension contains the
    sizes   = zeros(1,16*numSeqs); 
    timings = zeros(length(dehazers),16*numSeqs);
    
    n = 1;
    for i = 3:length(dateFolders)
        dateName = dateFolders(i).name;
        datePath = fullfile(path,dateName);
        if ~isfile(fullfile(datePath,'calib_cam_to_cam.txt')), continue; end

        sequenceFolders = dir(datePath);
        sequenceFolders = sequenceFolders([sequenceFolders.isdir]);
        
        calib = loadCalibrationCamToCam(fullfile(datePath,'calib_cam_to_cam.txt'));
        K = calib.P_rect{3};
        knowns.K = K;
        
        for system = dehazers
           system.newSequence(knowns); 
        end

        for j = 3:length(sequenceFolders)
            sequenceName = sequenceFolders(j).name;
            sequencePath = fullfile(datePath, sequenceName);
            
            fid = fopen(fullfile(sequencePath,'haze','props.json'));
            if fid<0
                fprintf("props.json not found for sequence %s. Skipping\n",sequenceName)
                continue;
            end
            
            IPath = fullfile(sequencePath,'haze','image');
            ISize = length(dir(IPath))-2;
            
            % Choose 1 frame randomly from each sequence.
            % Don't include first or last frame.
            idx = randi(ISize-1);
            
            I = imread(fullfile(IPath,num2str(idx,'%010.f')+".png"));
            I = im2double(I);
            % Start at 4x pixel size;
            I = imresize(I,2*sqrt(2));
            ICWin = centerCropWindow2d(size(I,[1 2]),[min(size(I,[1 2])) size(I,[1 2])]);
            IC = imcrop(I, ICWin);
            
            for szf = 1:8
                sizes(n)    = size(I,1)*size(I,2);
                sizes(n+1)  = size(IC,1)*size(IC,2); 
                d = 1;
                for system = dehazers
                    
                    %% Get time on full size
                    if system.FrameDelay
                        system.dehazeFrame(I); % Dummy to do first set up. Wont be counted toward time
                        timeTic = tic;
                        system.dehazeFrame(I);
                        time = toc(timeTic);
                        % Do a new sequence to reset
                        system.newSequence(knowns); 
                    else
                        timeTic = tic;
                        system.dehazeFrame(I);
                        time = toc(timeTic);
                    end
                    
                    timings(d,n) = time;
                    
                    %% Get time on cropped
                    if system.FrameDelay
                        system.dehazeFrame(IC); % Dummy to do first set up. Wont be counted toward time
                        timeTic = tic;
                        system.dehazeFrame(IC);
                        time = toc(timeTic);
                        % Do a new sequence to reset
                        system.newSequence(knowns); 
                    else
                        timeTic = tic;
                        system.dehazeFrame(IC);
                        time = toc(timeTic);
                    end
                    
                    timings(d,n+1) = time;
                    %% Increment
                    d = d+1;
                end
                    
                n = n + 2;
                
                if szf<8
                    I = imresize(I,1/sqrt(2));
                    ICWin = centerCropWindow2d(size(I,[1 2]),[min(size(I,[1 2])) size(I,[1 2])]);
                    IC = imcrop(I, ICWin);
                end
            end
            
        end
    end
    
end

