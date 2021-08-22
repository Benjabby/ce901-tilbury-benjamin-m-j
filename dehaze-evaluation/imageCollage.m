function imageCollage(dehazers,path,frameNum)
   
    framePath = num2str(frameNum,'%010.f')+".png";
    J = imread(fullfile(path,"image_02","data",framePath));
    J = im2double(J);
    I = imread(fullfile(path,"haze","image",framePath));
    I = im2double(I);
    
    framePath = num2str(frameNum+1,'%010.f')+".png";
    IPlus = imread(fullfile(path,"haze","image",framePath));
    IPlus = im2double(IPlus);
    
    rcalib = loadCalibrationCamToCam(fullfile(path,"..","calib_cam_to_cam.txt"));
    knowns.K = rcalib.P_rect{3};
    
    imgs = zeros([size(J) length(dehazers)+2]);
    imgs(:,:,:,1) = J;
    imgs(:,:,:,2) = I;
    
    i = 3;
    for dehazer = dehazers
        dehazer.newSequence(knowns);
        if dehazer.FrameDelay
            dehazer.dehazeFrame(I);
            imgs(:,:,:,i) = dehazer.dehazeFrame(IPlus);
        else
            imgs(:,:,:,i) = dehazer.dehazeFrame(I);
        end
        
        imwrite(imgs(:,:,:,i),append("T:\ce901_tilbury_benjamin_m_j\machineA-local-only\dehaze-evaluation\tests\figs\",dehazer.Name,".png"));
        
        i = i + 1;
    end
    
    figure;montage(imgs,'Size', [5, 2],'ThumbnailSize',[200,Inf]);
end

